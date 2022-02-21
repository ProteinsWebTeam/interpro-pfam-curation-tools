import sys, os, json, ssl, re
from urllib import request
from urllib.error import HTTPError
from time import sleep
import argparse


def has_swissprot(BASE_URL):
    # disable SSL verification to avoid config issues
    context = ssl._create_unverified_context()

    next = BASE_URL

    attempts = 0

    while next:
        try:
            req = request.Request(next, headers={"Accept": "application/json"})
            res = request.urlopen(req, context=context)
            # If the API times out due a long running query
            if res.status == 408:
                # wait just over a minute
                sleep(61)
                # then continue this loop with the same URL
                continue
            elif res.status == 204:
                # no data so leave loop
                break
            payload = json.loads(res.read().decode())
            next = ""
            attempts = 0

        except HTTPError as e:
            if e.code == 408:
                sleep(61)
                continue
            else:
                # If there is a different HTTP error, it wil re-try 3 times before failing
                if attempts < 3:
                    attempts += 1
                    sleep(61)
                    continue
                else:
                    sys.stderr.write("LAST URL: " + next)
                    raise e

        count_swissprot = payload["proteins"]["reviewed"]
        count_trembl = payload["proteins"]["unreviewed"]

    if count_swissprot > 10:
        return True
    else:
        # print(count_trembl)
        return False


def search_trembl_pmid(boring_pmids, BASE_URL):
    # disable SSL verification to avoid config issues
    context = ssl._create_unverified_context()

    next = BASE_URL
    attempts = 0

    list_pmid_acc = set()

    while next:
        try:
            req = request.Request(next, headers={"Accept": "application/json"})
            res = request.urlopen(req, context=context)
            # If the API times out due a long running query
            if res.status == 408:
                # wait just over a minute
                sleep(61)
                # then continue this loop with the same URL
                continue
            elif res.status == 204:
                # no data so leave loop
                break
            payload = json.loads(res.read().decode())
            next = payload["next"]
            attempts = 0

        except HTTPError as e:
            if e.code == 408:
                sleep(61)
                continue
            else:
                # If there is a different HTTP error, it wil re-try 3 times before failing
                if attempts < 3:
                    attempts += 1
                    sleep(61)
                    continue
                else:
                    sys.stderr.write("LAST URL: " + next)
                    raise e

        for i, item in enumerate(payload["results"]):
            # get UniProt accession
            accession = item["metadata"]["accession"]
            # search for list of PMIDs
            list_pmid = search_pmid(boring_pmids, accession, list_pmid_acc)
            if len(list_pmid) != 0:
                for pmid in list_pmid:
                    list_pmid_acc.add(pmid)

        # Don't overload the server, give it time before asking for more
        if next:
            sleep(1)

    # print(list_pmid_acc)
    return list_pmid_acc


def load_boring_pmids(boringfile):
    print("Loading boring PMIDs into memory")
    boring_pmids = []
    with open(boringfile, "r") as f:
        for line in f:
            pmid = line.strip("\n").split()[0]
            boring_pmids.append(pmid)

    return boring_pmids


def search_pmid(boring_pmids, accession, list_pmid_acc):
    # disable SSL verification to avoid config issues
    context = ssl._create_unverified_context()
    next = f"https://www.ebi.ac.uk/proteins/api/proteins/{accession}"
    attempts = 0

    while next:
        try:
            req = request.Request(next, headers={"Accept": "application/json"})
            res = request.urlopen(req, context=context)
            # If the API times out due a long running query
            if res.status == 408:
                # wait just over a minute
                sleep(61)
                # then continue this loop with the same URL
                continue
            elif res.status == 204:
                # no data so leave loop
                break
            payload = json.loads(res.read().decode())
            next = ""
            attemps = 0
        except HTTPError as e:
            if e.code == 408:
                sleep(61)
                continue
            else:
                # If there is a different HTTP error, it wil re-try 3 times before failing
                if attempts < 3:
                    attempts += 1
                    sleep(61)
                    continue
                else:
                    sys.stderr.write("LAST URL: " + next)
                    raise e

    list_pmids = set()
    pmid = ""

    # print(accession)
    for i, item in enumerate(payload["references"]):
        if "dbReferences" in item["citation"]:
            title = item["citation"]["title"]
            t = re.search(r"(gene|genome|Gene|Genome|Genomic|genomic|sequen|Sequen)", title)

            for j, ref in enumerate(item["citation"]["dbReferences"]):
                if ref["type"] == "PubMed":
                    pmid = ref["id"]
                    # print(accession, pmid)
                    if not t and pmid not in boring_pmids and pmid not in list_pmid_acc:
                        # print(accession, pmid, title)
                        list_pmids.add(pmid)

    return list_pmids


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--database", help="resource database to search data for", required=True
    )
    parser.add_argument(
        "-i",
        "--inputf",
        metavar="METHOD_FILE",
        help="file containing PANTHER accessions",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--outputf",
        metavar="OUTPUT_FILE",
        help="file to write output results in",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--boringpmidf",
        metavar="PMID_FILE",
        help="file containing boring PMID accessions",
        required=True,
    )

    args = parser.parse_args()
    if not os.path.isfile(args.inputf):
        parser.error(f"Error file not found '{args.inputf}'")
    if not os.path.isfile(args.boringpmidf):
        parser.error(f"Error file not found '{args.boringpmidf}'")

    # load list of boring pmids into memory
    boring_pmids = load_boring_pmids(args.boringpmidf)

    # get list of panther accession from file
    print(f"Searching data for {args.database} signatures")
    # list_signatures = dict()
    with open(args.inputf, "r") as f:
        with open(args.outputf, "w") as outf:
            for line in f:
                signature = line.strip("\n")
                # print(signature)
                url = f"https://www.ebi.ac.uk:443/interpro/api/protein/entry/{args.database}/{signature}/"
                if has_swissprot(url):
                    pass
                else:
                    url = f"https://www.ebi.ac.uk:443/interpro/api/protein/unreviewed/entry/{args.database}/{signature}/?page_size=200"
                    list_pmid = search_trembl_pmid(boring_pmids, url)
                    # sys.exit()
                    if len(list_pmid) != 0:
                        text = " | ".join(list_pmid)
                        outf.write(f"{signature},{text}\n")

# get list of panther signatures unintegrated without comments:
# select m.method_ac
# from interpro.method m
# left join interpro.entry2method e2m on m.method_ac=e2m.method_ac
# left join interpro.method_comment c on c.method_ac=m.method_ac
# where m.method_ac like 'PTHR%' and m.method_ac not like 'PTHR%:%' and c.status is null and e2m.method_ac is null
# ;
