import sys, os, json, ssl, re
from urllib import request
from urllib.error import HTTPError
from time import sleep
import argparse
from configparser import ConfigParser
from multiprocessing import Pool


class memberdb_pmid:
    def __init__(self, member_db, boringfile):
        self.load_boring_pmids(boringfile)
        self.database = member_db
        self.sign_in = list()

    def load_boring_pmids(self, boringfile):
        print("Loading boring PMIDs into memory")
        self.boring_pmids = list()
        if not os.path.isfile(boringfile):
            print(f"Error file not found '{boringfile}'")
            sys.exit(1)
        with open(boringfile, "r") as f:
            for line in f:
                pmid = line.strip("\n").split()[0]
                self.boring_pmids.append(pmid)

    def has_swissprot(self, signature):
        # disable SSL verification to avoid config issues
        context = ssl._create_unverified_context()

        next = f"https://www.ebi.ac.uk:443/interpro/api/protein/entry/{self.database}/{signature}/"

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
                        print("LAST URL: " + next)
                        print(e)
                        next = ""

            count_swissprot = 0
            if "reviewed" in payload["proteins"]:
                count_swissprot = payload["proteins"]["reviewed"]
            # count_trembl = payload["proteins"]["unreviewed"]

        if count_swissprot > 10:
            return True
        else:
            # print(signature, count_trembl)
            return False

    def process_sign(self, signature):
        url = f"https://www.ebi.ac.uk:443/interpro/api/protein/unreviewed/entry/{self.database}/{signature}/?page_size=200"
        list_pmid_acc = self.search_trembl_pmid(url)

        if len(list_pmid_acc) != 0:
            text_complete = f"{signature}, "
            for pmid, acc_list in list_pmid_acc.items():
                text_complete += f"{pmid}: "
                text = "; ".join(acc_list)
                text_complete += f"{text} | "
            text_complete = text_complete.strip(" | ")
            return text_complete
        return

    def search_trembl_pmid(self, BASE_URL):
        # disable SSL verification to avoid config issues
        context = ssl._create_unverified_context()

        next = BASE_URL
        attempts = 0

        list_pmid_acc = dict()

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
                        print("LAST URL: " + next)
                        print(e)
                        next = ""

            for i, item in enumerate(payload["results"]):
                # get UniProt accession
                accession = item["metadata"]["accession"]
                # search for list of PMIDs
                list_pmid = self.search_pmid(accession)

                if len(list_pmid) != 0:
                    for pmid in list_pmid:
                        try:
                            list_pmid_acc[pmid].append(accession)
                        except KeyError:
                            list_pmid_acc[pmid] = [accession]

            # Don't overload the server, give it time before asking for more
            if next:
                sleep(1)

        return list_pmid_acc

    def search_pmid(self, accession):
        # disable SSL verification to avoid config issues
        context = ssl._create_unverified_context()
        next = f"https://www.ebi.ac.uk/proteins/api/proteins/{accession}"
        attempts = 0

        while next:
            try:
                sleep(5)
                req = request.Request(next, headers={"Accept": "application/json"})
                # res = request.urlopen(req, context=context)
                with request.urlopen(req) as res:
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
                        print("LAST URL: " + next)
                        print(e)
                        next = ""
            except http.client.RemoteDisconnected as e:
                # http.client.RemoteDisconnected
                if attempts < 3:
                    attempts += 1
                    sleep(61)
                    continue
                else:
                    print("LAST URL: " + next)
                    print(e)
                    next = ""

        list_pmids = set()
        pmid = ""
        if payload:
            for i, item in enumerate(payload["references"]):
                if "dbReferences" in item["citation"]:
                    title = item["citation"]["title"]
                    t = re.search(r"(gene|genome|Gene|Genome|Genomic|genomic|sequen|Sequen)", title)

                    for j, ref in enumerate(item["citation"]["dbReferences"]):
                        if ref["type"] == "PubMed":
                            pmid = ref["id"]
                            if not t and pmid not in self.boring_pmids:
                                list_pmids.add(pmid)

        return list_pmids


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")

    args = parser.parse_args()
    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}'': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    database = config["files"]["member_db"]
    inputf = config["files"]["inputfile"]
    outputf = config["files"]["outputfile"]
    boringpmidf = config["files"]["boringpmidfile"]

    # init
    process = memberdb_pmid(database, boringpmidf)

    if not os.path.isfile(inputf):
        parser.error(f"Error file not found '{inputf}'")

    # get list of panther accession from file
    print(f"Searching data for {database} signatures")

    with open(inputf, "r") as f:
        count = 0
        for line in f:
            signature = line.strip("\n")

            if process.has_swissprot(signature):
                pass
            else:
                count += 1
                process.sign_in.append(signature)
            if count == 100:
                break

    print(len(process.sign_in))
    # results = []
    # with Pool(10) as p:
    #     results = p.map(process.process_sign, process.sign_in)

    # print("Writing results in file")
    # with open(outputf, "w") as outf:
    #     for item in results:
    #         # print(item)
    #         if item != None:
    #             outf.write(f"{item}\n")

    with open(outputf, "a") as outf:
        for signature in process.sign_in:
            print(f"Processing {signature}")
            results = process.process_sign(signature)
            if results:
                outf.write(f"{results}\n")


# get list of panther signatures unintegrated without comments:
# select m.method_ac
# from interpro.method m
# left join interpro.entry2method e2m on m.method_ac=e2m.method_ac
# left join interpro.method_comment c on c.method_ac=m.method_ac
# where m.method_ac like 'PTHR%' and m.method_ac not like 'PTHR%:%' and c.status is null and e2m.method_ac is null
# ;
