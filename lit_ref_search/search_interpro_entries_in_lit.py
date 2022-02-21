import sys, errno, re, json, ssl
import argparse
import os
from configparser import ConfigParser
from urllib import request
from urllib.error import HTTPError
from time import sleep


def get_list_articles_entry_id(entry_id):
    context = ssl._create_unverified_context()

    next = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={entry_id}&resultType=lite&cursorMark=*&pageSize=100&format=json&fromSearchPost=false"

    last_page = False

    attempts = 0
    list_pmid = set()

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

        for i, item in enumerate(payload["resultList"]["result"]):
            pmid = item["pmid"]
            link = f"https://europepmc.org/article/MED/{pmid}"

            list_pmid.add(link)

        # Don't overload the server, give it time before asking for more
        if next:
            sleep(1)

        return " | ".join(list_pmid)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "entries_file",
        metavar="FILE",
        help="file containing list of interpro entries to search literature references for",
    )
    parser.add_argument("output_file", metavar="FILE", help="outputfile path")

    args = parser.parse_args()

    if os.path.isfile(args.entries_file):
        with open(args.entries_file, "r") as fin, open(args.output_file, "w") as fout:
            for line in fin.readlines():
                entry = line.strip()
                pmid = get_list_articles_entry_id(entry)
                fout.write(f"{entry},{pmid}\n")
