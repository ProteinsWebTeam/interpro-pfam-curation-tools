from deduf_utils import pfam_duf
import sys, errno, re, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep


class pfam_litterature(pfam_duf):
    def __init__(self):
        super().__init__()

    def get_list_articles_pfamid(self, pfamid):
        context = ssl._create_unverified_context()

        next = f"https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByEntity?entity={pfamid}&filter=1&format=JSON&pageSize=4"

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
                if len(payload["articles"]) > 0:
                    if payload["nextCursorMark"] != -1:
                        cursor = payload["nextCursorMark"]
                        next = f"https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByEntity?entity={pfamid}&filter=1&format=JSON&pageSize=4&cursorMark={cursor}"
                    else:
                        next = ""
                else:
                    next = ""
                attempts = 0
                if not next:
                    last_page = True
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

            if len(payload["articles"]) != 0:
                for i, item in enumerate(payload["articles"]):
                    # print(item)
                    pmid = item["extId"]
                    # annotation = []
                    # for j, ann in enumerate(item["annotations"]):
                    #     prefix = ""
                    #     postfix = ""
                    #     try:
                    #         prefix = ann["prefix"]
                    #         postfix = ann["postfix"]
                    #         annotation.append(f"{prefix} {pfamid} {postfix}")
                    #     except KeyError as e:
                    #         annotation.append(
                    #             f"{pfamid} Term can be listed but not highlighted within the text"
                    #         )

                    # list_pmid
                    list_pmid.add(pmid)

                    # list_pmid.add(f"https://europepmc.org/article/MED/{pmid}")

            # Don't overload the server, give it time before asking for more
            if next:
                sleep(1)

        self.list_duf[pfamid]["articles_pfam"] = list_pmid

    def get_list_articles_duf(self, pfamid):

        dufid = self.list_duf[pfamid]["dufid"]

        # disable SSL verification to avoid config issues
        context = ssl._create_unverified_context()

        next = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={dufid}&cursorMark=*&pageSize=50&format=json"
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
                if payload["nextCursorMark"] != payload["request"]["cursorMark"]:
                    cursor = payload["nextCursorMark"]
                    next = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={dufid}&cursorMark={cursor}&pageSize=50&format=json"
                else:
                    next = ""
                attempts = 0
                if not next:
                    last_page = True
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

            if len(payload["resultList"]["result"]) != 0:
                for i, item in enumerate(payload["resultList"]["result"]):

                    # add name to set of swissprot names
                    if item["source"] == "MED":
                        # print(item)
                        pmid = item["pmid"]
                        list_pmid.add(pmid)
                        # list_pmid.add(f"https://europepmc.org/article/MED/{pmid}")

            # Don't overload the server, give it time before asking for more
            if next:
                sleep(1)
        self.list_duf[pfamid]["articles_duf"] = list_pmid

    def save_result_in_file(self, outputfile):
        with open(outputfile, "w") as outf:
            outf.write(
                "Pfam identifier,Pfam short name,InterPro identifier,count pmid,EuropePMC link\n"
            )
            for pfamid, item in self.list_duf.items():
                dufid = item["dufid"]
                ipr = item["ipr"]
                set_pfam_id = set()

                list_pfam = item["articles_pfam"]
                list_duf = item["articles_duf"]

                list_pmdi = set(list_pfam).union(set(list_duf))
                nb_pmid = len(list_pmdi)
                if nb_pmid > 0:
                    print(pfamid, nb_pmid)

                    link = f"https://europepmc.org/search?query={dufid}%20or%20{pfamid}"
                    outf.write(f"{pfamid},{dufid},{ipr},{nb_pmid},{link}\n")
                # for pmid, matches in item["articles_pfam"].items():
                #     set_pfam_id.add(pmid)
                #     outf.write(f"{pfamid},{dufid},{ipr},{pmid},{matches}\n")

                # missing_pmids = item["articles_duf"] - set_pfam_id
                # for pm in missing_pmids:
                #     outf.write(f"{pfamid},{dufid},{ipr},{pm}\n")

