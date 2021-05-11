from search_swissprot import pfam_swiss
import json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep


class pfam_go(pfam_swiss):
    def __init__(self):
        super().__init__()

    def get_list_duf_from_file(self, duf_file):
        with open(duf_file, "r") as f:
            for line in f:
                # print(line)
                line = line.split("\t")
                # print(line)
                pfamid = line[0]
                dufid = line[1]
                nb_dom = line[3]
                if nb_dom == "1":
                    # print(pfamid, dufid, nb_dom)
                    self.list_duf[pfamid] = {"dufid": dufid}

    def get_swissprot_accessions(self):
        print("Searching SwissProt matches")

        for pfamid in self.list_duf:
            BASE_URL = f"https://www.ebi.ac.uk:443/interpro/api/protein/reviewed/entry/pfam/{pfamid}/?page_size=200"
            swissprot_names, swissprot_acc, count_swissprot = self.output_list(BASE_URL)
            self.list_duf[pfamid]["swissprot_acc"] = swissprot_acc

    def get_go_annotation(self):
        print("Searching GO and Keyword annotations")
        for pfamid in self.list_duf:
            list_go = set()
            list_keywords = set()

            for acc in self.list_duf[pfamid]["swissprot_acc"]:
                BASE_URL = f"https://www.ebi.ac.uk/proteins/api/proteins/{acc}"

                attempts = 0

                # disable SSL verification to avoid config issues
                context = ssl._create_unverified_context()

                try:
                    req = request.Request(BASE_URL, headers={"Accept": "application/json"})
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
                            raise e

                for i, item in enumerate(payload["dbReferences"]):
                    if item["type"] == "GO":
                        name = item["properties"]["term"]
                        category = name.split(":")[0]
                        # only save GO terms related to Biological process or Molecular function
                        if category == "P" or category == "F":
                            term = f"{item['id']}_{name}"
                            list_go.add(term)

                if "keywords" in payload:
                    for i, item in enumerate(payload["keywords"]):
                        keyword = item["value"]
                        if keyword != "Reference proteome":
                            list_keywords.add(keyword)

            self.list_duf[pfamid]["go_terms"] = list_go
            self.list_duf[pfamid]["keywords"] = list_keywords
            print(pfamid, list_go, list_keywords)

    def save_results_in_file(self, outputfile):
        print(f"Saving results in {outputfile}")

        with open(outputfile, "w") as outf:
            outf.write("Status\tPfam identifier\tPfam short name\tList GO terms\tList keywords\n")

            for pfamid, content in self.list_duf.items():
                list_go = " | ".join(content["go_terms"])
                list_keywords = " | ".join(content["keywords"])
                if list_go != "" or list_keywords != "":
                    outf.write(f"\t{pfamid}\t{content['dufid']}\t{list_go}\t{list_keywords}\n")
