#!/usr/bin/env python3

from deduf_utils import pfam_duf
import sys, errno, re, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep
import os


class pfam_swiss(pfam_duf):
    def __init__(self):
        super().__init__()
        self.list_duf_unintegrated = dict()

    def get_pfam_duf_list_unintegrated(self):
        request = "select m.method_ac, m.name \
                    from interpro.method m \
                    left join interpro.entry2method e2m on e2m.method_ac=m.method_ac \
                    where m.method_ac like 'PF%' and m.name like 'DUF%' and e2m.entry_ac is null"

        self.cursor.execute(request)
        self.list_duf_unintegrated = {
            row[0]: {"dufid": row[1], "ipr": None, "name": None} for row in self.cursor
        }

    def get_list_swissprot_names(self):
        for pfamid in self.list_duf:
            self.get_list_swissprot_names_pfam(pfamid)

    def get_list_swissprot_names_pfam(self, pfamid):

        BASE_URL = f"https://www.ebi.ac.uk:443/interpro/api/protein/reviewed/entry/pfam/{pfamid}/?page_size=200"
        swissprot_names = self.output_list(BASE_URL)
        self.list_duf[pfamid]["swissprot"] = swissprot_names
        self.list_duf[pfamid]["swissprot_count"] = len(swissprot_names)

    def save_swissprot_in_file(self, outputdir):
        outputfile = os.path.join(outputdir, "duf_swiss_names.csv")
        with open(outputfile, "w") as outf:
            outf.write(
                "Status,Pfam identifier\tPfam short name\tInterPro identifier\tAVG number of domains per InterPro\tInterPro entry name\tSwissprot count\tList swissprot names\n"
            )

            for pfamid, content in self.list_duf.items():
                # keep DUF with SwissProt matches only
                if content["swissprot_count"] != 0:
                    # keep DUF with SwissProt matches not Uncharacterized only
                    swissprotnames = ""
                    for name in content["swissprot"]:
                        if not re.search("[Uu]ncharacterized", name) and not re.search("UPF", name):
                            swissprotnames += f"{name} | "
                    if swissprotnames != "":
                        ipr = content["ipr"]
                        ipr_count = self.count_interpro_dom[ipr]
                        outf.write(
                            f"\t{pfamid}\t{content['dufid']}\t{ipr}\t{ipr_count}\t{content['name']}\t{content['swissprot_count']}\t{swissprotnames}\n"
                        )

    def output_list(self, BASE_URL):
        # disable SSL verification to avoid config issues
        context = ssl._create_unverified_context()

        next = BASE_URL
        last_page = False

        attempts = 0
        list_names = set()

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

            for i, item in enumerate(payload["results"]):

                # add to dictionary {accession:name}
                # list_names[item["metadata"]["accession"]] = item["metadata"]["name"]

                # add name to set of swissprot names
                list_names.add(item["metadata"]["name"])

            # Don't overload the server, give it time before asking for more
            if next:
                sleep(1)

        return list_names
