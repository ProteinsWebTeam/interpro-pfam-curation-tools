from urllib import request
from urllib.error import HTTPError
from time import sleep
import os, sys, re, json, ssl


def is_fragment(protein):
    # disable SSL verification to avoid config issues
    context = ssl._create_unverified_context()

    next = f"https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/{protein}"
    payload = dict()

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
                    raise e

        if "metadata" in payload:
            fragment = payload["metadata"]["is_fragment"]
            if fragment:
                return True
            else:
                return False
    print("failed")
    return


if __name__ == "__main__":

    input_file = "sorted_output_matches_1.csv"
    output_file = "sorted_output_matches_1_no_fragment.csv"
    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        for line in fin.readlines():
            line = line.strip()
            protein = line.split(",")[1]

            if not is_fragment(protein):
                fout.write(f"{line}\n")
