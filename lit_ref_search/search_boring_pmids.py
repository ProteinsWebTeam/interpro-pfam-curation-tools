# A script to identify how many times each paper in UNiprot is
# referred to This identifies papers which are reffered by a lot of
# entries and are probably not much use for annotation.

import argparse
import os
import sys
import gzip
import re

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("rel_version", help="configuration file")
    args = parser.parse_args()

    # Open Uniprot files
    uniprot_file = f"/homes/agb/Work/UniProt_analysis/UniProt_release_dir/uniprot{args.rel_version}/uniprot_trembl.dat.gz"

    if not os.path.isfile(uniprot_file):
        print(f"Error file not found {uniprot_file}")
        sys.exit(1)

    pmid_dict = dict()
    boring_pmid = []

    print("Searching PMIDs")

    with gzip.open(uniprot_file, "rb") as f:
        for line in f:
            line = line.decode("utf-8").strip("\n")
            m = re.search(r"PubMed=(\d+)", line)
            if m:
                try:
                    pmid_dict[m.group(1)] += 1
                except KeyError:
                    pmid_dict[m.group(1)] = 1

    print("Looking for boring PMIDs")
    with open("boring_pmids_2", "w") as f:
        for pmid, count in pmid_dict.items():
            if count > 9:
                boring_pmid.append(pmid)
                f.write(f"{pmid}\t{count}\n")

