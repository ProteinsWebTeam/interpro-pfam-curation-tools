#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script searches marine species taxid in InterPro

@arguments [CONFIG_FILE]: file containing path to directories to process

"""

import argparse
import os
from configparser import ConfigParser
import sys
import requests
from time import sleep
from multiprocessing import Pool


def get_sequence(uniprot):
    BASE_URL = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot}"

    r = requests.get(BASE_URL, headers={"Accept": "text/x-fasta"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    return r.text


def search_sequences_from_file(filein, outputfile):
    if not os.path.isfile(outputfile) or os.path.getsize(outputfile) == 0:
        with open(filein, "r") as fin, open(outputfile, "w") as fout:
            for line in fin.readlines():
                uniprot = line.strip()
                fout.write(get_sequence(uniprot))
    return filein


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)
    inputdir = config["files"]["output_nomatches"]
    outputdir = config["files"]["output_fasta"]
    os.makedirs(outputdir, exist_ok=True)

    filesin = os.listdir(inputdir)

    list_files = []
    for f in filesin:
        fpath = os.path.join(inputdir, f)
        fout = f[18:-4]
        outputf = os.path.join(outputdir, f"{fout}.fasta")
        finout = (fpath, outputf)
        list_files.append((fpath, outputf))

    print("Searching sequences")
    with Pool(10) as p:
        results = p.starmap(search_sequences_from_file, list_files)
    print("Search completed successfully")

