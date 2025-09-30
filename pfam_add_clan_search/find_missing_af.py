"""
This script performs the following steps:

1.  Extract Pfam representatives from the output directory of chopped alphafold structures.
2.  Extract Pfam accessions from the Pfam SQL database.
3.  Compare the above lists and generate a list of Pfam accessions without a shopped AlphaFold structure.

Usage: python3 find_multi_af.py CONFIG_FILE
"""

from find_similar_struct import extract_pfams, select_representative_protein
from utils import sql_connection
from configparser import ConfigParser
import argparse
import os

def get_processed_pfams(file):
    pfams = set()
    with open(file, 'r') as f:
        for line in f.readlines():
            pfams.add(line.strip())
    return pfams

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="CONFIG_FILE", help="configuration file")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    #mysql connection info
    host = config["db"]["host"]
    user = config["db"]["user"]
    password = config["db"]["password"]
    database = config["db"]["database"]
    port = config["db"]["port"]
    dbinfo = [host, user, password, database, port]

    pfams = extract_pfams(dbinfo)
    processed_pfams = get_processed_pfams(config["files"]["processed"])

    all_pfams = set()    
    for pfam in pfams:
        pfam_acc = pfam['pfamA_acc']
        all_pfams.add(pfam_acc)

    unprocessed_pfams = all_pfams - processed_pfams
    print(len(unprocessed_pfams))

    with open(config["files"]["to_process"] , 'w') as f:
        for pfam in unprocessed_pfams:
            for pfam_acc in pfams:
                if pfam == pfam_acc['pfamA_acc']:
                    f.write(f"{pfam}\t{pfam_acc['type']}\t{pfam_acc['clan_acc']}\n")

if __name__ == "__main__":
    main()
