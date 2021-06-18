import sys, os, json, ssl, re
from urllib import request
from urllib.error import HTTPError
from time import sleep
import argparse
from configparser import ConfigParser
from get_memberdb_pmid import memberdb_pmid


class signature_pmid(memberdb_pmid):
    def __init__(self, member_db, boringfile, signature):
        super().__init__(member_db, boringfile)
        self.signature = signature


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="CONFIG_FILE", help="configuration file")
    parser.add_argument("signature", help="signature to process")
    args = parser.parse_args()
    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}'': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    database = config["files"]["member_db"]
    boringpmidf = config["files"]["boringpmidfile"]

    # init
    process = signature_pmid(database, boringpmidf, args.signature)
    if process.has_swissprot(process.signature):
        print(f"Signature '{process.signature}' also matching SwissProt")
    else:
        text = process.process_sign(process.signature)
        if text:
            print(f"PMID found for {process.signature}: {text}")
        else:
            print(f"No PMID found for {process.signature}")
