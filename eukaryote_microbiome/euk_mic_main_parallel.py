#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script searches marine species taxid in InterPro and produces fasta files when no matches are found

@arguments CONFIG_FILE: file containing the database credentials and output directories locations (see config.ini)

"""

import argparse
import os
from configparser import ConfigParser
import cx_Oracle
import traceback

import sys, re, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep
from multiprocessing import Pool
from get_sequences import get_sequence


class eukaryote_microbiome:
    def __init__(self, output_matches, outpout_no_matches, output_fasta, processed):
        self.signature_matches = dict()
        self.no_matches = set()
        self.processed_taxid = set()
        self.processed_file = processed
        self.outputm = output_matches
        self.outputnm = outpout_no_matches
        self.outputf = output_fasta
        self.withcomments = set()

    def getConnection(self, user, password, schema):
        """
		Set database connection
		"""

        connectString = "".join([user, "/", password, "@", schema])
        try:
            # subprocess.run(["source", "~oracle/ora112setup.sh"])
            self.connection = cx_Oracle.connect(connectString)
            self.cursor = self.connection.cursor()
        except:
            stackTrace = traceback.format_exc()
            print(stackTrace)
            if "invalid username" in stackTrace:
                print("Could not connect to {0} as user {1}".format(schema, user))
                print(
                    "NB if your oracle username contains the '$' character either escape it or surround it with quotes"
                )
                print('eg "ops$craigm" or ops\$craigm')
                print("Otherwise the shell will remove the '$' and all subsequent characters!")
                sys.exit(1)

    def close_connection(self):
        """
		Close database connection
		"""

        self.cursor.close()
        self.connection.close()

    def search_proteins(self, taxid):
        # disable SSL verification to avoid config issues
        context = ssl._create_unverified_context()

        next = f"https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/taxonomy/uniprot/{taxid}?page_size=200"
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
                        print(f"Couldn't find proteins for taxid {taxid}")
                        next = ""

            if "results" in payload:
                if payload["previous"] == None:
                    print(f"{payload['count']} proteins found for {taxid}")

                for i, item in enumerate(payload["results"]):

                    protein = item["metadata"]["accession"]
                    self.search_interpro_mapping(protein)

            # Don't overload the server, give it time before asking for more
            if next:
                sleep(1)

    def search_interpro_mapping(self, protein):
        if not self.is_fragment(protein):
            # disable SSL verification to avoid config issues
            context = ssl._create_unverified_context()

            next = f"https://www.ebi.ac.uk:443/interpro/api/entry/all/protein/UniProt/{protein}"
            payload = dict()
            interpro_acc = set()
            memberdb_acc = set()

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
                            raise e

                if "results" in payload:
                    for i, item in enumerate(payload["results"]):
                        acc = item["metadata"]["accession"]
                        # print(protein, acc)
                        if re.search("IPR", acc):
                            interpro_acc.add(acc)
                        else:
                            # do not take into account PANTHER subfamilies
                            if not re.search(r":SF", acc):
                                memberdb_acc.add(acc)

                # Don't overload the server, give it time before asking for more
                if next:
                    sleep(1)

            # protein not integrated in InterPro entry
            if len(interpro_acc) == 0:
                # match to signature found
                if len(memberdb_acc) != 0:
                    # print(f"unintegrated {protein} found in signatures {memberdb_acc}")
                    self.signature_matches[protein] = []
                    for signature in memberdb_acc:
                        if signature not in self.withcomments and not self.has_comment(signature):
                            count_uniprot, count_swissprot = self.search_memberdb_protein_count(
                                signature
                            )
                            if count_uniprot != 0:
                                result = f"{signature} ({count_uniprot} / {count_swissprot})"
                                self.signature_matches[protein].append(result)
                    if len(self.signature_matches[protein]) == 0:
                        del self.signature_matches[protein]

                # no signature matches
                else:
                    # print(f"no matches found for {protein}")
                    self.no_matches.add(protein)

    def is_fragment(self, protein):
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
        return

    def has_comment(self, signature):
        request = f"select mc.value\
                    from interpro.method m \
                    left join interpro.method_comment mc on mc.method_ac=m.method_ac \
                    where m.method_ac='{signature}' and mc.status='Y'"

        self.cursor.execute(request)
        comments = []
        for row in self.cursor:
            comments.append(row[0])
        if len(comments) > 0:
            self.withcomments.add(signature)
            return True

        return False

    def search_memberdb_protein_count(self, signature):
        # disable SSL verification to avoid config issues
        context = ssl._create_unverified_context()

        count_uniprot = 0
        count_swissprot = 0

        memberdb = ""
        if re.search(r"PTHR", signature):
            memberdb = "panther"
        elif re.search(r"G3DSA", signature):
            memberdb = "cathgene3d"
        elif re.search(r"SSF", signature):
            memberdb = "ssf"
        elif re.search(r"PF", signature):
            memberdb = "pfam"
        elif re.search(r"cd", signature):
            memberdb = "cdd"
        elif re.search(r"SM", signature):
            memberdb = "smart"
        elif re.search(r"SFLD", signature):
            memberdb = "sfld"
        elif re.search(r"TIGR", signature):
            memberdb = "tigrfams"
        elif re.search(r"MF", signature):
            memberdb = "hamap"
        elif re.search(r"PR", signature):
            memberdb = "prints"
        elif re.search(r"PS", signature):
            memberdb = "profile"
        elif re.search(r"PIRSF", signature):
            memberdb = "pirsf"

        next = f"https://www.ebi.ac.uk:443/interpro/api/protein/entry/{memberdb}/{signature}"
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
                    if memberdb == "profile":
                        memberdb = "prosite"
                        next = f"https://www.ebi.ac.uk:443/interpro/api/protein/entry/{memberdb}/{signature}"
                        continue
                    else:
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

            if "proteins" in payload:
                for key, value in payload["proteins"].items():
                    if key == "uniprot":
                        count_uniprot = value
                    elif key == "reviewed":
                        count_swissprot = value

            # Don't overload the server, give it time before asking for more
            if next:
                sleep(1)

        return count_uniprot, count_swissprot

    def write_in_file(self, taxid):
        print(f"Writting results in file for {taxid}")
        if len(self.signature_matches) > 0:
            output_matches = os.path.join(self.outputm, f"output_matches_{taxid}.csv")

            # columns in file: taxid, uniprot, memberdb (protein count / swissprot count)
            with open(output_matches, "w") as f:
                for key, members in self.signature_matches.items():
                    memberdb = " | ".join(members)
                    f.write(f"{taxid},{key},{memberdb}\n")

        if len(self.no_matches) > 0:
            output_no_matches = os.path.join(self.outputnm, f"output_no_matches_{taxid}.csv")
            output_fasta = os.path.join(self.outputf, f"{taxid}.fasta")

            with open(output_no_matches, "w") as f, open(output_fasta, "w") as fasta:
                for protein in self.no_matches:
                    f.write(f"{protein}\n")
                    seq = get_sequence(protein)
                    fasta.write(seq)

        else:
            print(f"No result found for {taxid}")

        with open(self.processed_file, "a") as f:
            f.write(f"{taxid}\n")

        self.signature_matches = dict()
        self.no_matches = set()


def process_taxid(
    outputmatches, outnomatches, outputfasta, processed_file, user, password, schema, taxid
):
    euk = eukaryote_microbiome(outputmatches, outnomatches, outputfasta, processed_file)

    if os.path.isfile(processed_file):
        with open(processed_file, "r") as f:
            for line in f.readlines():
                taxidtmp = line.strip()
                euk.processed_taxid.add(taxidtmp)

    if taxid not in euk.processed_taxid:
        print(f"Processing {taxid}")
        sys.stdout.flush()
        euk.getConnection(user, password, schema)
        # search proteins matches and find out if proteins are integrated in an InterPro entry or not
        euk.search_proteins(taxid)
        euk.write_in_file(taxid)
        euk.close_connection()
        return taxid

    return


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    processed_file = config["files"]["outputfile_processed"]
    user = config["database"]["user"]
    password = config["database"]["password"]
    schema = config["database"]["schema"]
    taxid_file = config["files"]["inputfile"]
    output_matches = config["files"]["output_matches"]
    output_no_matches = config["files"]["output_nomatches"]
    output_fasta = config["files"]["output_fasta"]

    if not os.path.isfile(taxid_file):
        print(f"Cannot open '{taxid_file}': " f"no such file or directory")
        sys.exit(1)

    os.makedirs(output_matches, exist_ok=True)
    os.makedirs(output_no_matches, exist_ok=True)
    os.makedirs(output_fasta, exist_ok=True)

    print("Searching proteins")
    listtaxid = []
    with open(taxid_file, "r") as f:
        for line in f.readlines():
            taxid = line.strip()
            tuple = (
                output_matches,
                output_no_matches,
                output_fasta,
                processed_file,
                user,
                password,
                schema,
                taxid,
            )
            listtaxid.append(tuple)

    # print(listtaxid)
    with Pool(10) as p:
        results = p.starmap(process_taxid, listtaxid)

    # need to concat output_matches files
    # find ./output/matches -name 'output_matches*' -exec cat {} + > output_matches_1.csv
    # remove lines with duplicate signatures matches
    # sort -u -k3 output_matches_1.csv > sorted_output_matches_1.csv
