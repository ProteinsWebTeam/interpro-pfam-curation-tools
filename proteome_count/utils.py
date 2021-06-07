#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script is a utility class for scripts on proteomes of key species (InterPro BBR grant)

"""

import cx_Oracle
import os
import sys
import traceback
import subprocess
from urllib import request
import re


class proteome:
    def __init__(self):
        self.connection = None
        self.cursor = None
        self.tax_id = None
        self.proteome = None
        self.dir = None
        self.old_proteinlist = dict()
        self.new_proteinlist = dict()

    def getConnection(self, user, password, schema):
        """
		Set database connection
		"""

        connectString = "".join([user, "/", password, "@", schema])
        try:
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

    def chunks(self, l, n):
        """Yield chunks of size n from iterable.

		Args:
			l (Iterable): The iterable to chunk
			n (int): Maximum number of items in chunk

		Yields:
			Iterable: Tuples of length n or less (final bucket will have only the remaining items)

		"""
        for i in range(0, len(l), n):
            # Create an index range for l of n items:
            yield l[i : i + n]

    def search_taxid(self, organism):
        """
		Search taxid for specified organism

		Args:
			configfile: Configuration file containing database connection credentials
			organism: Organism to search the taxid for

		Yields:
			taxid: taxonomic identifier for the organism of interest

		"""

        request = "Select tax_id from INTERPRO.ETAXI where scientific_name=:1"
        self.cursor.execute(request, (organism,))
        try:
            self.tax_id = self.cursor.fetchone()[0]
            print(f"Found taxid {self.tax_id} for {organism}")
        except TypeError:
            print(f"Error: Taxid not found for {organism}")
            sys.exit(1)

    def get_proteins(self):
        """
		Search protein accessions for the taxid

		Args: None

		Yields: 
			protein_list: list of proteins

		"""

        # print(f"Searching protein accessions for taxid {self.tax_id}")

        request = "SELECT UNIQUE P.PROTEIN_AC \
			FROM INTERPRO.PROTEIN P \
			JOIN INTERPRO.ETAXI ET ON P.TAX_ID = ET.TAX_ID \
			WHERE P.TAX_ID=:1"

        self.cursor.execute(request, (self.tax_id,))
        protein_list = set(row[0] for row in self.cursor)
        # print(f"Found {len(protein_list)} accessions")
        return protein_list

    def get_proteins_from_uniprot(self):
        """
		Search protein accessions for the reference proteome

		Args: None

		Yields: 
			protein_list: list of proteins

		"""
        uniprot_file = os.path.join(self.dir, f"uniprot-proteome_{self.proteome}.list")
        if not os.path.isfile(uniprot_file):
            print(f"Downloading protein file for {self.proteome} from UniProt)")
            proteome_file = (
                f"https://www.uniprot.org/uniprot/?query=proteome:{self.proteome}&format=list"
            )
            wget = subprocess.Popen(["wget", "-O", uniprot_file, proteome_file])  # ~5min
            wget.communicate()

        print(f"Loading list of proteins for proteome {self.proteome} into memory")
        with open(uniprot_file, "r") as f:
            proteinlist = set()
            for line in f:
                line = line.strip("\n")
                proteinlist.add(line)

        return proteinlist

    def get_protein2ipr_file(self, version):
        """
		Search protein accessions in InterPro entries

		Args: version: InterPro release version

		Yields: 
			protein_list: list of proteins

		"""
        protein_file = os.path.join(self.dir, f"proteinlist_{version}")
        compressed_file = os.path.join(self.dir, f"{version}_protein2ipr.dat.gz")
        if not os.path.isfile(protein_file) or os.path.getsize(protein_file) == 0:
            print(f"Downloading protein file for InterPro {version} from ftp")
            if not os.path.isfile(compressed_file) or os.path.getsize(compressed_file) == 0:
                interpro_file = (
                    f"ftp://ftp.ebi.ac.uk/pub/databases/interpro/{version}/protein2ipr.dat.gz"
                )
                wget = subprocess.Popen(["wget", "-O", compressed_file, interpro_file])  # ~5min
                wget.communicate()

            cmd = f"zcat {compressed_file} | cut -f1 | sort | uniq > {protein_file}"
            zcat = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE)  # ~ 1hour
            zcat.communicate()

        print(f"Loading list of integrated proteins for InterPro {version} into memory")
        with open(protein_file, "r") as f:
            proteinlist = set()
            for line in f:
                line = line.strip("\n")
                proteinlist.add(line)

        return proteinlist

    def close_connection(self):
        """
		Close database connection
		"""

        self.cursor.close()
        self.connection.close()
