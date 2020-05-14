#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script is a utility class for scripts on proteomes of key species (InterPro BBR grant)

"""

import cx_Oracle
import os
import sys
import traceback


class proteome:
    def __init__(self):
        self.connection = None
        self.cursor = None
        self.tax_id = None

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
        protein_list = [row[0] for row in self.cursor]
        # print(f"Found {len(protein_list)} accessions")
        return protein_list

    def close_connection(self):
         """
        Close database connection
        """
        
        self.cursor.close()
        self.connection.close()
