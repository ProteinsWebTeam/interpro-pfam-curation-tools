#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script gets the unintegrated Uniprot coverage for a given InterPro member database

@arguments [-u USER]: database user (if -r taxid)
           [-p PASSWORD]: database password for the user (if -r taxid)
           [-s SCHEMA]: database schema to use (if -r taxid)
           [-d DATABASE]: Member database to search coverage on

@usage Start interactive shell with 20K memory: bsub -M 20000 -R"rusage[mem=20000]" -n 5 -Is $SHELL
       source oracle 11.2: source ~oracle/ora112setup.sh
       run script: python search_unintegrated_cov.py -u interpro -p xxx -s IPPRO -d MEMBERDB

"""

import argparse
import traceback
import cx_Oracle
import sys
import os


class coverage:
    def __init__(self, user, password, schema, memberdb):
        self.getConnection(user, password, schema)
        self.memberdb = memberdb
        self.result_file = f"result_file_{self.memberdb}.csv"
        self.list_abrev = {
            "PANTHER": "PTHR",
            "CDD": "cd",
            "GENE3D": "G3DSA",
            "PFAM": "pf",
            "SFLD": "SFLD",
            "PROSITE profile": "PS",
            "PROSITE patterns": "PS",
            "SUPERFAMILY": "SF",
            "TIGRFAMs": "TIGR",
            "HAMAP": "MF",
            "PIRSF": "PIRSF",
            "PRINTS": "PR",
            "SMART": "SM",
        }

        if os.path.isfile(self.result_file):
            os.remove(self.result_file)
        with open(self.result_file, "w") as f:
            f.write("signature, unintegrated_proteins, total_proteins, total_swissprot\n")

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

    def close_connection(self):
        """
		Close database connection
		"""

        self.cursor.close()
        self.connection.close()

    def get_unintegrated_signatures(self):
        temp_text = ""
        if self.memberdb == "PANTHER":
            temp_text = "AND M.METHOD_AC NOT LIKE 'PTHR%:SF%' "
        request = f"""SELECT DISTINCT M.METHOD_AC, M.PROTEIN_AC, P.DBCODE
                    FROM INTERPRO.MATCH M
                    LEFT OUTER JOIN INTERPRO.ENTRY2METHOD E2M ON M.METHOD_AC = E2M.METHOD_AC
                    LEFT OUTER JOIN METHOD_COMMENT EC ON EC.METHOD_AC=M.METHOD_AC
                    JOIN INTERPRO.PROTEIN P ON P.PROTEIN_AC=M.PROTEIN_AC
                    WHERE (M.METHOD_AC LIKE '{self.list_abrev[self.memberdb]}%' {temp_text}AND E2M.METHOD_AC IS NULL) AND EC.METHOD_AC IS NULL"""

        self.cursor.execute(request)
        self.signatures_prot_list = dict()
        self.signatures_swiss = dict()

        for row in self.cursor:
            try:
                self.signatures_prot_list[row[0]].append(row[1])
            except:
                self.signatures_prot_list[row[0]] = [row[1]]

            if row[2] == "S":
                try:
                    self.signatures_swiss[row[0]] += 1
                except:
                    self.signatures_swiss[row[0]] = 1

        # print(len(self.signatures_prot_list))

    def get_unintegrated_proteins(self):
        # This request get the counts matching unintegrated signatures per protein
        request = """SELECT SUB1.PROTEIN_AC
        FROM (
        SELECT PROTEIN_AC, COUNT(DISTINCT METHOD_AC) AS ALL_COUNT
        FROM INTERPRO.MATCH
        GROUP BY PROTEIN_AC) SUB1
        JOIN (
        SELECT M.PROTEIN_AC, COUNT(DISTINCT M.METHOD_AC) AS UN_COUNT
        FROM INTERPRO.MATCH M
        LEFT JOIN INTERPRO.ENTRY2METHOD EM ON M.METHOD_AC = EM.METHOD_AC
        WHERE EM.ENTRY_AC IS NULL
        GROUP BY M.PROTEIN_AC) SUB2 ON SUB2.PROTEIN_AC=SUB1.PROTEIN_AC
        WHERE SUB1.ALL_COUNT-SUB2.UN_COUNT=0
        """
        self.cursor.execute(request)
        self.un_uniprot_list = {row[0]: 0 for row in self.cursor}

    def get_unintegrated_proteins_per_signatures(self):
        self.signatures_un_prot = dict()

        for signature, protein_list in self.signatures_prot_list.items():
            count_un = 0
            for protein in protein_list:
                if protein in self.un_uniprot_list:
                    count_un += 1
            if count_un > 0:
                swiss_count = (
                    self.signatures_swiss[signature] if signature in self.signatures_swiss else 0
                )
                self.signatures_un_prot[signature] = {
                    "un_prot": count_un,
                    "prot": len(protein_list),
                    "swiss": swiss_count,
                }
                self.print_to_file(signature, count_un, len(protein_list), swiss_count)

        # print(self.signatures_un_prot)

    def print_to_file(self, signature, count_un, totelcount, swisscount):
        with open(self.result_file, "a") as f:
            f.write(f"{signature}, {count_un}, {totelcount}, {swisscount}\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--user", help="username for database connection", required=True)
    parser.add_argument("-p", "--password", help="password for database connection", required=True)
    parser.add_argument("-s", "--schema", help="database schema to connect to", required=True)
    parser.add_argument(
        "-d",
        "--database",
        help="Member database to search coverage on",
        choices=[
            "PANTHER",
            "CDD",
            "GENE3D",
            "PFAM",
            "SFLD",
            "PROSITE profile",
            "PROSITE patterns",
            "SUPERFAMILY",
            "TIGRFAMs",
            "HAMAP",
            "PIRSF",
            "PRINTS",
            "SMART",
        ],
        required=True,
    )

    args = parser.parse_args()

    cov = coverage(args.user, args.password, args.schema, args.database)
    print(f"Searching data for {cov.memberdb}")
    cov.get_unintegrated_signatures()
    print("Searching unintegrated proteins")
    cov.get_unintegrated_proteins()
    print("Searching signatures with unintegrated proteins")
    cov.get_unintegrated_proteins_per_signatures()
