#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script searches unintegrated methods for CATH database from an input file given by the CATH team 
        and print the results in an output file 

@arguments [-u USER]: database user
           [-p PASSWORD]: database password for the user
           [-s SCHEMA]: database schema to use
           [-f INPUTFILE]: cathb file containing new data
           [-o OUTPUTFILE]: output file to write info to
           [-c COMMENTFILE]: [optional] file containing previous comments on the methods

"""

import argparse
import cx_Oracle
import os
import sys
import re
import traceback


class new_cathb:
    def __init__(self, user, password, schema):
        self.getConnection(user, password, schema)

        self.getOldNames()
        self.getProteinCounts()
        self.comments = dict()
        self.unintegrated = dict()

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

    def getOldNames(self):
        """
        Search for previous names given by CATH in the database for each method
        """
        self.descriptions = dict()
        request = """SELECT m.METHOD_AC, m.DESCRIPTION
                    FROM INTERPRO.METHOD m 
                    WHERE m.DESCRIPTION is not null and m.METHOD_AC like 'G3DSA%'
                """

        self.cursor.execute(request)
        results = self.cursor.fetchall()

        for row in results:
            self.descriptions[row[0]] = row[1]

    def getProteinCounts(self):
        """
        Search for the number of proteins for each method
        """
        self.protein_counts = dict()
        request = """SELECT METHOD_AC, count(PROTEIN_AC)
                    FROM INTERPRO.MATCH
                    WHERE METHOD_AC like 'G3DSA%'
                    GROUP BY METHOD_AC
                """

        self.cursor.execute(request)
        results = self.cursor.fetchall()

        for row in results:
            self.protein_counts[row[0]] = row[1]

    def getComments(self, commentfile):
        """
        Search for the curation comments previously assigned to each method
        """
        pattern = r"(G3DSA:[\d\.]+),(.+)"
        try:
            with open(commentfile, "r") as f:
                for line in f:
                    line = line.strip("\n")
                    print(line)
                    m = re.search(pattern, line)
                    sfid = m.group(1)
                    if not sfid in self.comments:
                        self.comments[sfid] = [m.group(2)]
                    else:
                        self.comments[sfid].append(m.group(2))
        except FileNotFoundError as e:
            print(e)
            sys.exit()

    def getUnintegrated(self):
        """
        Search for integrated methods in unchecked InterPro entries with the curators comments
        """
        unintegrated = dict()
        request = """SELECT m.METHOD_AC, e2m.ENTRY_AC
                    FROM METHOD m
                    LEFT JOIN INTERPRO.ENTRY2METHOD e2m on m.METHOD_AC=e2m.METHOD_AC
                    LEFT JOIN INTERPRO.ENTRY e on e2m.ENTRY_AC=e.ENTRY_AC
                    WHERE m.METHOD_AC like 'G3DSA%' and e2m.METHOD_AC is null
            UNION
                    SELECT m.METHOD_AC, e2m.ENTRY_AC
                    FROM METHOD m
                    LEFT JOIN INTERPRO.ENTRY2METHOD e2m on m.METHOD_AC=e2m.METHOD_AC
                    LEFT JOIN INTERPRO.ENTRY e on e2m.ENTRY_AC=e.ENTRY_AC
                    WHERE m.METHOD_AC like 'G3DSA%' and e.CHECKED='N'
        """

        self.cursor.execute(request)
        results = self.cursor.fetchall()

        for row in results:
            sfid = row[0]
            ipr = row[1]

            if sfid in self.comments:
                comments = self.comments[sfid]
            else:
                comments = []
            self.unintegrated[sfid] = {"comment": comments, "entry": ipr}

    def getSfInfo(self, inputfile):
        """
        Get CATH new names for unintegrated methods
        """
        self.sf_dict = dict()
        sfid = ""
        pattern = r"^\d+\.\d+\.\d+\.\d+"
        try:
            with open(inputfile, "r") as f:
                for line in f:
                    line = line.strip("\n").split("    ")
                    if re.search(pattern, line[0]):
                        # print(line)
                        sfid = f"G3DSA:{line[0]}"

                        if sfid in self.unintegrated:
                            tmp_name = line[2].split(":")[1]
                            if tmp_name != "":
                                self.unintegrated[sfid]["name"] = tmp_name
                            else:
                                del self.unintegrated[sfid]
        except FileNotFoundError as e:
            print(e)
            sys.exit()

    def printNewNames(self, outputfile):
        """
        print results in file
        """
        with open(outputfile, "w") as f:
            f.write(
                "Status,Gene3D,protein count (03/21),IPR,comment,Method name in InterPro,CATH method name\n"
            )
            for sfid in self.unintegrated:
                count = 0
                newname = self.unintegrated[sfid]["name"].replace(",", ";")

                if sfid in self.protein_counts:
                    count = self.protein_counts[sfid]
                if sfid in self.descriptions:
                    oldname = self.descriptions[sfid].replace(",", ";")

                if self.unintegrated[sfid]["entry"] != None:
                    entry = self.unintegrated[sfid]["entry"]
                else:
                    entry = ""

                if len(self.unintegrated[sfid]["comment"]) > 0:
                    comments = " | ".join(self.unintegrated[sfid]["comment"])
                else:
                    comments = ""

                f.write(f",{sfid},{count},{entry},{comments},{oldname},{newname}\n")

        print(f"Output file location: {outputfile}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--user", help="username for database connection", required=True)
    parser.add_argument("-p", "--password", help="password for database connection", required=True)
    parser.add_argument("-s", "--schema", help="database schema to connect to", required=True)

    parser.add_argument("-f", "--inputfile", help="CATHb file containing new data", required=True)

    parser.add_argument("-o", "--outputfile", help="output file containing new data", required=True)

    parser.add_argument("-c", "--commentfile", help="File containing previous comments")

    args = parser.parse_args()

    cathb = new_cathb(args.user, args.password, args.schema)

    if args.commentfile:
        cathb.getComments(args.commentfile)

    cathb.getUnintegrated()
    cathb.getSfInfo(args.inputfile)
    cathb.printNewNames(args.outputfile)
