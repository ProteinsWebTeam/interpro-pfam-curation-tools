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
        self.user, self.password = user, password
        self.schema = schema
        self.connection = self.getConnection()
        self.cursor = self.connection.cursor()
        self.unintegrated = self.getUnIntegratedInEntry()
        self.descriptions = self.getOldNames()
        self.protein_counts = self.getProteinCounts()
        self.integrated = self.getIntegrated()
        self.comments = dict()

    def getConnection(self):
        """
        Set database connection
        """

        connectString = "".join([self.user, "/", self.password, "@", self.schema])
        try:
            return cx_Oracle.connect(connectString)
        except:
            stackTrace = traceback.format_exc()
            print(stackTrace)
            if "invalid username" in stackTrace:
                print("Could not connect to {0} as user {1}".format(self.schema, self.user))
                print(
                    "NB if your oracle username contains the '$' character either escape it or surround it with quotes"
                )
                print("eg 'ops$craigm' or ops\$craigm")
                print("Otherwise the shell will remove the '$' and all subsequent characters!")
            sys.exit(1)

    def getOldNames(self):
        """
        Search for previous names given by CATH in the database for each method
        """
        descriptions = dict()
        request = """SELECT m.METHOD_AC, m.DESCRIPTION
                    FROM INTERPRO.METHOD m 
                    WHERE m.DESCRIPTION is not null and m.METHOD_AC like 'G3DSA%'
                """

        self.cursor.execute(request)
        results = self.cursor.fetchall()

        for row in results:
            descriptions[row[0]] = row[1]

        return descriptions

    def getProteinCounts(self):
        """
        Search for the number of proteins for each method
        """
        protein_counts = dict()
        request = """SELECT mm.METHOD_AC, mm.PROTEIN_COUNT
                    FROM INTERPRO.MV_METHOD_MATCH mm
                    WHERE mm.METHOD_AC like 'G3DSA%'
                """

        self.cursor.execute(request)
        results = self.cursor.fetchall()

        for row in results:
            protein_counts[row[0]] = row[1]

        return protein_counts

    def getComments(self, commentfile):
        """
        Search for the curation comments previously assigned to each method
        """
        pattern = r"([\d\.]+),(.+)"
        try:
            with open(commentfile, "r") as f:
                for line in f:
                    line = line.strip("\n")
                    m = re.search(pattern, line)
                    sf = m.group(1)
                    sfid = f"G3DSA:{sf}"
                    self.comments[sfid] = m.group(2)
        except FileNotFoundError as e:
            print(e)
            sys.exit()

    def getUnIntegratedInEntry(self):
        """
        Search for integrated methods in unchecked InterPro entries with the curators comments
        """
        unintegrated = dict()
        request = """SELECT e.ENTRY_AC, e2m.METHOD_AC, LISTAGG(mc.VALUE, '; ') WITHIN GROUP (ORDER BY mc.VALUE) COMMENTS
        FROM INTERPRO.ENTRY2METHOD e2m
        LEFT JOIN INTERPRO.ENTRY e on e2m.ENTRY_AC=e.ENTRY_AC
        LEFT JOIN INTERPRO.METHOD_COMMENT mc on e2m.METHOD_AC=mc.METHOD_AC
        WHERE e2m.METHOD_AC like 'G3DSA%' AND e.CHECKED='N'
        GROUP BY e.ENTRY_AC, e2m.METHOD_AC
        """

        self.cursor.execute(request)
        results = self.cursor.fetchall()

        for row in results:
            if row[2] == "N":
                self.unintegrated[row[1]] = {"comment": row[3], "entry": row[0]}

        return unintegrated

    def getIntegrated(self):
        """
        Search for integrated methods
        """

        request = """SELECT e2m.METHOD_AC 
                    FROM INTERPRO.ENTRY2METHOD e2m
                    JOIN INTERPRO.ENTRY e on e2m.ENTRY_AC=e.ENTRY_AC
                    WHERE e.CHECKED='Y'
                """
        self.cursor.execute(request)
        results = self.cursor.fetchall()

        return [row[0] for row in results]

    def getSfInfo(self, inputfile):
        """
        Get CATH new names for unintegrated methods
        """
        self.sf_dict = dict()
        sfid = ""
        try:
            with open(inputfile, "r") as f:
                for line in f:
                    line = line.strip("\n").split("     ")
                    if "ID" in line[0]:
                        sfid = f"G3DSA:{line[1].strip()}"
                        if sfid in self.unintegrated:
                            self.sf_dict[sfid] = {
                                "entry": self.unintegrated[sfid]["entry"],
                                "comment": self.unintegrated[sfid]["comment"],
                            }
                        else:
                            self.sf_dict[sfid] = {}

                    elif "NAME" in line[0]:
                        if line[1] != "-":
                            self.sf_dict[sfid]["name"] = line[1]
                        else:
                            del self.sf_dict[sfid]
        except FileNotFoundError as e:
            print(e)
            sys.exit()

    def printNewNames(self, outputfile):
        """
        print results in file
        """
        with open(outputfile, "w") as f:
            f.write(
                "Process,Gene3D,protein count (01/20),InterPro method name,CATH method name,comment,IPR\n"
            )
            for sfid in self.sf_dict:
                string_to_print = ""
                oldname = ""
                count = 0
                if sfid not in self.integrated:
                    if sfid in self.protein_counts:
                        count = self.protein_counts[sfid]
                    if sfid in self.descriptions:
                        oldname = self.descriptions[sfid]
                    if oldname != self.sf_dict[sfid]["name"]:
                        string_to_print = f",{sfid},{count},{oldname.replace(',',';')},{self.sf_dict[sfid]['name'].replace(',',';')}"

                        if "entry" in self.sf_dict[sfid]:
                            string_to_print += (
                                ",self.sf_dict[sfid]['entry'],N,self.sf_dict[sfid]['comment'],"
                            )
                        else:
                            string_to_print += ",,,,"
                        if self.comments and sfid in self.comments:
                            string_to_print += self.comments[sfid]
                        string_to_print += "\n"

                    f.write(string_to_print)

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

    cathb.getSfInfo(args.inputfile)
    cathb.printNewNames(args.outputfile)
