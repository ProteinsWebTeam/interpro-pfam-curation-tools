#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script searches DUF Pfam and tries to find evidence of characterisation

@arguments [-u USER]: database user
           [-p PASSWORD]: database password for the user
           [-s SCHEMA]: database schema to use
           [-o OPTION]: 1)search swissprot names for unintegrated pfam DUF  
                        2)search for DUF in literature
                        3)search GO term/keywords in Swissprot matches
           [-f INPUTFILE]: if -o 3 a file containing a subset of identifiers and their InterPro matches can be imported, by default the file result of -o 1 is used
"""

import argparse
import os
from search_swissprot import pfam_swiss
from search_literature import pfam_litterature
from search_go_keywords import pfam_go

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--user", help="username for database connection", required=True)
    parser.add_argument("-p", "--password", help="password for database connection", required=True)
    parser.add_argument("-s", "--schema", help="database schema to connect to", required=True)
    parser.add_argument(
        "-o",
        "--option",
        help="specify option to run 1)search swissprot names for unintegrated pfam DUF  2)search for DUF in literature",
        required=True,
        default=[1, 2, 3],
    )
    parser.add_argument(
        "-d",
        "--outputdir",
        help="output directory to write the data",
        required=False,
        default="/nfs/production/interpro/programmers/typhaine/interpro-curation-tools/results/duf",
    )
    parser.add_argument(
        "-f",
        "--input_file",
        help="tsv file containing DUF/IPR/SwissProt matches",
        required=False,
        default="/nfs/production/interpro/programmers/typhaine/interpro-curation-tools/results/duf/duf_swiss_names.tsv",
    )

    args = parser.parse_args()

    outfile = os.path.join(args.outputdir, "list_duf.csv")
    outfileipr = os.path.join(args.outputdir, "list_duf_with_ipr_name.csv")

    print(args.option)
    if args.option == "1":
        print("Searching for SwissProt names")
        # initialising
        pfam_pip = pfam_swiss()
        pfam_pip.getConnection(args.user, args.password, args.schema)

        pfam_pip.get_pfam_duf_list()
        # pfam_pip.get_pfam_duf_list_unintegrated()

        pfam_pip.get_list_swissprot_names()
        pfam_pip.get_nb_domain_per_interpro()
        pfam_pip.get_go_annotation()

        pfam_pip.save_duf_list_in_file(outfile, outfileipr)

        outputfile = os.path.join(args.outputdir, "duf_swiss_names.tsv")
        pfam_pip.save_swissprot_in_file(args.outputdir)

    elif args.option == "2":
        print("Searching DUF in literature")
        pfam_pip = pfam_litterature()
        pfam_pip.getConnection(args.user, args.password, args.schema)
        pfam_pip.get_pfam_duf_list()

        pfam_pip.save_duf_list_in_file(outfile, outfileipr)

        for pfamid in pfam_pip.list_duf:
            print(pfamid)
            pfam_pip.get_list_articles_duf(pfamid)
            pfam_pip.get_list_articles_pfamid(pfamid)

        outputfile = os.path.join(args.outputdir, "duf_litterature.csv")
        pfam_pip.save_result_in_file(outputfile)

    elif args.option == "3":
        print("Searching for GO and Keywords annotations in SwissProt entries")
        pfam_pip = pfam_go()
        try:
            pfam_pip.get_list_duf_from_file(args.input_file)
            pfam_pip.get_swissprot_accessions()
            pfam_pip.get_go_annotation()

            outputfile = os.path.join(args.outputdir, "duf_go_keywords.tsv")
            pfam_pip.save_results_in_file(outputfile)

        except FileNotFoundError:
            print(f"Error file {args.input_file} not found")

