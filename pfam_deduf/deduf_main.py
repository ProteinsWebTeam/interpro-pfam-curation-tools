#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script searches DUF Pfam and tries to find evidence of characterisation

@arguments [-u USER]: database user
           [-p PASSWORD]: database password for the user
           [-s SCHEMA]: database schema to use
           [-o OPTION]: 1)search swissprot names for unintegrated pfam DUF  
                        2)search for DUF in literature
"""
import sys
import argparse
import os
from search_swissprot import pfam_swiss
from search_literature import pfam_litterature

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
        default=[1, 2],
    )
    parser.add_argument(
        "-d",
        "--outputdir",
        help="output directory to write the data",
        required=False,
        default="/nfs/production/interpro/programmers/typhaine/interpro-curation-tools/results/duf",
    )

    args = parser.parse_args()

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

        pfam_pip.save_duf_list_in_file(args.outputdir)
        pfam_pip.save_swissprot_in_file(args.outputdir)

    elif args.option == "2":
        print("Searching DUF in literature")
        pfam_pip = pfam_litterature()
        pfam_pip.getConnection(args.user, args.password, args.schema)
        pfam_pip.get_pfam_duf_list()

        pfam_pip.save_duf_list_in_file(args.outputdir)
        for pfamid in pfam_pip.list_duf:
            print(pfamid)
            pfam_pip.get_list_articles_duf(pfamid)
            pfam_pip.get_list_articles_pfamid(pfamid)

            # print(pfam_pip.list_duf[pfamid])
            # break

        outputfile = os.path.join(args.outputdir, "duf_litterature.csv")
        pfam_pip.save_result_in_file(outputfile)

