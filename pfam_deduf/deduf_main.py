#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script searches DUF Pfam and tries to find evidence of characterisation

@arguments [CONFIG_FILE]: file containing database connection and files location (see config_pfam.ini)
           [-o OPTION]: 1)search swissprot names for unintegrated pfam DUF  
                        2)search for DUF in literature
                        3)search GO term/keywords in Swissprot matches
                        4)search predicted structures

"""

import argparse
import os
from configparser import ConfigParser

from search_swissprot import pfam_swiss
from search_literature import pfam_litterature
from search_go_keywords import pfam_go
from search_models import pfam_model

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")

    parser.add_argument(
        "-o",
        "--option",
        help="specify option to run 1)search swissprot names for unintegrated pfam DUF  2)search for DUF in literature 3)search GO terms/keywords in SwissProt 4)search predicted structures",
        required=True,
        default=[1, 2, 3, 4],
    )

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    outputdir = config["dir"]["outputdir"]
    outfile = os.path.join(outputdir, "list_duf.csv")
    outfileipr = os.path.join(outputdir, "list_duf_with_ipr_name.csv")
    outputswiss = os.path.join(outputdir, "duf_swiss_names.tsv")
    model_dir = config["dir"]["modeldir"]

    # database connection values
    user = config["database"]["user"]
    password = config["database"]["password"]
    schema = config["database"]["schema"]

    print(args.option)
    if args.option == "1":
        print("Searching for SwissProt names")
        # initialising
        pfam_pip = pfam_swiss()
        pfam_pip.getConnection(user, password, schema)

        pfam_pip.get_pfam_duf_list()
        # pfam_pip.get_pfam_duf_list_unintegrated()
        pfam_pip.close_connection()

        pfam_pip.get_list_swissprot_names()
        pfam_pip.get_nb_domain_per_interpro()
        pfam_pip.get_go_annotation()

        pfam_pip.save_duf_list_in_file(outfile, outfileipr)

        pfam_pip.save_swissprot_in_file(outputswiss)

    elif args.option == "2":
        print("Searching DUF in literature")
        pfam_pip = pfam_litterature()
        pfam_pip.getConnection(user, password, schema)
        pfam_pip.get_pfam_duf_list()
        pfam_pip.close_connection()

        pfam_pip.save_duf_list_in_file(outfile, outfileipr)

        for pfamid in pfam_pip.list_duf:
            print(pfamid)
            pfam_pip.get_list_articles_duf(pfamid)
            pfam_pip.get_list_articles_pfamid(pfamid)

        outputfile = os.path.join(outputdir, "duf_litterature.csv")
        pfam_pip.save_result_in_file(outputfile)

    elif args.option == "3":
        print("Searching for GO and Keywords annotations in SwissProt entries")
        pfam_pip = pfam_go()
        try:
            pfam_pip.get_list_duf_from_file(outputswiss)
            pfam_pip.get_swissprot_accessions()
            pfam_pip.get_go_annotation()

            outputfile = os.path.join(outputdir, "duf_go_keywords.tsv")
            pfam_pip.save_results_in_file(outputfile)

        except FileNotFoundError:
            print(f"Error file {outputswiss} not found, run option 1 first")

    elif args.option == "4":
        print("Searching for predicted structures ")
        pfam_pip = pfam_model(model_dir)
        pfam_pip.getConnection(user, password, schema)
        pfam_pip.get_pfam_duf_list()
        pfam_pip.close_connection()

        pfam_pip.search_model_families()

