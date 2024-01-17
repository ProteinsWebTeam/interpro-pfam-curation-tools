#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script counts newly created InterPro entries between 2 InterPro release versions from proteomes or taxid of key species (InterPro BBR grant)

@arguments [CONFIG_FILE]: file containing database connection and files location (see config.ini)
           [-r taxid | proteome]: Specify if search for a taxid or a reference proteome

@usage example: bsub -M 40000 -R"rusage[mem=40000]" python count_all_proteomes.py config.ini -r taxid | proteome
"""


import argparse
import sys, os
from configparser import ConfigParser

from count_integrated_proteome import statistics

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="FILE", help="configuration file")
    parser.add_argument(
        "-r",
        "--ref_type",
        help="Search for taxid or Reference proteome",
        choices=["taxid", "proteome"],
        required=True,
    )
    args = parser.parse_args()

    list_taxid = {
        161934: "Beta vulgaris",
        9913: "Bos taurus",
        9031: "Gallus gallus",
        183674: "Miscanthus x giganteus",
        8030: "Salmo salar",
        9823: "Sus scrofa",
        4565: "Triticum aestivum",
        4577: "Zea mays",
        3702: "Arabidopsis thaliana",
    }

    ref_proteomes = {
        161934: "UP000035740",
        9913: "UP000009136",
        9031: "UP000000539",
        9823: "UP000008227",
        8030: "UP000087266",
        4577: "UP000007305",
        4565: "UP000019116",
        183674: "None",
        3702: "UP000006548",
    }

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    # initialising
    stats = statistics()
    stats.dir = config["dir"]["proteindir"]

    # verifying all database parameters are provided
    if args.ref_type == "taxid":
        try:
            stats.getConnection(config["database"]["user"], config["database"]["password"], config["database"]["schema"])
        except Exception as e:
            print(f"Error {e}, invalid database connection credentials")
            sys.exit()

    # getting InterPro files
    stats.old_proteinlist = stats.get_protein2ipr_file(config["interpro"]["old_version"])
    stats.new_proteinlist = stats.get_protein2ipr_file(config["interpro"]["new_version"])

    print("Searching integrated proteins")

    count_integrated_all = dict()

    for taxid, proteome in ref_proteomes.items():
        if args.ref_type == "taxid":
            stats.tax_id = taxid
            print(f"searching integrated taxid {taxid}")
        else:
            if proteome != "None":
                print(taxid, proteome)
                stats.proteome = proteome
            else:
                continue

        count_integrated, protein_count, total_integrated = stats.count()
        percent_integrated = round(count_integrated * 100 / protein_count, 2)
        all_percent_integrated = round(total_integrated * 100 / protein_count, 2)

        count_integrated_all[taxid] = {
            "nb_integrated": count_integrated,
            "percent_integrated": percent_integrated,
            "total_integrated": total_integrated,
            "protein_count": protein_count,
            "all_percent_integrated": all_percent_integrated,
        }

    try:
        stats.close_connection()
    except:
        pass

    # print(
    #     f"Newly created InterPro entries between InterPro {args.old_version} and {args.new_version}:"
    # )
    for taxid in count_integrated_all:
        print(
            f"{list_taxid[taxid]}\t{count_integrated_all[taxid]['nb_integrated']} (+{count_integrated_all[taxid]['percent_integrated']}%)\t total integrated: {count_integrated_all[taxid]['total_integrated']} of {count_integrated_all[taxid]['protein_count']} ({count_integrated_all[taxid]['all_percent_integrated']}%)"
        )

