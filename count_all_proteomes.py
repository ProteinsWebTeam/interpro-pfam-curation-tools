#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script counts newly created InterPro entries between 2 InterPro release versions from proteomes or taxid of key species (InterPro BBR grant)

@arguments [-u USER]: database user (if -r taxid)
           [-p PASSWORD]: database password for the user (if -r taxid)
           [-s SCHEMA]: database schema to use (if -r taxid)
           [-b OLD_RELEASE_VERSION]: InterPro version to start the comparison
           [-e NEW_RELEASE_VERSION]: InterPro version to stop the comparison
           [-d PROTEIN_DIR]: Directory to download InterPro release files
           [-r taxid | proteome]: Specify if search for a taxid or a reference proteome

@usage: bsub -M 40000 -R"rusage[mem=40000]" python count_all_proteomes.py -b 76.0 -e 81.0 -d data_folder
"""


import argparse
import sys

from count_integrated_proteome import statistics

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--user", help="username for database connection")
    parser.add_argument("-p", "--password", help="password for database connection")
    parser.add_argument("-s", "--schema", help="database schema to connect to")
    parser.add_argument("-b", "--old_version", help="InterPro version to compare to", required=True)
    parser.add_argument("-e", "--new_version", help="Current InterPro release version")
    parser.add_argument("-d", "--proteindir", help="Protein file directory", required=True)
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
    }

    # initialising
    stats = statistics()
    stats.dir = args.proteindir

    # verifying all database parameters are provided
    if args.ref_type == "taxid":
        if args.user and args.schema and args.password:
            stats.getConnection(args.user, args.password, args.schema)
        else:
            print("Error, missing database connection credentials -u user -p password -s schema")
            sys.exit()

    # getting InterPro files
    stats.old_proteinlist = stats.get_protein2ipr_file(args.old_version)
    stats.new_proteinlist = stats.get_protein2ipr_file(args.new_version)

    print("Searching integrated proteomes")

    count_integrated_all = dict()

    for taxid, proteome in ref_proteomes.items():
        if args.ref_type == "taxid":
            stats.tax_id = taxid
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

