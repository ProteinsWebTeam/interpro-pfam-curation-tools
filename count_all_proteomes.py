#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script counts newly created InterPro entries between 2 dates from proteomes of key species (InterPro BBR grant)

@arguments [-u USER]: database user
           [-p PASSWORD]: database password for the user
           [-s SCHEMA]: database schema to use
           [-b OLD_RELEASE_VERSION]: InterPro version to start the comparison
           [-e NEW_RELEASE_VERSION]: InterPro version to stop the comparison
           [-d PROTEIN_DIR]: Directory to download InterPro release files

@usage: bsub -M 40000 -R"rusage[mem=40000]" python count_all_proteomes.py 
"""


import argparse

from count_integrated_proteome import statistics

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--user", help="username for database connection", required=True)
    parser.add_argument("-p", "--password", help="password for database connection", required=True)
    parser.add_argument("-s", "--schema", help="database schema to connect to", required=True)
    parser.add_argument("-b", "--old_version", help="InterPro version to compare to", required=True)
    parser.add_argument("-e", "--new_version", help="Current InterPro release version")
    parser.add_argument("-d", "--proteindir", help="Protein file directory", required=True)
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

    # initialising
    stats = statistics()

    # getting InterPro files
    stats.old_proteinlist = stats.get_protein2ipr_file(args.old_version, args.proteindir)
    stats.new_proteinlist = stats.get_protein2ipr_file(args.new_version, args.proteindir)

    print("Searching integrated proteomes")

    # database connection
    stats.getConnection(args.user, args.password, args.schema)

    count_integrated_all = dict()

    for taxid in list_taxid:
        # print(f"Getting data for {list_taxid[taxid]} ({taxid})")
        stats.tax_id = taxid
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

    print(
        f"Newly created InterPro entries between InterPro {args.old_version} and {args.new_version}:"
    )
    for taxid in count_integrated_all:
        print(
            f"{list_taxid[taxid]}\t{count_integrated_all[taxid]['nb_integrated']} (+{count_integrated_all[taxid]['percent_integrated']}%)\t total integrated: {count_integrated_all[taxid]['total_integrated']} of {count_integrated_all[taxid]['protein_count']} ({count_integrated_all[taxid]['all_percent_integrated']}%)"
        )

    stats.close_connection()
