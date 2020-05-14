#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script counts newly created InterPro entries between 2 dates from proteomes of key species (InterPro BBR grant)

@arguments [-u USER]: database user
           [-p PASSWORD]: database password for the user
           [-s SCHEMA]: database schema to use
           [-b BEGIN_DATE]: date of creation to start count from
           [-e END_DATE]: date of creation to stop counting [optional]
"""


import argparse
import sys
from datetime import date

from count_integrated_proteome import statistics

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-u", "--user", help="username for database connection", required=True
    )
    parser.add_argument(
        "-p", "--password", help="password for database connection", required=True
    )
    parser.add_argument(
        "-s", "--schema", help="database schema to connect to", required=True
    )
    parser.add_argument(
        "-b", "--begin_date", help="start date for the search ", required=True
    )
    parser.add_argument("-e", "--end_date", help="end date for the search ")
    args = parser.parse_args()

    list_taxid = {
        161934: "Beta vulgaris",
        9913: "Bos taurus",
        9031: "Gallus gallus",
        183674: "Miscanthus X giganteus",
        8030: "Salmo salar",
        9823: "Sus scrofa",
        4565: "Triticum aestivum",
        4577: "Zea mays",
    }

    # initialising
    stats = statistics()

    end_date = args.end_date if args.end_date else str(date.today())

    if stats.verify_date_format(end_date) and stats.verify_date_format(args.begin_date):

        print("Searching integrated proteomes")

        # database connection
        stats.getConnection(args.user, args.password, args.schema)

        count_integrated = dict()

        for taxid in list_taxid:
            stats.tax_id = taxid
            count_integrated[taxid] = len(stats.count(args.begin_date, end_date))

        print(
            f"Newly created InterPro entries between ${args.begin_date}, ${end_date}:"
        )
        for taxid in count_integrated:
            print(f"{list_taxid[taxid]}\t{count_integrated[taxid]}")

        stats.close_connection()
    else:
        print("Error: wrong date format, usage YYYY-MM-DD")
        sys.exit()
