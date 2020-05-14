"""
@author T. Paysan-Lafosse

@brief This script counts newly created InterPro entries between 2 dates for a given organism or taxid (InterPro BBR grant)

@arguments [-u USER]: database user
           [-p PASSWORD]: database password for the user
           [-s SCHEMA]: database schema to use
           [-o ORGANISM or -t TAXID]: organism (scientific name) or taxid to look for
           [-b BEGIN_DATE]: date of creation to start count from
           [-e END_DATE]: date of creation to stop counting
"""

import argparse
import sys
import re
from datetime import date
from utils import proteome


class statistics(proteome):
    def get_integrated(self, protein_list, start_date, end_date):
        """
        Search integrated proteins between 2 dates

        Args:
            protein_list: list containing proteins to search for
            start_date: lower limit for the search
            end_date: upper limit for the search

        Yields:
            list_integrated: list of proteins integrated in InterPro entries

        """

        # print("Searching integrated proteins")
        uniprot_chunks = list(self.chunks(protein_list, 1000))

        list_integrated = set()
        for chunk in uniprot_chunks:
            protein_list_quote = [f"'{row}'" for row in chunk]
            request = f"SELECT UNIQUE E.ENTRY_AC \
                    FROM INTERPRO.MV_ENTRY2PROTEIN E2P \
                    JOIN INTERPRO.ENTRY E ON E.ENTRY_AC=E2P.ENTRY_AC \
                    WHERE E2P.PROTEIN_AC IN ({','.join(protein_list_quote)}) \
                    AND E.CHECKED='Y' AND E.CREATED>=TO_DATE('{start_date}', 'YYYY-MM-DD') \
                    AND E.CREATED<=TO_DATE('{end_date}', 'YYYY-MM-DD')"
            self.cursor.execute(request)

            list_integrated.update(set([row[0] for row in self.cursor]))
        return list_integrated

    def verify_date_format(self, date):
        pattern = r"\d\d\d\d-\d\d-\d\d"
        if re.match(pattern, date):
            return True
        else:
            print("Wrong date format, usage: yyyy-mm-dd")
            sys.exit()

    def count(self, begin_date, end_date):

        # search the proteome
        protein_list = self.get_proteins()

        list_integrated = self.get_integrated(protein_list, begin_date, end_date)

        return list_integrated


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--user", help="username for database connection", required=True)
    parser.add_argument("-p", "--password", help="password for database connection", required=True)
    parser.add_argument("-s", "--schema", help="database schema to connect to", required=True)
    parser.add_argument("-b", "--begin_date", help="start date for the search ", required=True)
    parser.add_argument("-e", "--end_date", help="end date for the search ")
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-o", "--organism", help="Scientific name of the organism to get the conservation score for"
    )
    group.add_argument(
        "-t", "--taxid", help="Taxid of the organism to get the conservation score for"
    )

    args = parser.parse_args()

    if args.organism:
        taxid = args.organism
    elif args.tax_id:
        taxid = args.tax_id
    else:
        print("Error no organism or taxid provided")
        sys.exit(1)

    # initialising
    stats = statistics()

    end_date = args.end_date if args.end_date else str(date.today())

    if stats.verify_date_format(end_date) and stats.verify_date_format(args.begin_date):

        stats.getConnection(args.user, args.password, args.schema)

        # initialise tax_id value
        pattern = r"\d+"
        if re.search(pattern, taxid):
            stats.tax_id = taxid
        else:
            print(f"Searching taxid for {taxid}")
            stats.search_taxid(taxid)

        # search for integrated proteins
        list_integrated = stats.count(args.begin_date, end_date)

        print(
            f"Newly created InterPro entries for {stats.tax_id} between {args.begin_date}, {end_date}: {len(list_integrated)}"
        )

        stats.close_connection()
    else:
        print("Error: wrong date format, usage YYYY-MM-DD")
        sys.exit()

