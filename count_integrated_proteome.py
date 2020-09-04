"""
@author T. Paysan-Lafosse

@brief This script counts newly created InterPro entries between 2 dates for a given organism or taxid (InterPro BBR grant)

@arguments [-u USER]: database user
           [-p PASSWORD]: database password for the user
           [-s SCHEMA]: database schema to use
           [-o ORGANISM or -t TAXID]: organism (scientific name) or taxid to look for
           [-b OLD_RELEASE_VERSION]: InterPro version to start the comparison
           [-e NEW_RELEASE_VERSION]: InterPro version to stop the comparison

@usage: bsub -M 40000 -R"rusage[mem=40000]" python count_all_proteomes.py
"""

import argparse
import sys
from utils import proteome


class statistics(proteome):
    def get_integrated(self, taxon_proteinlist, old_proteinlist, new_proteinlist):

        """ Description:  Search integrated proteins between 2 InterPro releases

        Args:
            taxon_proteinlist: list containing proteins to search for
            old_proteinlist: list of proteins integrated in old InterPro version
            end_date: list of proteins integrated in new InterPro version

        Yields:
            count_integrated: number proteins newly integrated in InterPro entries
            len(taxon_proteinlist): total number of proteins for the current taxon
            total_integrated: total number of proteins integrated in InterPro entries

        """

        integrated_new = 0
        integrated_old = 0
        for protein in taxon_proteinlist:
            if protein in old_proteinlist:
                integrated_old += 1
            if protein in new_proteinlist:
                integrated_new += 1

        return integrated_new - integrated_old, integrated_new

    def count(self):

        # search the proteome
        taxon_proteinlist = self.get_proteins()

        count_integrated, total_integrated = self.get_integrated(
            taxon_proteinlist, self.old_proteinlist, self.new_proteinlist
        )

        return count_integrated, len(taxon_proteinlist), total_integrated


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--user", help="username for database connection", required=True)
    parser.add_argument("-p", "--password", help="password for database connection", required=True)
    parser.add_argument("-s", "--schema", help="database schema to connect to", required=True)
    parser.add_argument("-b", "--old_version", help="InterPro version to compare to", required=True)
    parser.add_argument("-e", "--new_version", help="Current InterPro release version")
    parser.add_argument("-d", "--proteindir", help="Protein file directory", required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-o", "--organism", help="Scientific name of the organism to get the conservation score for"
    )
    group.add_argument(
        "-t", "--taxid", help="Taxid of the organism to get the conservation score for"
    )

    args = parser.parse_args()

    # initialising
    stats = statistics()

    # getting InterPro files
    stats.old_proteinlist = stats.get_protein2ipr_file(args.old_version, args.proteindir)
    stats.new_proteinlist = stats.get_protein2ipr_file(args.new_version, args.proteindir)

    stats.getConnection(args.user, args.password, args.schema)

    # initialise tax_id value
    if args.organism:
        print(f"Searching taxid for {args.organism}")
        stats.search_taxid(args.organism)
    elif args.taxid:
        stats.tax_id = args.taxid
    else:
        print("Error no organism or taxid provided")
        sys.exit(1)

    # search for integrated proteins
    count_integrated, protein_count, total_integrated = stats.count()
    percent_integrated = round(count_integrated * 100 / protein_count, 2)

    print(
        f"Newly created InterPro entries for {stats.tax_id} between InterPro {args.old_version} and InterPro {args.new_version}: {count_integrated} (+{percent_integrated}%), total integrated: {total_integrated} of {protein_count}"
    )

    stats.close_connection()

