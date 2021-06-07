"""
@author T. Paysan-Lafosse

@brief This script counts newly created InterPro entries between 2 InterPro releases for a given reference proteome, taxid or Scientific name (InterPro BBR grant)

@arguments [-u USER]: database user (if -t or -o)
           [-p PASSWORD]: database password for the user (if -t or -o)
           [-s SCHEMA]: database schema to use (if -t or -o)
           [-b OLD_RELEASE_VERSION]: InterPro version to start the comparison
           [-e NEW_RELEASE_VERSION]: InterPro version to stop the comparison
           [-d PROTEIN_DIR]: Directory to download InterPro release files
           [-r REF_PROTEOME | -t TAXID | -o ORGANISM] : UniProt reference proteome identifier, taxid or Scientific name 

@usage: bsub -M 40000 -R"rusage[mem=40000]" python count_all_proteomes.py -b 76.0 -e 81.0 -d protein_folder -o UPXXXX
"""

import argparse
import sys
from utils import proteome


class statistics(proteome):
    def get_integrated(self, taxon_proteinlist):

        """ Description:  Search integrated proteins between 2 InterPro releases

        Args:
            taxon_proteinlist: list containing proteins for a particular taxon/proteome to search for

        Yields:
            count_integrated: number proteins newly integrated in InterPro entries
            total_integrated: total number of proteins integrated in InterPro entries

        """

        integrated_new = 0
        integrated_old = 0
        # search if protein for taxon/proteome is integrated
        for protein in taxon_proteinlist:
            if protein in self.old_proteinlist:  # old interpro version
                integrated_old += 1
            if protein in self.new_proteinlist:  # new interpro version
                integrated_new += 1

        return integrated_new - integrated_old, integrated_new

    def count(self):
        """
		Initialise list of proteins and count newly integrated

		Args: None

		Yields: 
			count_integrated: number of newly integrated proteins
            len(taxon_proteinlist): total number of proteins
            total_integrated: total number of proteins integrated

		"""

        # search proteins
        if self.tax_id:  # taxid
            taxon_proteinlist = self.get_proteins()
        else:  # proteome
            taxon_proteinlist = self.get_proteins_from_uniprot()

        count_integrated, total_integrated = self.get_integrated(taxon_proteinlist)

        return count_integrated, len(taxon_proteinlist), total_integrated


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--user", help="username for database connection")
    parser.add_argument("-p", "--password", help="password for database connection")
    parser.add_argument("-s", "--schema", help="database schema to connect to")
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-o", "--organism", help="Scientific name of the organism to get the conservation score for"
    )
    group.add_argument(
        "-t", "--taxid", help="Taxid of the organism to get the conservation score for"
    )
    group.add_argument(
        "-r", "--proteome", help="Reference proteome to search the information for", required=True
    )
    parser.add_argument("-b", "--old_version", help="InterPro version to compare to", required=True)
    parser.add_argument("-e", "--new_version", help="Current InterPro release version")
    parser.add_argument("-d", "--proteindir", help="Protein file directory", required=True)

    args = parser.parse_args()

    # initialising
    stats = statistics()
    stats.dir = args.proteindir

    # initialise tax_id value
    if args.organism:
        if args.user and args.schema and args.password:
            stats.getConnection(args.user, args.password, args.schema)
            print(f"Searching taxid for {args.organism}")
            stats.search_taxid(args.organism)
        else:
            print("Error, missing database connection credentials -u user -p password -s schema")
            sys.exit()
    elif args.taxid:
        if args.user and args.schema and args.password:
            stats.getConnection(args.user, args.password, args.schema)
            stats.tax_id = args.taxid
        else:
            print("Error, missing database connection credentials -u user -p password -s schema")
            sys.exit()
    else:
        stats.proteome = args.proteome

    # getting InterPro files
    stats.old_proteinlist = stats.get_protein2ipr_file(args.old_version)
    stats.new_proteinlist = stats.get_protein2ipr_file(args.new_version)

    # search for integrated proteins
    count_integrated, protein_count, total_integrated = stats.count()
    percent_integrated = round(count_integrated * 100 / protein_count, 2)

    if stats.proteome:
        print(
            f"Newly created InterPro entries for {stats.proteome} between InterPro {args.old_version} and InterPro {args.new_version}: {count_integrated} (+{percent_integrated}%), total integrated: {total_integrated} of {protein_count}"
        )
    else:
        print(
            f"Newly created InterPro entries for {stats.tax_id} between InterPro {args.old_version} and InterPro {args.new_version}: {count_integrated} (+{percent_integrated}%), total integrated: {total_integrated} of {protein_count}"
        )
        stats.close_connection()

