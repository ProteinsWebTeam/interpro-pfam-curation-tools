"""
@author T. Paysan-Lafosse

@brief This script counts newly created InterPro entries between 2 InterPro releases for a given reference proteome, taxid or Scientific name (InterPro BBR grant)

@arguments [CONFIG_FILE]: file containing database connection and files location (see config.ini)
           [-r REF_PROTEOME | -t TAXID | -o ORGANISM] : UniProt reference proteome identifier, taxid or Scientific name 

@usage example: bsub -M 40000 -R"rusage[mem=40000]" python count_all_proteomes.py config.ini -o UPXXXX
"""

import argparse
import sys, os
from configparser import ConfigParser
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
    parser.add_argument("config", metavar="FILE", help="configuration file")
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

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    # initialising
    stats = statistics()
    stats.dir = config["dir"]["proteindir"]

    # initialise tax_id value
    if args.organism:
        try:
            stats.getConnection(config["database"]["user"], config["database"]["password"], config["database"]["schema"])
            print(f"Searching taxid for {args.organism}")
            stats.search_taxid(args.organism)
        except Exception as e:
            print(f"Error {e}, invalid database connection credentials")
            sys.exit()
    elif args.taxid:
        try:
            stats.getConnection(config["database"]["user"], config["database"]["password"], config["database"]["schema"])
            stats.tax_id = args.taxid
        except Exception as e:
            print(f"Error {e}, invalid database connection credentials")
            sys.exit()
    else:
        stats.proteome = args.proteome

    # getting InterPro files
    stats.old_proteinlist = stats.get_protein2ipr_file(config["interpro"]["old_version"])
    stats.new_proteinlist = stats.get_protein2ipr_file(config["interpro"]["new_version"])

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

    try:
        stats.close_connection()
    except:
        pass

