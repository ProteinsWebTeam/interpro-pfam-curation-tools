"""
@author T. Paysan-Lafosse

@brief This script counts newly created InterPro entries between 2 InterPro releases for a given reference proteome (InterPro BBR grant)

@arguments [-b OLD_RELEASE_VERSION]: InterPro version to start the comparison
           [-e NEW_RELEASE_VERSION]: InterPro version to stop the comparison
           [-d PROTEIN_DIR]: Directory to download InterPro release files
           [-o PROTEOME]: UniProt reference proteome identifier

@usage: bsub -M 40000 -R"rusage[mem=40000]" python count_all_proteomes.py -b 76.0 -e 81.0 -d protein_folder -o UPXXXX
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
        """
		Initialise list of proteins and count newly integrated

		Args: None

		Yields: 
			count_integrated: number of newly integrated proteins
            len(taxon_proteinlist): total number of proteins
            total_integrated: total number of proteins integrated

		"""

        # search the proteome
        # taxon_proteinlist = self.get_proteins()
        taxon_proteinlist = self.get_proteins_from_uniprot()

        count_integrated, total_integrated = self.get_integrated(
            taxon_proteinlist, self.old_proteinlist, self.new_proteinlist
        )

        return count_integrated, len(taxon_proteinlist), total_integrated


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--old_version", help="InterPro version to compare to", required=True)
    parser.add_argument("-e", "--new_version", help="Current InterPro release version")
    parser.add_argument("-d", "--proteindir", help="Protein file directory", required=True)
    parser.add_argument(
        "-o", "--proteome", help="Reference proteome to search the information for", required=True
    )

    args = parser.parse_args()

    # initialising
    stats = statistics()
    stats.dir = args.proteindir
    stats.proteome = args.proteome

    # getting InterPro files
    stats.old_proteinlist = stats.get_protein2ipr_file(args.old_version)
    stats.new_proteinlist = stats.get_protein2ipr_file(args.new_version)

    # search for integrated proteins
    count_integrated, protein_count, total_integrated = stats.count()
    percent_integrated = round(count_integrated * 100 / protein_count, 2)

    print(
        f"Newly created InterPro entries for {stats.proteome} between InterPro {args.old_version} and InterPro {args.new_version}: {count_integrated} (+{percent_integrated}%), total integrated: {total_integrated} of {protein_count}"
    )

