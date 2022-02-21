#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script searches unintegrated proteins for a given organism or taxid 
        in InterPro signatures, if they are not found in signatures, they are clustered based on UniRef clusters

@arguments [CONFIG_FILE]: config file containing database credentials and output dir (see config.ini)
           [-o ORGANISM or -t TAXID]: organism (scientific name) or taxid to look for
"""

import argparse
import os, re, sys
import requests
from configparser import ConfigParser

from utils import proteome


class protein_pipeline(proteome):
    def __init__(self):
        super().__init__()
        self.uniref50 = dict()
        self.clusters = dict()

    def get_integrated(self, protein_list):
        """
        Search integrated proteins

        Args:
            protein_list: list containing proteins to search for

        Yields:
            list_integrated: list of proteins integrated in InterPro entries

        """

        print("Searching for integrated proteins")
        uniprot_chunks = self.chunks(list(protein_list), 1000)

        list_integrated = set()
        for chunk in uniprot_chunks:
            protein_list_quote = [f"'{row}'" for row in chunk]
            request = f"SELECT P.PROTEIN_AC \
                    FROM INTERPRO.MV_ENTRY2PROTEIN E2P \
                    JOIN INTERPRO.PROTEIN P ON E2P.PROTEIN_AC=P.PROTEIN_AC \
                    WHERE E2P.PROTEIN_AC IN ({','.join(protein_list_quote)}) \
                    AND P.FRAGMENT='N'"
            self.cursor.execute(request)

            list_integrated.update(set([row[0] for row in self.cursor]))
        return list_integrated

    def get_count_signature_taxid(self, list_signatures):
        """
            Search for protein counts for a list of InterPro signatures

        Args:
            list_signatures: list of InterPro signatures

        Yields:
            count_prot_signatures: dictionnary with signature as key and protein_count as value

        """
        count_prot_signatures = dict()

        signature_chunks = list(self.chunks(list(list_signatures), 1000))
        for chunk in signature_chunks:
            signature_list_quote = [f"'{row}'" for row in chunk]
            request = f"SELECT M2P.METHOD_AC,COUNT(P.PROTEIN_AC) \
                FROM INTERPRO.PROTEIN P \
                JOIN INTERPRO.MV_METHOD2PROTEIN M2P ON P.PROTEIN_AC = M2P.PROTEIN_AC \
                JOIN INTERPRO.ETAXI ET ON P.TAX_ID = ET.TAX_ID \
                WHERE ET.TAX_ID=:1 \
                AND M2P.METHOD_AC IN ({','.join(signature_list_quote)}) \
                AND P.FRAGMENT='N' \
                GROUP BY M2P.METHOD_AC"

            self.cursor.execute(request, (self.tax_id,))
            count_prot_signatures.update({row[0]: row[1] for row in self.cursor})

        return count_prot_signatures

    def get_accession_in_signature(self, folder, protein_list):
        """
        Search for proteins found in InterPro signatures but not integrated
        Write the results in a csv file with each row corresponding to a protein/signature pair (protein,dbcode,organism,signature,total_prot_count,count_proteome,comment)

        Args:
            folder: output directory
            protein_list: list containing proteins to search for

        Yields:
            list of proteins found in unintegrated signatures

        """

        print("Searching for unintegrated proteins in signature")
        uniprot_chunks = list(self.chunks(list(protein_list), 1000))
        list_signatures_processed = set()
        # list_proteins_with_signature = dict()
        nbprot_in_signature = 0
        unintegrated_file = os.path.join(
            folder, f"unintegrated_prot_in_signatures_{self.tax_id}.csv"
        )
        with open(unintegrated_file, "w") as outf:
            outf.write("signature,total_prot_count,count_swiss_prot\n")

            for chunk in uniprot_chunks:
                # if comments needed in future:  C.VALUE, LISTAGG(MC.VALUE, '; ') WITHIN GROUP (ORDER BY MC.VALUE) COMMENTS
                protein_list_quote = [f"'{row}'" for row in chunk]
                request = f"SELECT M2P.METHOD_AC \
                            FROM INTERPRO.PROTEIN P \
                            JOIN INTERPRO.MV_METHOD2PROTEIN M2P ON P.PROTEIN_AC = M2P.PROTEIN_AC \
                            LEFT JOIN INTERPRO.ENTRY2METHOD E2M ON E2M.METHOD_AC = M2P.METHOD_AC \
                            WHERE P.PROTEIN_AC IN ({','.join(protein_list_quote)}) \
                            AND P.FRAGMENT='N' \
                            AND M2P.METHOD_AC not like '%:SF%'\
                            AND E2M.METHOD_AC is null"

                # print(request)
                self.cursor.execute(request)

                results = self.cursor.fetchall()
                nbprot_in_signature += len(results)
                for row in results:
                    signature = row[0]

                    if signature not in list_signatures_processed:
                        list_signatures_processed.add(signature)
                        # ony consider CDD and PANTHER family signatures
                        if re.search("cd", signature) or re.search("PTHR.+", signature):
                            if not self.has_comment(signature):
                                swissprot_count = self.get_swissprot_count(signature)

                                # print(signature, protein_count, swissprot_count)
                                if swissprot_count > 0:
                                    protein_count = self.get_protein_count(signature)
                                    outf.write(f"{signature},{protein_count},{swissprot_count}\n")
                                # list_proteins_with_signature[signature] = [
                                #     protein_count,
                                #     swissprot_count,
                                # ]

                # for signature, proteins in list_proteins_with_signature.items():
                #     if proteins[1] != 0:
                #         outf.write(f"{signature},{proteins[0]},{proteins[1]},\n")
        return nbprot_in_signature

    def has_comment(self, signature):
        request = f"select mc.value\
                    from interpro.method m \
                    left join interpro.method_comment mc on mc.method_ac=m.method_ac \
                    where m.method_ac='{signature}' and mc.status='Y'"

        self.cursor.execute(request)
        comments = []
        for row in self.cursor:
            comments.append(row[0])
        if len(comments) > 0:
            return True

        return False

    def get_protein_count(self, signature):
        request = f"SELECT COUNT(UNIQUE M.PROTEIN_AC) \
                    FROM INTERPRO.MATCH M \
                    JOIN INTERPRO.PROTEIN P ON M.PROTEIN_AC = P.PROTEIN_AC \
                    WHERE M.METHOD_AC='{signature}' \
                    AND P.FRAGMENT='N' "

        self.cursor.execute(request)
        return self.cursor.fetchone()[0]

    def get_swissprot_count(self, signature):
        request = f"SELECT COUNT(M.PROTEIN_AC) \
                FROM INTERPRO.MATCH M \
                INNER JOIN INTERPRO.PROTEIN P ON M.PROTEIN_AC = P.PROTEIN_AC \
                JOIN INTERPRO.MV_METHOD2PROTEIN M2P ON P.PROTEIN_AC = M2P.PROTEIN_AC\
                WHERE P.DBCODE = 'S' \
                AND P.FRAGMENT = 'N' \
                AND M.METHOD_AC = M2P.METHOD_AC \
                AND M2P.METHOD_AC='{signature}'"
        self.cursor.execute(request)
        return self.cursor.fetchone()[0]

    def search_uniprotid_in_uniref(self, uniprotid):
        """
        Search if the uniprotid is already referenced in the uniref50 dictionnary to avoid querying UniProt multiple times

        Args:
            uniprotid: UniProt accession to search for

        Yields:
            uniref: UniRef cluster found
            False: uniprotid not found

        """

        for uniref, accessions in self.uniref50.items():
            if uniprotid in accessions:
                return uniref
        return False

    def get_cluster(self, protein_list):
        """
        Search clustering information in UniRef from UniProt for a given UniProt accession

        Args:
            None

        """
        print("Clustering UniProt accessions unintegrated with no signature using Uniref50")
        for uniprotid in protein_list:
            uniref_cluster = self.search_uniprotid_in_uniref(uniprotid)

            if uniref_cluster:
                self.clusters.setdefault(uniref_cluster, []).append(uniprotid)
            else:
                url = f"https://www.uniprot.org/uniref/?query={uniprotid}&fil=identity:0.5&columns=id,members&format=tab"
                response = requests.get(url)
                data = response.text

                if response.status_code != 200:
                    print(f"FAILURE::{url}")

                uniref_all = data.split("\n")[1:]

                for uniref_info in uniref_all:
                    if uniref_info:
                        name, accessions = uniref_info.split("\t")
                        accessions = accessions.split("; ")

                        if name not in self.uniref50:
                            self.uniref50[name] = accessions
                        self.clusters.setdefault(name, []).append(uniprotid)


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

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    user = config["database"]["user"]
    password = config["database"]["password"]
    schema = config["database"]["schema"]
    outputdir = config["files"]["outputdir"]

    # initialising
    protein_pip = protein_pipeline()
    protein_pip.getConnection(user, password, schema)

    # create output directory if it doesn't exist
    os.makedirs(outputdir, exist_ok=True)

    # initialise tax_id value
    if args.organism:
        print(f"Searching taxid for {args.organism}")
        protein_pip.search_taxid(args.organism)
    elif args.taxid:
        protein_pip.tax_id = args.taxid
    else:
        print("Error no organism or taxid provided")
        sys.exit(1)

    # search the proteome
    print(f"Searching list of proteins for {protein_pip.tax_id}")
    protein_list = protein_pip.get_proteins()

    # search for integrated proteins
    list_integrated = protein_pip.get_integrated(protein_list)
    print(f"UniProt accessions integrated: {len(list_integrated)}")

    # list of unintegrated proteins
    unintegrated_subset = set(protein_list).difference(list_integrated)
    print(f"UniProt accessions unintegrated: {len(unintegrated_subset)}")

    # search for proteins in unintegrated InterPro signatures
    list_in_signature = protein_pip.get_accession_in_signature(outputdir, unintegrated_subset)

    # list_in_signature = set(
    #     protein_pip.get_accession_in_signature(args.folder, unintegrated_subset)
    # )
    print(f"UniProt accession unintegrated matching signature: {list_in_signature}")

    # list of unintegrated proteins not found in InterPro signatures
    # list_not_in_signature = unintegrated_subset.difference(list_in_signature)
    list_not_in_signature = len(unintegrated_subset) - list_in_signature
    print(f"UniProt accession unintegrated with no signature: {list_not_in_signature}")

    # close database connection
    protein_pip.connection.close()

    # # clustering unintegrated proteins
    # protein_pip.get_cluster(list_not_in_signature)
    # print(f"{len(protein_pip.clusters)} clusters found")

    # # write clustering results in file
    # cluster_file = os.path.join(args.folder, f"clusters_proteome_taxid_{protein_pip.tax_id}.csv")
    # with open(cluster_file, "w") as f:
    #     f.write("cluster_id,accessions\n")
    #     for cluster, accessions in protein_pip.clusters.items():
    #         f.write(f"{cluster},{'; '.join(accessions)}\n")

    # uniref50_cluster_file = os.path.join(
    #     args.folder, f"all_clusters_taxid_{protein_pip.tax_id}.csv"
    # )
    # with open(uniref50_cluster_file, "w") as f:
    #     f.write("cluster_id,count proteome matches,accessions\n")
    #     for cluster, accessions in protein_pip.uniref50.items():
    #         f.write(f"{cluster},{len(protein_pip.clusters[cluster])},{'; '.join(accessions)}\n")

