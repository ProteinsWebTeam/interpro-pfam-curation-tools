#!/usr/bin/env python3

"""
@author T. Paysan-Lafosse

@brief This script searches unintegrated proteins for a given organism or taxid 
        in InterPro signatures, if they are not found in signatures, they are clustered based on UniRef clusters

@arguments [-u USER]: database user
           [-p PASSWORD]: database password for the user
           [-s SCHEMA]: database schema to use
           [-o ORGANISM or -t TAXID]: organism (scientific name) or taxid to look for
           [-f FOLDER]: output folder
"""

import argparse
import configparser
import os
import sys
import traceback
from pathlib import Path

import cx_Oracle
import requests


class protein_pipeline:
    def __init__(self, user, password, schema):
        self.user, self.password = user, password
        self.schema = schema
        self.connection = self.getConnection()
        self.cursor = self.connection.cursor()
        self.tax_id = None
        self.uniref50 = dict()
        self.clusters = dict()

    def getConnection(self):
        connectString = "".join([self.user, "/", self.password, "@", self.schema])
        try:
            return cx_Oracle.connect(connectString)
        except:
            stackTrace = traceback.format_exc()
            print(stackTrace)
            if "invalid username" in stackTrace:
                print(
                    "Could not connect to {0} as user {1}".format(
                        self.schema, self.user
                    )
                )
                print(
                    "NB if your oracle username contains the '$' character either escape it or surround it with quotes"
                )
                print('eg "ops$craigm" or ops\$craigm')
                print(
                    "Otherwise the shell will remove the '$' and all subsequent characters!"
                )
            sys.exit(1)

    def chunks(self, l, n):
        """Yield chunks of size n from iterable.

        Args:
            l (Iterable): The iterable to chunk
            n (int): Maximum number of items in chunk

        Yields:
            Iterable: Tuples of length n or less (final bucket will have only the remaining items)

        """
        for i in range(0, len(l), n):
            # Create an index range for l of n items:
            yield l[i : i + n]

    def search_taxid(self, organism):
        """
        Search taxid for specified organism

        Args:
            configfile: Configuration file containing database connection credentials
            organism: Organism to search the taxid for

        Yields:
            taxid: taxonomic identifier for the organism of interest

        """

        request = "Select tax_id from INTERPRO.ETAXI where scientific_name=:1"
        self.cursor.execute(request, (organism,))
        self.tax_id = self.cursor.fetchone()[0]

        print(f"Found taxid {self.tax_id} for {organism}")

    def get_proteins(self):
        """
        Search protein accessions for the taxid

        Args: None

        Yields: 
            protein_list: list of proteins

        """

        print(f"Searching for protein accessions for taxid {self.tax_id}")

        request = "select unique p.protein_ac \
            from INTERPRO.PROTEIN p \
            join INTERPRO.ETAXI et on p.tax_id = et.tax_id \
            where p.tax_id=:1"

        self.cursor.execute(request, (self.tax_id,))
        protein_list = [row[0] for row in self.cursor]
        print(f"Found {len(protein_list)} accessions")

        return protein_list

    def get_integrated(self, protein_list):
        """
        Search integrated proteins

        Args:
            protein_list: list of proteins

        Yields:
            list_integrated: list of proteins integrated in InterPro entries

        """

        print("Searching for integrated proteins")
        uniprot_chunks = list(self.chunks(protein_list, 1000))

        list_integrated = set()
        for chunk in uniprot_chunks:
            protein_list_quote = [f"'{row}'" for row in chunk]
            request = f"select p.protein_ac \
                    from INTERPRO.MV_ENTRY2PROTEIN e2p \
                    join Interpro.protein p on e2p.PROTEIN_AC=p.PROTEIN_AC \
                    where e2p.protein_ac in ({','.join(protein_list_quote)})"
            self.cursor.execute(request)

            list_integrated.update(set([row[0] for row in self.cursor]))
        return list_integrated

    def get_count_signature_taxid(self, list_signatures):
        """
            Search for protein counts for a list of InterPro signatures

        Args:
            list_signatures: list of InterPro signatures
        
        Yields:
            dictionnary with signature as key and protein_count as value

        """

        signature_chunks = list(self.chunks(list(list_signatures), 1000))
        for chunk in signature_chunks:
            signature_list_quote = [f"'{row}'" for row in chunk]
        request = f"select m2p.method_ac,count(p.protein_ac) \
                from INTERPRO.PROTEIN p \
                join INTERPRO.MV_METHOD2PROTEIN m2p on p.PROTEIN_AC = m2p.PROTEIN_AC \
                join INTERPRO.ETAXI et on p.tax_id = et.tax_id \
                where et.tax_id=:1 and m2p.METHOD_AC in ({','.join(signature_list_quote)}) \
                group by m2p.method_ac"
        self.cursor.execute(request, (self.tax_id,))
        return {row[0]: row[1] for row in self.cursor}

    def get_accession_in_signature(self, folder, protein_list):
        """
        Search for proteins found in InterPro signatures but not integrated
        Write the results in a csv file with each row corresponding to a protein/signature pair (protein,dbcode,organism,signature,total_prot_count,count_taxa,comment)

        Args:
            protein_list: list of proteins
        
        Yields:
            list of proteins found in unintegrated signatures

        """

        print("Searching for unintegrated proteins in signature")
        uniprot_chunks = list(self.chunks(list(protein_list), 1000))
        list_signatures = set()
        list_proteins_with_signature = dict()

        for chunk in uniprot_chunks:
            protein_list_quote = [f"'{row}'" for row in chunk]
            request = f"select p.protein_ac, p.dbcode, et.scientific_name, m2p.METHOD_AC, mm.PROTEIN_COUNT, LISTAGG(mc.VALUE, '; ') WITHIN GROUP (ORDER by mc.value) COMMENTS \
                    from INTERPRO.PROTEIN p \
                    join INTERPRO.ETAXI et on p.tax_id = et.tax_id \
                    join INTERPRO.MV_METHOD2PROTEIN m2p on p.PROTEIN_AC = m2p.PROTEIN_AC \
                    join INTERPRO.MV_METHOD_MATCH mm on mm.method_ac = m2p.method_ac \
                    left join INTERPRO.METHOD_COMMENT mc on mc.method_ac = m2p.method_ac \
                    where p.protein_ac in ({','.join(protein_list_quote)}) \
                    group by p.protein_ac, p.dbcode, et.scientific_name, m2p.METHOD_AC, mm.PROTEIN_COUNT"
            self.cursor.execute(request)

            results = self.cursor.fetchall()
            for row in results:
                protein = row[0]
                signature = row[3]
                list_signatures.add(signature)
                try:
                    list_proteins_with_signature[protein][signature] = [
                        row[1],
                        row[2],
                        row[4],
                        row[5],
                    ]
                except KeyError:
                    list_proteins_with_signature[protein] = dict()
                    list_proteins_with_signature[protein][signature] = [
                        row[1],
                        row[2],
                        row[4],
                        row[5],
                    ]

        count_prot_signatures = self.get_count_signature_taxid(list_signatures)

        unintegrated_file = os.path.join(
            folder, f"unintegrated_prot_in_signatures_{self.tax_id}.csv"
        )
        with open(unintegrated_file, "w") as outf:
            outf.write(
                "protein,dbcode,organism,signature,total_prot_count,count_taxa,comment\n"
            )
            for protein, signatures in list_proteins_with_signature.items():
                for signature, values in signatures.items():
                    outf.write(
                        f"{protein},{values[0]},{values[1]},{signature},{values[2]},{count_prot_signatures[signature]},{values[3] if values[3]!='None' else ''}\n"
                    )
        return list_proteins_with_signature.keys()

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
            protein_list: list of UniProt accessions to cluster

        """
        print(
            "Clustering UniProt accessions unintegrated with no signature using Uniref50"
        )
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
    parser.add_argument(
        "-u", "--user", help="username for database connection", required=True
    )
    parser.add_argument(
        "-p", "--password", help="password for database connection", required=True
    )
    parser.add_argument(
        "-s", "--schema", help="database schema to connect to", required=True
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-o",
        "--organism",
        help="Scientific name of the organism to get the conservation score for",
    )
    group.add_argument(
        "-t", "--taxid", help="Taxid of the organism to get the conservation score for"
    )
    parser.add_argument(
        "-f", "--folder", help="folder directory to write output files", required=True
    )

    args = parser.parse_args()

    # initialise the class
    protein_pip = protein_pipeline(args.user, args.password, args.schema)

    # create output directory if it doesn't exist
    Path(args.folder).mkdir(parents=True, exist_ok=True)

    # initialise tax_id value
    if args.organism:
        print(f"Searching taxid for {args.organism}")
        protein_pip.search_taxid(args.organism)
    elif args.taxid:
        protein_pip.tax_id = args.taxid
    else:
        print("Error no organism or taxid provided")
        exit()

    # search the proteome
    protein_list = protein_pip.get_proteins()

    # search for integrated proteins
    list_integrated = protein_pip.get_integrated(protein_list)
    print(f"UniProt accessions integrated: {len(list_integrated)}")

    # list of unintegrated proteins
    unintegrated_subset = set(protein_list).difference(list_integrated)
    print(f"UniProt accessions unintegrated: {len(unintegrated_subset)}")

    # search for proteins in unintegrated InterPro signatures
    list_in_signature = set(
        protein_pip.get_accession_in_signature(args.folder, unintegrated_subset)
    )
    print(
        f"UniProt accession unintegrated matching signature: {len(list_in_signature)}"
    )

    # list of unintegrated proteins not found in InterPro signatures
    list_not_in_signature = unintegrated_subset.difference(list_in_signature)
    print(
        f"UniProt accession unintegrated with no signature: {len(list_not_in_signature)}"
    )

    # close database connection
    protein_pip.connection.close()

    # clustering unintegrated proteins
    protein_pip.get_cluster(list_not_in_signature)
    print(f"{len(protein_pip.clusters)} clusters found")

    # write clustering results in file
    cluster_file = os.path.join(
        args.folder, f"clusters_proteome_taxid_{protein_pip.tax_id}.csv"
    )
    with open(cluster_file, "w") as f:
        f.write("cluster_id,accessions\n")
        for cluster, accessions in protein_pip.clusters.items():
            f.write(f"{cluster},{'; '.join(accessions)}\n")

    uniref50_cluster_file = os.path.join(
        args.folder, f"all_clusters_taxid_{protein_pip.tax_id}.csv"
    )
    with open(uniref50_cluster_file, "w") as f:
        f.write("cluster_id,size,accessions\n")
        for cluster, accessions in protein_pip.uniref50.items():
            f.write(f"{cluster},{len(accessions)},{'; '.join(accessions)}\n")
