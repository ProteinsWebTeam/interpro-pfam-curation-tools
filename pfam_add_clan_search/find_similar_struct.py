"""
This script performs the following steps:

1. Extract Pfam entries of type 'domain' or 'family' with clans from the Pfam SQL database.
2. Select a representative protein from the SEED for a given Pfam accession (highest AF plDDT score).
3. Split AlphaFold structure to Pfam boundaries.
4. Build foldseek database.
5. Run foldseek to find similar protein structures (easy-search default parameters).
6. Process the foldseek results and write them to a TSV file.

Usage: python3 find_similar_struct.py CONFIG_FILE [-p PFAM] [-f PFAM_LIST]
"""

import requests
import subprocess
import time, re, os
import argparse
import shutil
from configparser import ConfigParser
from split_af import split_alphafold
from foldseek import run_foldseek_all, create_foldseek_pfamdb, process_foldseek_results_pfam
from utils import sql_connection, write_results_to_tsv, extract_pfam_representatives_from_file
from similarity_graph import analyse_foldseek_results

def extract_pfams(dbinfo):
    """Extract Pfam entries of type 'domain' or 'family' from the database."""
    conn = sql_connection(dbinfo)
    cursor = conn.cursor(dictionary=True)
    cursor.execute("SELECT p.pfamA_acc, p.type, clan_acc FROM pfamA p left join clan_membership cl on cl.pfamA_acc=p.pfamA_acc and p.type in ('Family', 'Domain')")
    pfams = cursor.fetchall()
    cursor.close()
    conn.close()
    return pfams


def select_representative_protein(pfam_acc, dbinfo):
    """Select a representative protein from the SEED for a given pfam_acc."""
    conn = sql_connection(dbinfo)
    cursor = conn.cursor(dictionary=True)
    cursor.execute("select pfamA_acc, pfamseq_acc, seq_version, seq_start, seq_end from pfamA_reg_seed where pfamA_acc=%s", (pfam_acc,))
    proteins = cursor.fetchall()
    cursor.close()
    conn.close()
    representative = ""
    max_plddt = 0
    for row in proteins:
        protein = row['pfamseq_acc']
        plddt = get_af_plddt(protein)
        if plddt > max_plddt:
            representative = row
            max_plddt = plddt

    return representative


def get_af_plddt(protein):
    """Get the AF PLDDT score for a given protein."""
    url = f"https://alphafold.ebi.ac.uk/api/uniprot/summary/{protein}.json"
    r = requests.get(url)
    if r.status_code == 200:
        data = r.json()
        for structure in data.get('structures', []):
            score = structure['summary'].get('confidence_avg_local_score')
        return score
    return 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="CONFIG_FILE", help="configuration file")
    parser.add_argument("-p", "--pfam", help="Pfam accession")
    parser.add_argument("-f", "--file", help="List of Pfam accessions to process")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"Cannot open '{args.config}': " f"no such file or directory")

    config = ConfigParser()
    config.read(args.config)

    #mysql connection info
    host = config["db"]["host"]
    user = config["db"]["user"]
    password = config["db"]["password"]
    database = config["db"]["database"]
    port = config["db"]["port"]
    dbinfo = [host, user, password, database, port]

    output_format = config["files"]["output_format"]
    output_dir = config["files"]["file_dir"]
    chopped_struct_dir = os.path.join(output_dir, "chopped_cif")
    tmp_dir = config["files"]["tmp_dir"]

    afdb_dir = config["files"]["afdbpfam_dir"]

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(chopped_struct_dir, exist_ok=True)
    
    if args.pfam:
        pfams = [{'pfamA_acc':args.pfam}]
        pfam_id = args.pfam
    elif args.file:
        pfams = []
        with open(args.file, 'r') as f:
            for line in f.readlines():
                pfam_acc, pfam_type, pfam_clan = line.strip().split('\t')
                pfams.append({'pfamA_acc':pfam_acc, 'type':pfam_type, 'clan_acc':pfam_clan})
        pfam_id = None
    else:
        pfams = extract_pfams(dbinfo)
        pfam_id = None

    chopped_pfams = set()
    chopped_pfams, out_file = extract_pfam_representatives_from_file(chopped_struct_dir, pfam_id)

    for pfam in pfams:
        pfam_acc = pfam['pfamA_acc']
        print(f"Processing {pfam_acc}")
        # search representative protein accession
        if pfam_acc in chopped_pfams:
            protein_info = chopped_pfams[pfam_acc]
        else:
            protein_info = select_representative_protein(pfam_acc, dbinfo)

            # chop AF structure to Pfam boundaries
            if protein_info:
                out_file = split_alphafold(protein_info['pfamseq_acc'], protein_info['seq_version'], pfam_acc, chopped_struct_dir, int(protein_info['seq_start']), int(protein_info['seq_end']), output_format)
            else:
                out_file = None

    parent_dir = os.path.dirname(afdb_dir)
    try:
        print(f"Deleting old {parent_dir} directory")
        shutil.rmtree(parent_dir)
    except FileNotFoundError:
        pass
    os.makedirs(parent_dir, exist_ok=True)
    create_foldseek_pfamdb(chopped_struct_dir, afdb_dir, os.path.join(output_dir, "create_foldseek_db.log"))

    if args.pfam:
        output_file_foldseek = os.path.join(output_dir, f"foldseek_{args.pfam}.out")
        log_file = os.path.join(output_dir, f"foldseek_{args.pfam}.log")
        output_map_file = os.path.join(output_dir, f"foldseek_{args.pfam}_map.tsv")
        if not os.path.isfile(output_file_foldseek) or os.path.getsize(output_file_foldseek) == 0:
            foldseek_result_file = run_foldseek_all(output_file_foldseek, log_file, out_file, afdb_dir, tmp_dir)
            process_foldseek_results_pfam(output_file_foldseek, output_map_file, dbinfo)
            analyse_foldseek_results(output_map_file)

    else:
        output_file_foldseek = os.path.join(output_dir, f"foldseek_all.out")
        log_file = os.path.join(output_dir, f"foldseek_all.log")
        output_map_file = os.path.join(output_dir, f"foldseek_all_map.tsv")

        if not os.path.isfile(output_file_foldseek) or os.path.getsize(output_file_foldseek) == 0:
            foldseek_result_file = run_foldseek_all(output_file_foldseek, log_file, chopped_struct_dir, afdb_dir, tmp_dir)

            process_foldseek_results_pfam(output_file_foldseek, output_map_file, dbinfo)
            analyse_foldseek_results(output_map_file)

if __name__ == "__main__":
    main()
