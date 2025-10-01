"""
This script performs the following steps:

1.  Extract Pfam representatives from the output directory of chopped alphafold structures.
2.  Compare the AF boundaries of the Pfam representatives with the Pfam boundaries in the database.
3.  If the boundaries have changed, update the AF boundaries in the output directory.

> Note: this script should be run before running the main pipeline.

Usage: python3 update_af.py CONFIG_FILE
"""

from utils import sql_connection, extract_pfam_representatives_from_file
from split_af import split_alphafold
import argparse
from configparser import ConfigParser
import os


def update_af_boundaries(pfam_acc, protein_info, dbinfo):
    conn = sql_connection(dbinfo)
    cursor = conn.cursor(dictionary=True)

    pfam_seq_acc = protein_info['pfamseq_acc']

    cursor.execute("select pfamA_acc, pfamseq_acc, seq_version, seq_start, seq_end from pfamA_reg_seed where pfamA_acc=%s and pfamseq_acc=%s", (pfam_acc, pfam_seq_acc))
    for row in cursor.fetchall():
        if row['seq_start'] != protein_info['seq_start'] or row['seq_end'] != protein_info['seq_end']:
            print(f"Changed boundaries for {pfam_acc}")
            af_to_update.append(row)
            file_name = os.path.join(input_dir, f"{pfam_acc}_{pfam_seq_acc}_{protein_info['seq_version']}_res{protein_info['seq_start']}-{protein_info['seq_end']}.cif")
            print(f"deleting {file_name}")
            os.remove(file_name)
            out_file = split_alphafold(protein_info['pfamseq_acc'], protein_info['seq_version'], pfam_acc, chopped_struct_dir, int(protein_info['seq_start']), int(protein_info['seq_end']), output_format)

    cursor.close()
    conn.close()

    return None


def update_af_with_pfam_changed_boundaries(dbinfo, input_dir, file):
    results, outfile = extract_pfam_representatives_from_file(input_dir_chopped_af)

    for pfam_acc, protein_info in results.items():
        print(protein_info["pfamA_acc"])
        update_af_boundaries(pfam_acc, protein_info, dbinfo)

            
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("config", metavar="CONFIG_FILE", help="configuration file")

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

    input_dir = config["files"]["file_dir"]
    input_dir_chopped_af = os.path.join(input_dir, "chopped_cif")
    afdb_dir = config["files"]["afdbpfam_dir"]

    update_af_with_pfam_changed_boundaries(dbinfo, input_dir_chopped_af)
