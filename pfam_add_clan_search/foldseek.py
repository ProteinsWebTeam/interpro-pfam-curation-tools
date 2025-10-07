"""
This script performs the following steps:

1.  If foldseek output file does not exist, run foldseek to find similar protein structures. Otherwise, skip.
2.  Process the foldseek results i.e. find Pfam clans and write them to a TSV file.

Note: if you'd like to run foldseek and not just process the results, 
delete the foldseek output file (foldseekpfam_all.out) before running this script.

Usage: python3 foldseek.py CONFIG_FILE

Thresholds applied to the foldseek results:
- Protein overlap: 60%
- E-value: 1e-3
"""

import subprocess
import os
import argparse
from configparser import ConfigParser
from utils import write_results_to_tsv, sql_connection

def run_foldseek(output_dir, pfam_acc, input_file, afdb_dir, tmp_dir):
    """Run foldseek to find similar protein structures."""
    print("Running foldseek...")
    output_file = os.path.join(output_dir, f"foldseek_{pfam_acc}.out")
    cmd = ["foldseek", "easy-search", input_file, afdb_dir, output_file, tmp_dir]
    with open(os.path.join(output_dir,f"log_{pfam_acc}.txt"),"wb") as f:
        popen = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, shell=False)
        return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

    return output_file

def run_foldseek_all(output_file, log_file, input_dir, afdb_dir, tmp_dir):
    """Run foldseek to find similar protein structures."""
    print("Running foldseek...")
    cmd = ["foldseek", "easy-search", input_dir, afdb_dir, output_file, tmp_dir]
    with open(log_file,"wb") as f:
        popen = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, shell=False)
        return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

    return output_file

def create_foldseek_pfamdb(input_dir, afpfamdb_dir, log_file):
    """Create Foldseek database with AF structures chopped to Pfam boundaries."""
    print("Creating Foldseek database...")
    cmd = ["foldseek", "createdb", input_dir, afpfamdb_dir]
    with open(log_file,"wb") as f:
        popen = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, shell=False)
        return_code = popen.wait()
    if return_code:
        print(return_code)
        raise subprocess.CalledProcessError(return_code, cmd)


def get_pfam_with_clan(dbinfo):
    """Extract Pfam entries of type 'domain' or 'family' with clans from the database."""
    conn = sql_connection(dbinfo)
    cursor = conn.cursor(dictionary=True)
    cursor.execute("SELECT p.pfamA_acc, clan_acc FROM pfamA p join clan_membership cl on cl.pfamA_acc=p.pfamA_acc")
    
    pfams = {}
    for row in cursor.fetchall():
        pfam_acc = row['pfamA_acc']
        clan_acc = row['clan_acc']
        pfams[pfam_acc] = clan_acc

    cursor.close()
    conn.close()
    return pfams

def process_foldseek_results_pfam(foldseek_file, output_file, dbinfo):
    """Process the foldseek results and write them to a TSV file."""
    pfam_with_clan = get_pfam_with_clan(dbinfo)

    with open(foldseek_file, "r") as f, open(output_file, "w") as out:
        current_pfam = None
        for line in f.readlines():
            line = line.strip().split()
            pfam_initial_data = line[0]
            pfam_acc = pfam_initial_data.split("_")[0]
            
            pfam_match_data = line[1]
            pfammatch_acc = pfam_match_data.split("_")[0]

            fident = line[2]

            start, end = pfam_initial_data.split("_")[-1].split("-")
            start = int(start.replace("res", ""))
            end = int(end)
            length = end - start

            startm, endm = pfam_match_data.split("_")[-1].split("-")
            startm = int(startm.replace("res", ""))
            endm = int(endm)
            lengthm = endm - startm
            
            alnlen = int(line[3])
            evalue = float(line[10])
            bitscore = int(line[11])

            # if pfam_acc != pfammatch_acc and pfam_acc not in pfam_with_clan and pfammatch_acc in pfam_with_clan and alnlen > 0.6*length and alnlen > 0.6*lengthm and evalue <= 1e-3:
            if pfam_acc != pfammatch_acc and pfam_acc not in pfam_with_clan and alnlen > 0.6*length and alnlen > 0.6*lengthm and evalue <= 1e-3:
                print(pfam_acc)
                try:
                    clan_acc = pfam_with_clan[pfammatch_acc]
                except KeyError:
                    clan_acc = "NoClan"
                out.write(f"{pfam_initial_data}\t{pfam_match_data}\t{clan_acc}\t{fident}\t{alnlen}\t{evalue}\t{bitscore}\n")


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

    output_dir = config["files"]["file_dir"]
    tmp_dir = config["files"]["tmp_dir"]

    afdb_dir = config["files"]["afdbpfam_dir"]
    
    # initial_pfam_data = {'pfamA_acc': 'PF00032', 'pfamseq_acc': 'P42792', 'seq_version': 3, 'seq_start': 265, 'seq_end': 366}
    output_file_foldseek = os.path.join(output_dir, "foldseekpfam_all.out")
    log_file = os.path.join(output_dir, "foldseekpfam_all.log")
    output_map_file = os.path.join(output_dir, "foldseekpfam_map_all.tsv")

    if os.path.isfile(output_file_foldseek) and os.path.getsize(output_file_foldseek) > 0:
        print(f"already ran foldseek, skipping")
        if os.path.isfile(output_map_file):
            os.remove(output_map_file)

        process_foldseek_results_pfam(output_file_foldseek, output_map_file, dbinfo)

    else:
        output_file_foldseek = run_foldseek_all(output_file_foldseek, log_file, output_dir, afdb_dir, tmp_dir)
        process_foldseek_results_pfam(output_file_foldseek, output_map_file, dbinfo)