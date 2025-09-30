"""
This script performs the following steps:

1.  Download AlphaFold structure for the given UniProt accession from the AlphaFold website.
2.  Split the AlphaFold structure to Pfam boundaries.

Usage: python3 split_af.py pfam_acc uniprot_acc uniprot_version output_dir [--format {cif,pdb}] start end
"""

import sys
import os
import requests
from Bio.PDB import MMCIFParser, MMCIFIO, PDBIO, Select
import argparse
import configparser
import time


class ResidueRangeSelect(Select):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def accept_residue(self, residue):
        res_id = residue.get_id()[1]
        return self.start <= res_id <= self.end

def download_alphafold_model(uniprot_acc, out_path):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_acc}-F1-model_v4.cif"

    r = requests.get(url)

    if r.status_code == 200:
        with open(out_path, "wb") as f:
            f.write(r.content)
        print(f"Downloaded AlphaFold model for {uniprot_acc} to {out_path}")
        return True

    print(f"Failed to download AlphaFold model for {uniprot_acc} (HTTP {r.status_code})")
    return False
    

def split_structure(input_cif, boundaries, output_prefix, out_format="cif"):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('model', input_cif)
    if out_format == "pdb":
        io = PDBIO()
        ext = "pdb"
    else:
        io = MMCIFIO()
        ext = "cif"
    for i, (start, end) in enumerate(boundaries, 1):
        select = ResidueRangeSelect(start, end)
        out_file = f"{output_prefix}_res{start}-{end}.{ext}"
        if os.path.exists(out_file):
            print(f"Already exists {out_file}")
            os.remove(input_cif)
            return None
        io.set_structure(structure)
        io.save(out_file, select=select)
        # print(f"Saved: {out_file}")
        os.remove(input_cif)
        return out_file
    return None


def split_alphafold(uniprot_acc, seq_version, pfam_acc, output_dir, start, end, output_format="cif"):
    output_prefix = os.path.join(output_dir, f"{pfam_acc}_{uniprot_acc}_{seq_version}")
    cif_file = os.path.join(output_dir,f"AF-{uniprot_acc}-F1-model_v4.cif")
    if not os.path.exists(cif_file):
        if not download_alphafold_model(uniprot_acc, cif_file):
            return ""

    return split_structure(cif_file, [(start, end)], output_prefix, output_format)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Split AlphaFold structure by residue boundaries.")
    parser.add_argument("pfam_acc", help="UniProt accession")
    parser.add_argument("uniprot_acc", help="UniProt accession")
    parser.add_argument("uniprot_version", help="UniProt version")
    parser.add_argument("output_dir", help="Directory for output files")
    parser.add_argument("--format", choices=["cif", "pdb"], default="cif", help="Output format")
    parser.add_argument("start", type=int, help="Residue start: start")
    parser.add_argument("end", type=int, help="Residue end: end")

    args = parser.parse_args()

    out_file = split_alphafold(args.uniprot_acc, args.uniprot_version, args.pfam_acc, args.output_dir, args.start, args.end, args.format)