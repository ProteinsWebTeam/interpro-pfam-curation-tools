Finding Pfam clans using structural homology
============================================

This directory contains scripts for finding Pfam clans using structural homology.

Developed by Typhaine Paysan-Lafosse (September 2025)

## Requirements

- Python 3.11 or later
- Pandas
- Numpy
- NetworkX
- mysql.connector

Use a screen session and virtual environment to run the scripts in HPC cluster.

## Usage

### Main pipeline 
``` bash
python3 find_similar_struct.py CONFIG_FILE [-p PFAM] [-f PFAM_LIST]
```

> Note: Run on HPC cluster with the following machine configuration (for foldseek50 option): srun --gpus=1 -t 2:00:00 --mem=100G --pty bash

This script performs the following steps:

1. Extract Pfam entries of type 'domain' or 'family' with clans from the Pfam SQL database.
2. Select a representative protein from the SEED for a given Pfam accession (highest AF plDDT score).
3. Split AlphaFold structure to Pfam boundaries.
4. Build foldseek database.
5. Run foldseek to find similar protein structures (easy-search default parameters).
6. Process the foldseek results and write them to a TSV file (Thresholds applied to the foldseek results: Protein overlap: 60%, E-value: 1e-3).


### Identifying significant Pfam-Pfam homologues

``` bash
python3 similarity_graph.py outputdir/foldseekpfam_map_all.tsv
```

This script performs the following steps:

1.  Read the foldseek results and build a Pfam-Pfam similarity network.
2.  Suggest additions of Pfams to existing clans (based on neighbour agreement) *pfam_edges_all.csv*
3.  Discover *de novo* groups among Pfams with no clan (potential new clans): *proposed_groups_all.csv*
4.  Write the results to a csv file.

> Note: the *pfam_edges_all.csv* can be imported into https://cosmograph.app/ for visualization.

### Identifying Pfam without a chopped AlphaFold structure

```bash
python3 find_multi_af.py CONFIG_FILE
```

This script performs the following steps:

1.  Extract Pfam representatives from the output directory of chopped alphafold structures.
2.  Extract Pfam accessions from the Pfam SQL database.
3.  Compare the above lists and generate a list of Pfam accessions without a shopped AlphaFold structure.

### Updating chopped alphafold structures

``` bash
python3 update_af.py CONFIG_FILE
```

This script performs the following steps:

1.  Extract Pfam representatives from the output directory of chopped alphafold structures.
2.  Compare the AF boundaries of the Pfam representatives with the Pfam boundaries in the database.
3.  If the boundaries have changed, update the AF boundaries in the output directory.

> Note: this script should be run before running the main pipeline.


### Running foldseek only
``` bash
python3 foldseek.py CONFIG_FILE
```

This script performs the following steps:

1.  If foldseek output file does not exist, run foldseek to find similar protein structures. Otherwise, skip.
2.  Process the foldseek results i.e. find Pfam clans and write them to a TSV file.

> Note: if you'd like to run foldseek and not just process the results, delete the foldseek output file (foldseekpfam_all.out) before running this script.


### Chopping an AF structure for a given Pfam/list of Pfams
``` bash
python3 split_af.py pfam_acc uniprot_acc uniprot_version output_dir [--format {cif,pdb}] start end
```

This script performs the following steps:

1.  Download AlphaFold structure for the given UniProt accession from the AlphaFold website.
2.  Split the AlphaFold structure to Pfam boundaries.