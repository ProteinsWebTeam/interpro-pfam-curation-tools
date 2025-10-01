# Identifying InterPro entries and proteins that Pfam is missing

This pipeline is used to identify proteins matched by an InterPro entry, but not by a Pfam entry.

## Step 1: Finding the InterPro entries from the InterPro database

Connect to the InterPro Oracle database and run the `entry_protein_no_pfam.sql` script. It will geneate a file called `entry_no_pfam.csv` in the current directory.

> **Note:** This script will take a while to run, it is best being run on a screen session and a virtual machine (not login node).

## Step 2: Gathering the number of proteins for an InterPro entries and the overlap count to Pfam signatures

```bash
python3 search_entry_no_pfam.py config.ini
```

This script will create a file called `entry_no_pfam_counts.csv` in the `output` directory.
This file contains 4 columns:

1. InterPro entry accession
2. Number of proteins matched by the InterPro entry
3. Number of proteins matched by the InterPro entry and overlapping with a Pfam signature
4. Number of proteins matched by the InterPro entry and not overlapping with a Pfam signature

This script will also print a summary of the number of Initial InterPro entries, entries excluded (because a Pfam signature is included in them) and total number of entries to process.

> **Note:** This script will take a while to run, it is best being run on a screen session and a virtual machine (not login node).

## Step 3: Identifying proteins matched by an InterPro entry, but not by a Pfam entry

```bash
python3 parse_results.py config.ini
```

This script will create two files for each InterPro entry where the overlap count to Pfam signatures is equal to 0:

1. `IPRXXXXX.tsv`: contains the proteins matched by the InterPro entry, but not by a Pfam signature and the start and end positions of the InterPro entry
2. `IPRXXXXX.desc`: contains the description of the InterPro entry, including the InterPro name (DE), short name (ID), references and description (CC) when available, to help with the creation of the Pfam DESC file.

This script will also print a summary of the InterPro entries with overlap count to Pfam signatures by range (e.g 0 to 9, 10 to 19, 20 to 29 ... 60+).