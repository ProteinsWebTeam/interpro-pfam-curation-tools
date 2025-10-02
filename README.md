# interpro-pfam-curation-tools
Scripts to help InterPro and Pfam curation, this project has multiple modules

## Requirements
Python 3.6 and above

## Pfam DUF
The scripts are available under the **pfam_deduf** subdirectory

### List of Pfam entries that could potentially be deDUFFED
Execute queries in the Pronto postgres database (`find_deduf.sql`) to find Pfam entries that could potentially be deDUFFED.
- Entries that have PDB strucutres
- Entries that have SwissProt matches

### Search Pfam DUF in InterPro entries and count swissprot matches
Usage: `python pfam_duf/deduf_main.py config_pfam.ini [-o OPTION]: 1)search swissprot names for unintegrated pfam DUF 2)search for DUF in literature 3)search GO term/keywords in Swissprot matches 4)search predicted structures`

## Identifying Pfam to add to clans
See README file in the **pfam_add_clan_search** subdirectory for instructions.

## Identify InterPro entries that don't have Pfam coverage
See README file in the **interpro_no_pfam** subdirectory for instructions.

## COGs analysis
The scripts are available under the **cogs_analysis** subdirectory. 
More details are provided in the corresponding README file.

## Create folders to build potential new Pfam families from smORF encoded proteins
See README file in the **smORF_new_Pfam_families** subdirectory for instructions.