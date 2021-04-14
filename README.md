# interpro-curation-tools
Scripts to help InterPro and Pfam curation, this project has multiple modules

## Requirements
Python 3.6 and above

## Proteome analysis
The scripts are available under the **proteome_count** subdirectory

### Search proteome
Pipeline to search proteins for a given organism or taxid in InterPro signatures.
If unintegrated: search in unintegrated InterPro signatures or cluster proteins using UniRef50 (UniProt)
Usage: `python proteome_count/search_proteome.py -u USER -p PASSWORD -s DB_SCHEMA [-o ORGANISM | -t TAXID] -f FOLDER`

### Count integrated entries for an organism
Script to search the number of integrated InterPro entries for a given organism or taxid in a given period of time
Usage: `python proteome_count/count_all_proteome.py -u USER -p PASSWORD -s DB_SCHEMA [-o ORGANISM | -t TAXID] -b begin_date [-e end_date]`

### Count all integrated proteomes
Script to search the number of integrated InterPro entries corresponding to the key proteomes from the InterPro/Pfam BBR grant in a given period of time
Usage: `python proteome_count/count_all_proteome.py -u USER -p PASSWORD -s DB_SCHEMA -b begin_date [-e end_date]`


## Pfam DUF
The scripts are available under the **pfam_deduf** subdirectory

### Search Pfam DUF in InterPro entries and count swissprot matches
Usage: `python pfam_duf/deduf_main.py -u interpro -p ***REMOVED*** -s IPPRO -d results/duf -o 1`

### Search Pfam DUF in litterature
Usage: `python pfam_duf/deduf_main.py -u interpro -p ***REMOVED*** -s IPPRO -d results/duf -o 2`


## COGs analysis
The scripts are available under the **cogs_analysis** subdirectory. 
More details are provided in the corresponding README file.

## Search CATH names
This script searches for assigned names to CATH-Gene3D signatures not yet integrated in InterPro.
Usage: `python search_new_CATH_names.py -u USER -p PASSWORD -s DB_SCHEMA -f INPUTFILE -o OUTPUTFILE [-c COMMENTFILE]`