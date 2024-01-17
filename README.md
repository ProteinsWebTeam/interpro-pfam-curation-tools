# interpro-pfam-curation-tools
Scripts to help InterPro and Pfam curation, this project has multiple modules

## Requirements
Python 3.6 and above

## Proteome analysis
The scripts are available under the **proteome_count** subdirectory

### Search proteome
Pipeline to search proteins for a given organism or taxid in InterPro signatures.
If unintegrated: search in unintegrated InterPro signatures or cluster proteins using UniRef50 (UniProt)
Usage: `python proteome_count/search_proteome.py config.ini [-o ORGANISM | -t TAXID]`

### Count integrated entries for an organism
Script to search the number of integrated InterPro entries for a given reference proteome, organism or taxid in a given period of time
Usage: `python proteome_count/count_all_proteome.py config.ini [-r REF_PROTEOME | -t TAXID | -o ORGANISM]`

### Count all integrated proteomes
Script to search the number of integrated InterPro entries corresponding to the key proteomes from the InterPro/Pfam BBR grant in a given period of time
Usage: `python proteome_count/count_all_proteome.py config.ini [-r taxid | proteome]`


## Pfam DUF
The scripts are available under the **pfam_deduf** subdirectory

### Search Pfam DUF in InterPro entries and count swissprot matches
Usage: `python pfam_duf/deduf_main.py config_pfam.ini [-o OPTION]: 1)search swissprot names for unintegrated pfam DUF 2)search for DUF in literature 3)search GO term/keywords in Swissprot matches 4)search predicted structures`

## COGs analysis
The scripts are available under the **cogs_analysis** subdirectory. 
More details are provided in the corresponding README file.

## Search CATH names
This script searches for assigned names to CATH-Gene3D signatures not yet integrated in InterPro.
Usage: `python search_new_CATH_names.py -u USER -p PASSWORD -s DB_SCHEMA -f INPUTFILE -o OUTPUTFILE [-c COMMENTFILE]`