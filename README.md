# interpro-curation-tools
Scripts to help InterPro curation

## Requirements
Python 3.6 and above

## Search proteome
Pipeline to search proteins for a given organism or taxid in InterPro signatures.
If unintegrated: search in unintegrated InterPro signatures or cluster proteins using UniRef50 (UniProt)
Usage: `python search_proteome -u USER -p PASSWORD -s DB_SCHEMA [-o ORGANISM | -t TAXID] -f FOLDER`

## Count integrated entries for an organism
Script to search the number of integrated InterPro entries for a given organism or taxid in a given period of time
Usage: `python count_all_proteome -u USER -p PASSWORD -s DB_SCHEMA [-o ORGANISM | -t TAXID] -b begin_date [-e end_date]`

## Count all integrated proteomes
Script to search the number of integrated InterPro entries corresponding to the key proteomes from the InterPro/Pfam BBR grant in a given period of time
Usage: `python count_all_proteome -u USER -p PASSWORD -s DB_SCHEMA -b begin_date [-e end_date]`
