# interpro-curation-tools
Scripts to help InterPro curation

##Requirements
Python 3.4 and above

## Search proteome
Pipeline to search proteins for a given organism or taxid in InterPro signatures.
If unintegrated: search in unintegrated InterPro signatures or cluster proteins using UniRef50 (UniProt)
Usage: `python search_proteome -u USER -p PASSWORD -s DB_SCHEMA [-o ORGANISM | -t TAXID] -f FOLDER`
