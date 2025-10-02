Create folders to build potential new Pfam families from smORF encoded proteins mentioned in scientific literature
=========================================================

Beatriz Lazaro Pinto (2024)

### Main pipeline 
``` bash
python3 smORFs_main.py config.ini
```

This script performs the following steps:
1. gets Uniprot IDs of small protein sequences that are linked to a EuropePMC query specified in the config.ini and are not in the Pfam version used in the Uniprot website
2. creates folders for each with a SEED file in Stockolm format and a dummy file with the corresponding PMID
3. For each folder:
    - checks if the SEED sequence is in pfamseq
    - creates an HMM, and when this job is finished:
        - get SwissProt matches
        - add pdb references
        - get taxonomy distribution
        - amends the DESC file
            - adds author, PMID and common text in CC lines
            - add ID, DE lines and text in the CC lines if possible with information obtained from Uniprot
        - check overlaps:
            - if the percentage of overlaps is > 90%, renames the folder to OVERLAPS_x