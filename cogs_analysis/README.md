# JACCARD index calculation between COGs and integrated InterPro entries

## Requirements:

Python 3.6+, with package `cx_Oracle`, `redis`

## Configuration

This pipeline expect an INI config file.
Load Oracle environment: `source ~oracle/ora112setup.sh`

## Usage

The process is divided into 5 different steps, that can be run individually.

To run the full process (except the last step) use the following command line `python run_cog_analysis.py config.ini -1 -2 -3 -4`


1) The representative list of proteins from UNIPROT is loaded in a redis hash

2) For each protein in the representative list:
	
	For each integrated InterPro entry, we calculated the midpoint of the different fragments and count the number of protein matching
    
	All the data are saved in redis hashes:
	- protein_hash[protein_acc]={"IPR1":[midpoint1,midpoint2], "IPR2":[midpoint1]}
	- ipr_count[ipr]=count_protein

3) Get the list of COGs into a redis queue which will be used by 8 different jobs submitted in parallel for step 4

4) For each COG: 
	
	For each protein:
	- Get IPR and midpoints
	- Verify if intersection between COG and IPR

	Save the results in csv file format:
	cog, ipr, ji, jc1, jc2, intersect, union, totalprotein_cog, totalprotein_ipr

5) Concat all the results files into a unique file
`cat <OUTPUT_FILE_DIR>/* > COGs_results`




