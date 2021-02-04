# JACCARD index calculation between COGs and integrated InterPro entries

The process is divided into 5 different steps

1) The representative list of proteins from UNIPROT is loaded in a redis hash (`get_rp_proteins.py`)

2) For each protein in the representative list (`IPR_RP_AC.py`):
	For each integrated InterPro entry, we calculated the midpoint of the different fragments and count the number of protein matching
    
	All the data are saved in redis hashes:
	- protein_hash[protein_acc]={"IPR1":[midpoint1,midpoint2], "IPR2":[midpoint1]}
	- ipr_count[ipr]=count_protein

3) Get the list of COGs into a redis queue which will be used by 8 different jobs submitted in parallel for step 4 (`get_cog_list.py`)

4) For each COG (`search_intersect.py`):

		For each protein:
	
    	- Get IPR and midpoints
    	- Verify if intersection between COG and IPR

	Save the results in csv file format:
	cog, ipr, ji, jc1, jc2, intersect, union, totalprotein_cog, totalprotein_ipr

5) Concat all the results files into a unique file




