This directory is for the analysis of all SNPs from GWAS Catalog.

Data:
All associations in the GWAS Catalog with added ontology annotations (v1.0.1)

Data url:
https://www.ebi.ac.uk/gwas/api/search/downloads/alternative

Download date:
02-09-2016, 2:29.

Downloaded GWAS associations file: 
data/downloaded/gwas_catalog_v1.0.1-associations_e85_r2016-08-25.tsv

Extraction of SNPs:
script: scripts/extract_gwas_snps.py.  
Output: data/all_gwas.bed (20,783 non-repeating SNPs from associations)
	data/trait_snps (.txt files for mapped traits SNPs)
	data/error_log.txt (unresolved data)

