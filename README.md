# Py-Wpg
Platform for assembly, annotation and visualization of bacterial genomes

**Exercise 1**: A python (with embedded biopython) script to query NCBI to download E coli assemblies in fasta format. 
Five arguments need to be specified for running the script: your valid email address, your NCBI api key, scientific annotation of the organism (genus and species name, e.g. Salmonella enterica) and number of assemblies you want to download.

If you dont have an NCBI api key, create an account at https://www.ncbi.nlm.nih.gov/account/ using your email. Once you sign in, go to the top right corner where your email id displayed. Click on it and you will be directed to a new page where you can find your api key. Copy the key to use it as an argument for the script. 

Usage: python fetch_contigs.py "Valid email id" "NCBI api_key" "Genus name" "Species Name" "Number of assemblies" 


**Exercise 2**: Python script for analyzing fasta assemblies and calculating N50 contig length

Given a directory of zipped (.gz) assembly files, the script reads each assembly file and embedded fasta formatted contigs. The script subsequently calculates contig lengths and assembly specific N50 values. 

Usage: python calculating_n50_assemblies.py

The program produces two output files. "*summary_concise.txt*" and "*summary_verbose.txt*". *summary_concise.txt* prints assembly name with N50 values:

N50 length for assembly RSZF02 is 91963   
N50 length for assembly RTBZ02 is 71579   
N50 length for assembly RTFI02 is 93615   
.......................................   
.......................................   

*summary_verbose.txt* prints file name, assembly name and every contig's name and length in the following form:  

Opened file GCA_004160775.2_PDT000319866.2.fna.gz  
File has assembly named AAAAOK02 with following contigs and respective lengths:   
	AAAAOK020000001.1       237814   
	AAAAOK020000002.1       174736   
	AAAAOK020000003.1       165434   
	..............................   
	..............................   
Assembly AAAAOK02 has 268 contigs   

This information is printed for each assembly file   
