# Py-Wpg
Platform for assembly, annotation and visualization of bacterial genomes

Exercise 1: A python (with embedded biopython) script to query NCBI to download E coli assemblies in fasta format.
Five arguments need to be specified for running the script: your valid email address, your NCBI api key, scientific annotation of the organism (genus and species name, e.g. Salmonella enterica) and number of assemblies you want to download.
If you dont have an NCBI api key, create an account at https://www.ncbi.nlm.nih.gov/account/ using your email. Once you sign in, go to the top right corner where your email id displayed. Click on it and you will be directed to a new page where you can find your api key. Copy the key to use it as an argument for the script. 

Usage: python fetch_contigs.py "Valid email id" "NCBI api_key" "Genus name" "Species Name" "Number of assemblies" 
