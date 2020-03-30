# Py-Wpg
Platform for assembly, annotation and visualization of bacterial genomes

**Exercise 1**: A python (with embedded biopython) script to query NCBI to download E coli assemblies in fasta format. 
Five arguments need to be specified for running the script: your valid email address, your NCBI api key, scientific annotation of the organism (genus and species name, e.g. Salmonella enterica) and number of assemblies you want to download.

If you dont have an NCBI api key, create an account at https://www.ncbi.nlm.nih.gov/account/ using your email. Once you sign in, go to the top right corner where your email id displayed. Click on it and you will be directed to a new page where you can find your api key. Copy the key to use it as an argument for the script. 

Usage: python fetch_contigs.py "Valid email id" "NCBI api_key" "Genus name" "Species Name" "Number of assemblies" 


**Exercise 2**: Python script for analyzing fasta assemblies and calculating N50 contig length

Given a directory of fasta assembly files, the script reads each assembly file and embedded fasta formatted contigs to calculate N50 values. 

Usage: python calculating_n50_assemblies.py INPUT_DIRECTORY OUTPUT_DIRECTORY  

The program needs two arguments corresponding to input directory containing fasta assembly files and output directory where output file "N50.txt" is written with assembly name and N50 lengths in tab separated value format:   
   
RSZF02  91963   
RTBZ02  71579   
RTFI02  93615   
..............
.............. 
