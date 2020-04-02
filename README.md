# Py-Wpg
Platform for assembly, annotation and visualization of bacterial genomes

**Exercise 1**: A python (with embedded biopython) script to query NCBI to download E coli assemblies in fasta format. 
Six arguments need to be specified for running the script: your valid email address, your NCBI api key, scientific annotation of the organism (genus and species name, e.g. Salmonella enterica), the number of assemblies you want to download and an output directory to download the assembly files.

If you dont have an NCBI api key, create an account at https://www.ncbi.nlm.nih.gov/account/ using your email. Once you sign in, go to the top right corner where your email id displayed. Click on it and you will be directed to a new page where you can find your api key. Copy the key to use it as an argument for the script. 

Usage: python fetch_contigs.py "Valid email id" "NCBI api_key" "Genus name" "Species Name" "Number of assemblies" "Output directory"  
Output/s: zipped fasta assembly files inside the output directory

**Exercise 2**: Python script for analyzing fasta assemblies and calculating N50 lengths, summary statistics and plotting their distribution  

Given a directory of fasta assembly files, the script reads each assembly file and embedded fasta formatted contigs to calculate N50 values. Subsequently, calculates summary statistics using *pandas.Series.Describe* function and generates a histogram of N50 assembly lengths using *matplotlib* 

Usage: python calculating_n50_assemblies.py INPUT_DIRECTORY OUTPUT_DIRECTORY  
Output/s: *summary_statistics.txt* file containing summary statistics and *hist.pdf* plot showing histogram   

