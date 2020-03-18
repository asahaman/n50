import sys
import os
import urllib.request
from Bio import Entrez

#I am now heading to exercise two 
#I have a pseudocode at the moment which I would ideally like to implement in python
"""
while reading zipped assembly file line by line:
	if line matches "^>":
		word after that (e.g. RKZR02000001.1) = contig_id
		length_contig_id = 0
		append contig_id to contig_list
	else:
	#or alternatively elif line matches "^[ATGC]":
		split line by characters (nucleotides)
		length_contig_id += length of number of characters

for each element in contig_list:
	print contig_id , length_contig_id
	append length_contig_id to length_contig_list



#N50 calculation

sum_contig_lengths = 0

for each element in length_contig_list:
	sum_contig_lengths += length_contig_id

L50_length = sum_contig_lengths/2
temp_sum = 0
reverse sort length_contig_list (decreasing order, longest first)

for n (0..end_of_length_contig_list):
	temp_sum += length_contig_list[n]
	if (temp_sum >= L50_length):
		print "N50 = length_contig_list[n]"
"""

def main(email = None, api_key = None, genus = None, species = None, retmax = None):
    Entrez.email = email
    Entrez.api_key = api_key
    handle = Entrez.esearch(db= "assembly", term= genus+" "+species+"[Organism] AND contig[Assembly Level]", retmax = num)
    record = Entrez.read(handle)

    if len(record['IdList']) == 0:
        sys.exit("Incorrect organism name. Check for spelling, indentation.")

    else:
        handle1 = Entrez.esummary(db="assembly", id = ",".join(record['IdList']))
        record1 = Entrez.read(handle1)
        assembly_count = 0

        for i in range(0,int(num)):
            url = record1['DocumentSummarySet']['DocumentSummary'][i]['FtpPath_GenBank']
            file_name = os.path.basename(url)
            link = os.path.join(url,file_name+'_genomic.fna.gz')
            print(f'{file_name}.fna.gz' +" downloaded")
            assembly_count += 1
            urllib.request.urlretrieve(link, f'{file_name}.fna.gz')

        print ("Total {} assemblies downloaded".format(assembly_count))            

if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit("""
Command usage: python fetch_contigs.py EMAIL NCBI_API_KEY GENUS_NAME SPECIES_NAME NUM_CONTIGS
Need to pass 5 arguments corresponding to your email, your ncbi api key, a genus and species name (e.g. Salmonella enterica) and the number of contig assemblies you want to download.
If you are not sure how to obtain an NCBI api key, please refer to the README document. 
""")
    
    email = sys.argv[1]
    api_key = sys.argv[2]
    genus = sys.argv[3]
    species = sys.argv[4]
    num = sys.argv[5]

    main(email, api_key, genus, species, num)
