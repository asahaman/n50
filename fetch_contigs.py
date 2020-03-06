import sys
#import Entrez library from biopython, letting ncbi know who I am
from Bio import Entrez

def main(email=None, api_key=None):
	Entrez.email = email
	Entrez.api_key = api_key

	#use esearch function to retrieve primary IDs
	#by default, returns 20 items which can be adjusted using "retmax" parameter
	#I set it to 25 as example (will set to 1000 later) and generate an xml object
	handle = Entrez.esearch(db="assembly", term = "E. coli[Organism] AND contig[Assembly Level]", retmax = 25)

	#read function parses XML results to a python object (dictionary)
	record = Entrez.read(handle)



	for i in record['IdList']:
		handle1 = Entrez.esummary(db="assembly", id = i)
		record1 = Entrez.read(handle1)
		j = record1['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
		print(j)

		#j.filepart = j.strip(/) last element of list
		#j.string =  join(j, j.filepart,"_genomic.fna.gz")
		#wget j.string
		#gunzip j.string
		#The contigs within a particular assembly are within the GCA*_genomic.fna file

if __name__ == "__main__":
	print("In Main")
	email = sys.argv[1]
	api_key = sys.argv[2]

	print(email)
	print(api_key)

	main(email, api_key)