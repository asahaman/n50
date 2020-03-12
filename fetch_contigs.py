import sys
import os
import urllib.request
from Bio import Entrez

def main(email=None, api_key=None, retmax=None):
    Entrez.email = email
    Entrez.api_key = api_key

    handle = Entrez.esearch(db="assembly", term = "E. coli[Organism] AND contig[Assembly Level]", retmax = num)

    record = Entrez.read(handle)
    assembly_count = 0
    for i in record['IdList']:
        handle1 = Entrez.esummary(db="assembly", id = i)
        record1 = Entrez.read(handle1)
        url = record1['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        file_name = os.path.basename(url)
        link = os.path.join(url,file_name+'_genomic.fna.gz')
        print(f'{file_name}.fna.gz' +" downloaded")
        assembly_count += 1
        urllib.request.urlretrieve(link, f'{file_name}.fna.gz')

    print ("Total {} assemblies downloaded".format(assembly_count))            

if __name__ == "__main__":
    print("In Main")
    email = sys.argv[1]
    api_key = sys.argv[2]
    num = sys.argv[3]

    print(email)
    print(api_key)

    main(email, api_key, num)
