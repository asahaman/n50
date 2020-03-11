import sys
import os
import urllib.request
import gzip
import re
from Bio import Entrez

def main(email=None, api_key=None):
    Entrez.email = email
    Entrez.api_key = api_key

    handle = Entrez.esearch(db="assembly", term = "E. coli[Organism] AND contig[Assembly Level]", retmax = 1)

    record = Entrez.read(handle)

    for i in record['IdList']:
        handle1 = Entrez.esummary(db="assembly", id = i)
        record1 = Entrez.read(handle1)
        url = record1['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        file_name = os.path.basename(url)
        link = os.path.join(url,file_name+'_genomic.fna.gz')
        print(link)
        urllib.request.urlretrieve(link, f'{file_name}.fna.gz')
        with gzip.open(f'{file_name}.fna.gz','r') as f:
            x = re.search("^>", f)
            print (x)
        
if __name__ == "__main__":
    print("In Main")
    email = sys.argv[1]
    api_key = sys.argv[2]

    print(email)
    print(api_key)

    main(email, api_key)
