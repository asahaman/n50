import sys
import os
import urllib.request
from Bio import Entrez

def valid_organism(genus = None, species = None, num = None, email = None, api_key = None):
    Entrez.email = email
    Entrez.api_key = api_key
    handle = Entrez.esearch(db= "assembly", term= genus+" "+species+"[Organism] AND contig[Assembly Level]", retmax = num)
    record_list = Entrez.read(handle)['IdList']
    return record_list

def main():
    if len(sys.argv) != 7:
        sys.exit("""
Command usage: contig_dnld EMAIL NCBI_API_KEY GENUS_NAME SPECIES_NAME NUM_CONTIGS OUTPUT_DIRECTORY
Need to pass 6 arguments corresponding to your email, your ncbi api key, a genus and species name (e.g. Salmonella enterica), the number of contig assemblies you want to download and an output directory to download the assembly files.
If you are not sure how to obtain an NCBI api key, please refer to the README document. 
""")

    else:
        email = sys.argv[1]
        api_key = sys.argv[2]
        genus = sys.argv[3]
        species = sys.argv[4]
        num = sys.argv[5]
        out_dir = sys.argv[6]

        if len(valid_organism(genus, species, num, email, api_key)) == 0:
            sys.exit("""
            Incorrect organism name or ordering of arguments. 
            Check for spelling, indentation of the organism and refer to README for correct ordering of arguments.
            """)

        else:
            if os.path.exists(out_dir) == False:
                os.mkdir(out_dir)
            handle1 = Entrez.esummary(db="assembly", id = ",".join(valid_organism(genus, species, num, email, api_key)))
            record1 = Entrez.read(handle1)
            assembly_count = 0

            for i in range(0,int(num)):
                url = record1['DocumentSummarySet']['DocumentSummary'][i]['FtpPath_GenBank']
                file_name = os.path.basename(url)
                link = os.path.join(url,file_name+'_genomic.fna.gz')
                full_file_path = os.path.join(out_dir,f'{file_name}.fna.gz')
                print(f'{file_name}.fna.gz' +" downloaded")
                assembly_count += 1
                urllib.request.urlretrieve(link, full_file_path)

            print ("Total {} assemblies downloaded to {}".format(assembly_count, out_dir))

if __name__ == "__main__":
    main()