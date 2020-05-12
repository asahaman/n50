import sys
import os
import urllib.request
from Bio import Entrez


def entrez_api(
    genus=None, species=None, num=None, email=None, api_key=None
):
    """
    Making Entrez API calls to NCBI with valid parameters of
    organism name, number of contigs, email and api key.
    Returning a list of valid assembly IDs.
    """
    Entrez.email = email
    Entrez.api_key = api_key
    handle = Entrez.esearch(
        db="assembly",
        term=genus+" "+species+"[Organism] AND contig[Assembly Level]",
        retmax=num)
    record = Entrez.read(handle)
    id_list = record['IdList']
    return id_list


# Using assembly ID list to generate array of downloadable assemblies
def urls_array(id_list=None, num=None):
    handle = Entrez.esummary(db="assembly", id=",".join(id_list))
    record = Entrez.read(handle)
    urls_list = []
    for i in range(0, int(num)):
        url = (
            record['DocumentSummarySet']['DocumentSummary']
            [i]['FtpPath_GenBank']
        )
        urls_list.append(url)
    return urls_list


# Generating yes/no download switch for every entry in urls list
# Switch returns true if file does not exist or size does not match
# Switch returns false otherwise
def download_needed(full_file_path=None, file_size=None):
    if not os.path.exists(full_file_path) or \
      file_size != os.path.getsize(full_file_path):
        return True
    else:
        return False


# Downloading fasta assemblies specified by url array
def url_download(urls_list=None, out_dir=None):
    for i in range(0, len(urls_list)):
        file_name = os.path.basename(urls_list[i]) + '_genomic.fna.gz'
        url_link = os.path.join(urls_list[i], file_name)
        full_file_path = os.path.join(out_dir, file_name)
        url_open = urllib.request.urlopen(url_link)
        file_size = int(url_open.headers['Content-Length'])
        file_count = i + 1
        if file_count == 1:
            rank_appendix = 'st'
        elif file_count == 2:
            rank_appendix = 'nd'
        elif file_count == 3:
            rank_appendix = 'rd'
        elif file_count > 3:
            rank_appendix = 'th'
        if download_needed(full_file_path, file_size):
            urllib.request.urlretrieve(url_link, full_file_path)
            print(str(file_count) + rank_appendix + ' file downloaded')
        else:
            print(str(file_count) + rank_appendix + ' file already exists')


def main():
    if len(sys.argv) != 7:
        sys.exit("""
        Command usage: contig_dnld EMAIL NCBI_API_KEY GENUS_NAME SPECIES_NAME
        NUM_CONTIGS OUTPUT_DIRECTORY. Need to pass 6 arguments corresponding
        to your email, your ncbi api key, a genus and species
        name (e.g. Salmonella enterica), the number of contig assemblies
        you want to download and an output directory to download
        the assembly files. If you are not sure how to obtain an NCBI api key,
        please refer to the README document.
        """)
    else:
        email = sys.argv[1]
        api_key = sys.argv[2]
        genus = sys.argv[3]
        species = sys.argv[4]
        num = sys.argv[5]
        out_dir = sys.argv[6]
        api_list = entrez_api(genus, species, num, email, api_key)
        if not api_list:
            sys.exit("""
            Incorrect organism name or ordering of arguments.
            Check for spelling, indentation of the organism and
            refer to README for correct ordering of arguments.
            """)
        else:
            urls_list = urls_array(api_list, num)
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            url_download(urls_list, out_dir)


if __name__ == "__main__":
    main()
