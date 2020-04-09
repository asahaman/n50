import sys
import fetch_contigs
import calculating_n50_assemblies

def main(email=None, api_key=None, genus=None, species=None, num=None, out1_dir=None, out2_dir=None):
    fetch_contigs.main(email, api_key, genus, species, num, out1_dir)
    calculating_n50_assemblies.main(out1_dir, out2_dir, genus, species)

if __name__ == "__main__":
    if len(sys.argv) != 8:
        sys.exit("""
Need 7 arguments: email, api key, genus, species, number of contigs, 
directory 1 (output of fetch_contigs), directory 2 (output of calculating_n50_assemblies). Directory
1 and 2 can be same
""")
    email = sys.argv[1]
    api_key = sys.argv[2]
    genus = sys.argv[3]
    species = sys.argv[4]
    num = sys.argv[5]
    out1_dir = sys.argv[6]
    out2_dir = sys.argv[7]
    main(email, api_key, genus, species, num, out1_dir, out2_dir)