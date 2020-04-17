import re
import gzip
import operator
import sys
import os
import fnmatch
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) != 5:
        sys.exit("""
Command usage: contig_n50 INPUT_DIRECTORY OUTPUT_DIRECTORY GENUS SPECIES
Need to pass 4 arguments corresponding to input directory containing fasta assembly files, custom output directory and genus and species name of the organism.
""")

    else:
        in_dir = sys.argv[1] 
        out_dir = sys.argv[2]
        genus = sys.argv[3]
        species = sys.argv[4]
        
        if os.path.exists(in_dir) == False: 
            sys.exit("Input directory or path incorrect. Exiting..")
        if os.path.exists(out_dir) == False: 
            os.mkdir(out_dir)
        output_file = open(os.path.join(out_dir,'summary_statistics.txt'),'w')
        fasta_ext = ['fasta','fna','ffn','faa','frn']
        n50_array = [] 
        n50_array_log = []
        for f_name in os.listdir(in_dir):
            if fnmatch.fnmatch(f_name, '*.*'):
                array_f = f_name.split('.')
                file = os.path.join(in_dir,f_name)
                if (array_f[-1] == 'gz' and array_f[-2] in fasta_ext):
                    with gzip.open(file, mode = 'rt') as open_file:
                        n50_calc(open_file, n50_array, n50_array_log)
                elif array_f[-1] in fasta_ext:
                    with open(file, mode = 'r') as open_file:
                        n50_calc(open_file, n50_array, n50_array_log)
                    
        if not n50_array:
            sys.exit("No zipped or unzipped fasta assemblies are in the directory. Exiting..")
                    
        else:
            n50_summary(n50_array, n50_array_log, output_file, out_dir, genus, species)

def n50_calc(open_file = None, n50_array = None, n50_array_log = None):
    contig_length_dict = {}
    for line in open_file:
        x = re.findall(r'>(.+?)\s.+', line)
        if len(x) > 0:
            contig_name = x[0]
            contig_length = 0
        else:
            contig_length += len(line.rstrip('\n'))
            contig_length_dict[contig_name] = contig_length

    total_assembly_len = 0
    for key in contig_length_dict:
        total_assembly_len += contig_length_dict[key]
    l50_len = total_assembly_len/2
    rev_sorted_contig_length_dict = sorted(contig_length_dict.items(), key=operator.itemgetter(1),reverse=True)
    temp_sum = 0
    for item in rev_sorted_contig_length_dict:
        temp_sum += item[1]
        if temp_sum >= l50_len: 
            n50_array.append(item[1])
            n50_array_log.append(math.log(item[1],10))
            break

def n50_summary(n50_array = None, n50_array_log = None, output_file = None, out_dir = None, genus = None, species = None):    
    s = pd.Series(n50_array)
    output_file.write("Number of assemblies is {}\n".format(s.count()))
    output_file.write("Min and Max N50 values are {} and {}\n".format(s.min(),s.max()))
    output_file.write("Mean of N50 distribution is {:.3f}\n".format(s.mean()))
    output_file.write("Median of N50 distribution is {}\n".format(s.median()))
    output_file.write("Standard deviation of N50 distribution is {:.3f}\n".format(s.std()))
    plt.hist(n50_array_log, bins = 100, color='green')
    plt.xlabel('Log (base 10) N50 values')
    plt.ylabel('Counts')
    plt.title('Histogram of '+ genus + ' '+ species + ' assembly lengths')
    plt.savefig(os.path.join(out_dir,'hist.pdf'))

if __name__ == "__main__":
    main()