import re
import gzip
import operator
import sys
import os
import statistics
import numpy as np
import pandas as pd

def main(in_dir = None, out_dir = None):
    if os.path.exists(in_dir) == False: 
        sys.exit("Input directory or path incorrect. Exiting..")
    if os.path.exists(out_dir) == False: 
        os.mkdir(out_dir)
    output_file1 = open(os.path.join(out_dir,'N50.txt'),'w')
    output_file2 = open(os.path.join(out_dir,'summary_statistics.txt'),'w')
    fasta_ext = ['fasta','fna','ffn','faa','frn']
    n50_array = []
    switch = 'off'
    for f_name in os.listdir(in_dir):
        if f_name.endswith('.gz'):
            array_f = f_name.split('.')
            if array_f[-2] in fasta_ext:
                file = os.path.join(in_dir,f_name)
                with gzip.open(file, mode = 'rt') as open_file:
                    n50_calc(open_file, output_file1, n50_array)
                switch = 'on'
        else:
            array_f = f_name.split('.')
            if array_f[-1] in fasta_ext:
                file = os.path.join(in_dir,f_name)
                with open(file, mode = 'r') as open_file:
                    n50_calc(open_file, output_file1, n50_array)
                switch = 'on'
    
    if switch == 'on':
        n50_summary(n50_array, output_file2)

    else:
        sys.exit("No zipped or unzipped fasta assemblies are in the directory. Exiting..")
    

def n50_calc(open_file = None, output_file1 = None, n50_array = None):
    contig_length_dict = {}
    for line in open_file:
        assembly = re.findall(r'^>([A-Z]+?0\d).+',line)
        if len(assembly) > 0:
            assembly_name = assembly[0] 
        x = re.findall(r'^>(.+?)\s', line)
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
            output_file1.write("{}\t{}\n".format(assembly_name, item[1]))
            break

def n50_summary(n50_array = None, output_file2 = None):
    s = pd.Series(n50_array)
    output_file2.write("Number of N50 data points is {}\n".format(s.count()))
    output_file2.write("Min and Max N50 values are {} and {}\n".format(s.min(),s.max()))
    output_file2.write("Mean of N50 distribution is {:.3f}\n".format(s.mean()))
    output_file2.write("Median of N50 distribution is {}\n".format(s.median()))
    output_file2.write("Standard deviation of N50 distribution is {:.3f}\n".format(s.std()))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("""
Command usage: python calculating_n50_assemblies.py INPUT_DIRECTORY OUTPUT_DIRECTORY
Need to pass 2 arguments corresponding to input directory containing fasta assembly files
and output directory where output file "N50.txt" is written.
""")
    in_dir = sys.argv[1]; out_dir = sys.argv[2]
    main(in_dir, out_dir)