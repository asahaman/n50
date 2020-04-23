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

def input_dir_check(in_dir = None):
    if os.path.exists(in_dir):
        return True
    else:
        return False
    
def output_dir_check(out_dir = None):
    if os.path.exists(out_dir): 
        return True
    else:
        os.mkdir(out_dir)
        return False

def input_dir_empty_check(in_dir = None):
    if os.listdir(in_dir):
        return True
    else:
        return False

def fasta_file_check(in_dir = None):
    f_name_list = []
    fasta_ext = ['fasta','fna','ffn','faa','frn']
    for f_name in os.listdir(in_dir):
        if fnmatch.fnmatch(f_name, '*.*'):
            array_f = f_name.split('.')
            if (array_f[-1] in fasta_ext or array_f[-2] in fasta_ext):
                f_name_list.append(f_name)
    return f_name_list

def n50_wrapper(in_dir = None, f_name_list = None):
    n50_array = []
    n50_array_log = []
    for file_entry in f_name_list:
        file = os.path.join(in_dir,file_entry)
        array_f = file_entry.split('.')
        if (array_f[-1] == 'gz'):
            with gzip.open(file, mode = 'rt') as open_file:
                #if good_fasta_check(open_file):
                    n50_calc(open_file, n50_array, n50_array_log)
        else:
            with open(file, mode = 'r') as open_file:
                #if good_fasta_check(open_file):
                    n50_calc(open_file, n50_array, n50_array_log)
    return n50_array, n50_array_log

def good_fasta_check(open_file = None):
    fasta_check = []
    for line in open_file:
        if re.search(r'^>.+', line) != None or re.search(r'^[ATGCN-]+$', line) != None:
            pass
        else:
            fasta_check.append(1)
    if not fasta_check:
        return True
    else:
        return False

def n50_calc(open_file = None, n50_array = None, n50_array_log = None):
    contig_length_dict = {}
    for line in open_file:
        x = re.findall(r'>(.+?)\s.*', line)
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

def n50_summary(n50_array = None, n50_array_log = None, out_dir = None, genus = None, species= None):
    output_file = open(os.path.join(out_dir,'summary_statistics.txt'),'w')
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
    return True

def main():
    if len(sys.argv) != 5:
        sys.exit("""
Command usage: contig_n50 INPUT_DIRECTORY OUTPUT_DIRECTORY GENUS SPECIES
Need to pass 4 arguments corresponding to input directory containing fasta assembly files, custom output directory and genus and species name of the organism.
""")

    in_dir = sys.argv[1]
    out_dir = sys.argv[2]
    genus = sys.argv[3]
    species = sys.argv[4]

    if input_dir_check(in_dir):
        print ("Input directory exists, lets look inside...")
        if input_dir_empty_check(in_dir):
            print ("Input directory has some files. Lets look at the files...")
            if fasta_file_check(in_dir):
                print("Input directory has fasta files. We can proceed...")
                file_list = fasta_file_check(in_dir)
                if output_dir_check(out_dir):
                    print ("Output dir exists, files can be written")
                else:
                    print ("Output dir does not exist. new created...")
                if file_list:
                    arr1, arr2 = n50_wrapper(in_dir,file_list)
                    if (arr1, arr2):
                        n50_summary(arr1, arr2, out_dir, genus, species)
                        if n50_summary(arr1, arr2, out_dir, genus, species):
                            print ("Hooray! you have calculated n50 summary statistics and plotted a histogram..")
                        else:
                            print ("n50 values have been calculated but statistics cant be computed")
                    else:
                        sys.exit("Something is wrong with the fasta files...exiting..")
            else:
                sys.exit("There isnt any fasta file in input directory to work with..exiting..")
        else: 
            sys.exit("Input directory is empty..exiting..")
    else:
        sys.exit ("Incorrect input directory..exiting..")

if __name__ == "__main__":
    main()