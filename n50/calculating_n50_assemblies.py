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

def fasta_extn_check(in_dir = None):
    f_name_list = []
    fasta_ext = ['fasta','fna','ffn','faa','frn']
    for f_name in os.listdir(in_dir):
        if fnmatch.fnmatch(f_name, '*.*'):
            array_f = f_name.split('.')
            if (array_f[-1] in fasta_ext or array_f[-2] in fasta_ext):
                f_name_list.append(f_name)
    return f_name_list

def contigs_len_calc(open_file = None, contig_dict_ok = None):
    contig_dict = {}
    for line in open_file:
        x = re.match(r'^>(.+?)\s.*', line)
        y = re.match(r'^[ATGCN-]+$', line)
        if x != None:
            contig = x.group(1)
            contig_dict[contig] = 0
        elif y != None:
            try:
                contig_dict[contig] += len(y[0].rstrip('\n'))
            except UnboundLocalError:
                '''
                This is a case of orphaned nucleotide sequences
                without a legit > fasta header. Those sequences
                will be ignored and not added to contig dictionary
                '''
                pass
    for key, value in contig_dict.items():
        if value != 0:
            contig_dict_ok[key] = value
    if contig_dict_ok:
        return contig_dict_ok

def assembly_len_calc(in_dir = None, f_name_list = None):
    ass_len_list = []
    for file_entry in f_name_list:
        file = os.path.join(in_dir,file_entry)
        array_f = file_entry.split('.')
        contig_dict_ok = {}
        if (array_f[-1] == 'gz'):
            with gzip.open(file, mode = 'rt') as open_file:
                if contigs_len_calc(open_file, contig_dict_ok):
                    ass_len_list.append(contigs_len_calc(open_file, contig_dict_ok)) 
        else:
            with open(file, mode = 'r') as open_file:
                if contigs_len_calc(open_file, contig_dict_ok):
                    ass_len_list.append(contigs_len_calc(open_file, contig_dict_ok))
    return ass_len_list
            
def n50_calc(ass_len_list = None):
    n50_array = []
    n50_array_log = []
    for ass in ass_len_list:
        total_assembly_len = 0
        for key in ass:
            total_assembly_len += ass[key]
        l50_len = total_assembly_len/2
        rev_sorted_contig_dict = sorted(ass.items(), key=operator.itemgetter(1),reverse=True)
        temp_sum = 0
        for item in rev_sorted_contig_dict:
            temp_sum += item[1]
            if temp_sum >= l50_len: 
                n50_array.append(item[1])
                n50_array_log.append(round(math.log(item[1],10), 3))
                break
    return n50_array, n50_array_log

def n50_summary(n50_array = None, n50_array_log = None, out_dir = None, genus = None, species= None):
    
    output_file = open(os.path.join(out_dir,'summary_statistics.txt'),'w')
    s = pd.Series(n50_array)
    summary_list = [ 
        s.count(), 
        s.min(), s.max(), 
        round(s.mean(),3),
        s.median(), 
        round(s.std(),3)
    ]
    
    output_file.write("Number of assemblies is {}\n".format(s.count()))
    output_file.write("Min and Max N50 values are {} and {}\n".format(s.min(),s.max()))
    output_file.write("Mean of N50 distribution is {}\n".format(round(s.mean(),3)))
    output_file.write("Median of N50 distribution is {}\n".format(s.median()))
    output_file.write("Standard deviation of N50 distribution is {}\n".format(round(s.std(),3)))

    if plt.hist(n50_array_log, bins = 100, color='green'):
        plt.xlabel('Log (base 10) N50 values')
        plt.ylabel('Counts')
        plt.title('Histogram of '+ genus + ' '+ species + ' assembly lengths')
        plt.savefig(os.path.join(out_dir,'hist.pdf'))
        summary_list.append('histogram '+ genus + ' ' + species + ' generated')
     
    return summary_list

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
            if fasta_extn_check(in_dir):
                print("Input directory has fasta files. We can proceed...")
                file_list = fasta_extn_check(in_dir)
                print(file_list)
                ass_len_list = assembly_len_calc(in_dir,file_list)
                if ass_len_list:
                    print(ass_len_list)
                    print("At least one fasta assembly has valid header/s and non-zero nucleotide sequences")
                    arr1, arr2 = n50_calc(ass_len_list)
                    print(arr1, arr2)
                    print("n50 lengths of assemblies have been calculated")
                    if output_dir_check(out_dir):
                        print ("Output dir exists, files can be written")
                    else:
                        print ("Output dir does not exist. new created...")
                    n50_summary(arr1, arr2, out_dir, genus, species)
                    print(n50_summary(arr1, arr2, out_dir, genus, species))
                    print ("Hooray! you have calculated n50 summary statistics and plotted a histogram..")
                else:
                    sys.exit("No fasta file is valid. Absent headers or nucleoties..exiting..")
            else:
                sys.exit("There isnt any fasta file in input directory to work with..exiting..")
        else: 
            sys.exit("Input directory is empty..exiting..")
    else:
        sys.exit ("Incorrect input directory..exiting..")

if __name__ == "__main__":
    main()