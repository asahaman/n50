import re
import gzip
import operator
import sys
import os

def main(in_dir = None, out_dir = None):
    if os.path.exists(in_dir) == False: 
        sys.exit("Input directory or path incorrect. Exiting..")
    if os.path.exists(out_dir) == False: 
        os.mkdir(out_dir)
    global file_o
    file_o = open(os.path.join(out_dir,'N50.txt'),'w')
    fasta_ext = ['fasta','fna','ffn','faa','frn']

    switch = 'off'
    global f
    for f_name in os.listdir(in_dir):
        if f_name.endswith('.gz'):
            array_f = f_name.split('.')
            if array_f[-2] in fasta_ext:
                file = os.path.join(in_dir,f_name)
                with gzip.open(file, mode = 'rt') as f:
                    n50_calc()
                switch = 'on'
        else:
            array_f = f_name.split('.')
            if array_f[-1] in fasta_ext:
                file = os.path.join(in_dir,f_name)
                with open(file, mode = 'r') as f:
                    n50_calc()
                switch = 'on'

    if switch == 'off':
        sys.exit("No zipped or unzipped fasta assemblies are in the directory. Exiting..")

def n50_calc():
        contig_length_dict = {}
        for line in f:
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
                file_o.write("{}\t{}\n".format(assembly_name, item[1]))
                break

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("""
Command usage: python calculating_n50_assemblies.py INPUT_DIRECTORY OUTPUT_DIRECTORY
Need to pass 2 arguments corresponding to input directory containing fasta assembly files
and output directory where output file "N50.txt" is written.
""")
    in_dir = sys.argv[1]; out_dir = sys.argv[2]
    main(in_dir, out_dir)