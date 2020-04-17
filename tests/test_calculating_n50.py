from n50.calculating_n50_assemblies import n50_calc, n50_summary
import sys
import os
import fnmatch
import gzip
import numpy as np

in_dir = sys.argv[1] 
out_dir = sys.argv[2]

n50_expected = [7,4]

if os.path.exists(in_dir) == False: 
    sys.exit("Input directory or path incorrect. Exiting..")
if os.path.exists(out_dir) == False: 
    os.mkdir(out_dir)
output_file = open(os.path.join(out_dir,'summary_statistics.txt'),'w')
fasta_ext = ['fasta','fna','ffn','faa','frn']
n50_obtained = [] 
n50_obtained_log = []

def test_n50_calc():
    for f_name in os.listdir(in_dir):
        if fnmatch.fnmatch(f_name, '*.*'):
            array_f = f_name.split('.')
            file = os.path.join(in_dir,f_name)
            if (array_f[-1] == 'gz' and array_f[-2] in fasta_ext):
                with gzip.open(file, mode = 'rt') as open_file:
                    n50_calc(open_file, n50_obtained, n50_obtained_log)
            elif array_f[-1] in fasta_ext:
                with open(file, mode = 'r') as open_file:
                    n50_calc(open_file, n50_obtained, n50_obtained_log)
    
    assert np.array_equal(n50_expected,n50_obtained)