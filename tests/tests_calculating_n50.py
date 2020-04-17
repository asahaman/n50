from n50.calculating_n50_assemblies import n50_calc, n50_summary
import sys
import os
import fnmatch
import numpy as np

n50_expected = [7,4]

n50_obtained = []
n50_obtained_log = []
for f_name in os.listdir():
    if fnmatch.fnmatch(f_name, '*.fasta'):
        file = os.path.join(f_name)
        with open(file, mode = 'r') as open_file:
            n50_calc(open_file, n50_obtained, n50_obtained_log)

assert n50_expected == n50_obtained