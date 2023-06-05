'''
Copyright:

University of Manitoba & National Microbiology Laboratory, Canada, 2020

Written by: Arnab Saha Mandal

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
'''

import re
import gzip
import operator
import sys
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Declaring global variables based on zipped or unzipped files
GZ_TRUE = 1
GZ_FALSE = 0


# Checking files for valid fasta extension, passing to dictionary
def fasta_extn_check(in_dir=None):
    f_name_dic = {}
    fasta_ext = ['fasta', 'fna', 'ffn', 'faa', 'frn']

    '''Iterating through file list'''
    for f_name in os.listdir(in_dir):

        '''splitting file names by . to capture extension'''
        array_f = f_name.split('.')
        try:

            '''Checking for zipped fasta file. Append file name
            and value to dictionary'''
            if (array_f[-2] in fasta_ext and array_f[-1] == 'gz'):
                f_name_dic[f_name] = GZ_TRUE

                '''Checking for unzipped fasta file. Append file name
                and value to dictionary'''
            elif (array_f[-1] in fasta_ext):
                f_name_dic[f_name] = GZ_FALSE

            '''Raising exception if no fasta files exist'''
        except IndexError:
            pass

    return f_name_dic


# Calculating contig lengths and passing to dictionary
def contigs_len_calc(open_file=None, contig_dict_ok=None):
    contig_dict = {}

    '''Filehandle traversing through every line to capture
    either contig header or nucleotide sequences via regex'''
    for line in open_file:
        x = re.match(r'^>(.+?)\s.*', line)
        y = re.match(r'^[ATGCN-]+$', line)

        if x is not None:
            contig = x.group(1)

            '''Contig header defined as dictionary key, its value
            initialized to zero'''
            contig_dict[contig] = 0

        elif y is not None:
            try:
                '''Length of nucleotide stretches added as value
                for contig header key'''
                contig_dict[contig] += len(y[0].rstrip('\n'))

                '''Ignoring scenario where contig header and value
                cannot be defined'''
            except UnboundLocalError:
                pass

    '''Checking for non zero contig lengths and
    returning as new dictionary'''
    for key, value in contig_dict.items():
        if value != 0:
            contig_dict_ok[key] = value

    if contig_dict_ok:
        return contig_dict_ok


# Appending contig length dictionary per assembly file to a list
def assembly_len_calc(in_dir=None, f_name_dic=None):
    assem_len_list = []

    '''Iterating through input files dictionary'''
    for file_entry in f_name_dic:
        file = os.path.join(in_dir, file_entry)
        contig_dict_ok = {}

        '''Almost identical operations for zipped/unzipped files'''
        if (f_name_dic[file_entry] == GZ_TRUE):
            with gzip.open(file, mode='rt') as open_file:

                '''filehandle uses function to calculate contigs'
                lengths and append the {name:length} to a list'''
                if contigs_len_calc(open_file, contig_dict_ok):

                    assem_len_list.append(
                        contigs_len_calc(open_file, contig_dict_ok)
                    )

        else:
            with open(file, mode='r') as open_file:
                if contigs_len_calc(open_file, contig_dict_ok):
                    assem_len_list.append(
                        contigs_len_calc(open_file, contig_dict_ok)
                    )

    return assem_len_list


# Calculating n50 per assembly file
def n50_calc(assem_len_list=None):
    n50_array = []
    n50_array_log = []

    '''Iterating through dictionaries in assembly length list'''
    for assem in assem_len_list:
        total_assembly_len = 0

        '''Summing up all contig lengths for an assembly'''
        for key in assem:
            total_assembly_len += assem[key]

        '''Defining 50% of total genome length'''
        l50_len = total_assembly_len/2

        '''Reverse sorting contigs by lengths'''
        rev_sorted_contig_dict = sorted(
            assem.items(),
            key=operator.itemgetter(1),
            reverse=True
            )

        '''Adding contig lengths in decreasing order
        till shortest contig at 50% genome length is identified'''
        temp_sum = 0
        for item in rev_sorted_contig_dict:
            temp_sum += item[1]

            if temp_sum >= l50_len:

                '''Contig length identified as N50. Appended to a list'''
                n50_array.append(item[1])

                '''log base 10 of N50 rounded to 3 digits is appended
                to another list'''
                n50_array_log.append(round(math.log(item[1], 10), 3))

                '''Loop exits and moves to next assembly dictionary'''
                break

    return n50_array, n50_array_log


# Calculating n50 statistical summaries
def n50_stat_summary(n50_array=None, out_dir=None, genus=None, species=None):
    output_file = open(os.path.join(out_dir, 'summary_statistics.txt'), 'w')

    '''N50 list is read as pandas dataframe for
    calculation of summary statistics'''
    s = pd.Series(n50_array)
    summary_list = [
        s.count(),
        s.min(), s.max(),
        round(s.mean(), 3),
        s.median(),
        round(s.std(), 3)
    ]

    '''For special cases of 0 and 1 items in N50 list,
    corresponding summary statistic values are defined as NaN's'''
    if s.count() == 0:
        for num in range(1, len(summary_list)):
            summary_list[num] = np.nan
    if s.count() == 1:
        summary_list[-1] = np.nan

    '''Writing to output file in neat sentences'''
    output_file.write("Number of {} {} fasta assemblies is {}\n \
    ".format(genus, species, summary_list[0]))
    output_file.write("Min and Max N50 values are {} and {}\n \
    ".format(summary_list[1], summary_list[2]))
    output_file.write("Mean of N50 distribution is {}\n \
    ".format(summary_list[3]))
    output_file.write("Median of N50 distribution is {}\n \
    ".format(summary_list[4]))
    output_file.write("Standard deviation of N50 distribution \
    is {}\n".format(summary_list[-1]))

    return summary_list


# Plotting log of N50 histogram
def n50_fig_summary(
    n50_array_log=None, out_dir=None, genus=None, species=None
):
    '''Using N50 log as data, passing parameters to matplotlib
    and generating histogram'''
    plt.hist(n50_array_log, bins=100, color='green')
    plt.xlabel('Log (base 10) N50 values')
    plt.ylabel('Counts')
    plt.title('Histogram of ' + genus + ' ' + species + ' assembly lengths')
    plt.savefig(os.path.join(out_dir, 'hist.png'))

    return os.path.join(out_dir, 'hist.png')


# Checking if the user input organism name is ok
def valid_organism_check(genus=None, species=None):
    '''Using regular expressions to define simple checks
    for genus and species combinations'''
    genus_chk = re.match(r'^[A-Z]([a-z]+$|$|\.$)', genus)
    sp_chk = re.match(r'^[a-z]+$', species)

    if genus_chk is None or sp_chk is None:
        return False

    else:
        return True


def main():
    '''Initiating main function with mandatory number of parameters = 4'''
    if len(sys.argv) != 5:
        sys.exit('''
        Command usage: contig_n50 INPUT_DIRECTORY OUTPUT_DIRECTORY GENUS
        SPECIES. Need to pass 4 arguments corresponding to input directory
        containing fasta assembly files, custom output directory, genus and
        species name of the organism
        ''')

    in_dir = sys.argv[1]
    out_dir = sys.argv[2]
    genus = sys.argv[3]
    species = sys.argv[4]

    if not os.path.exists(in_dir):
        sys.exit("Incorrect input directory..exiting..")

    '''Creating dictionary of file names and zipped/unzippped status'''
    file_dic = fasta_extn_check(in_dir)

    '''Creating list of dictionaries of contig lengths'''
    assem_len_list = assembly_len_calc(in_dir, file_dic)

    '''Creating tuple of lists of N50 values and their logs'''
    n50, n50_log = n50_calc(assem_len_list)

    '''Program exits if input organism is not valid'''
    if not valid_organism_check(genus, species):
        sys.exit('''
        User input organism name is incorrect.
        Either input complete scientific name (e.g. Homo sapiens) or
        closely related abbreviation (e.g. H sapiens, H. sapiens are
        fine too). For an approximately complete list of organism names,
        visit NCBI's genome list catalogue at
        https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/
        ''')

    '''Creating output directory if doesn't exist'''
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    '''Calculating N50 summary statistics and generating histogram'''
    n50_stat_summary(n50, out_dir, genus, species)
    n50_fig_summary(n50_log, out_dir, genus, species)

    sys.stdout.write(
        """
        Hooray! you have calculated n50 summary statistics
        and plotted a histogram..
        """
    )


if __name__ == "__main__":
    main()
