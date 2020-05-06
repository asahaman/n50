import re
import gzip
import operator
import sys
import os
import math
import pandas as pd
import matplotlib.pyplot as plt

# Declaring global variables based on zipped or unzipped files
GZ_TRUE = 1
GZ_FALSE = 0


# Checking files for valid fasta extension, passing to dictionary
def fasta_extn_check(in_dir=None):
    f_name_dic = {}
    fasta_ext = ['fasta', 'fna', 'ffn', 'faa', 'frn']
    for f_name in os.listdir(in_dir):
        array_f = f_name.split('.')
        try:
            if (array_f[-2] in fasta_ext and array_f[-1] == 'gz'):
                f_name_dic[f_name] = GZ_TRUE
            elif (array_f[-1] in fasta_ext):
                f_name_dic[f_name] = GZ_FALSE
        except IndexError:
            pass
    return f_name_dic


# Calculating contig lengths and passing to dictionary
def contigs_len_calc(open_file=None, contig_dict_ok=None):
    contig_dict = {}
    for line in open_file:
        x = re.match(r'^>(.+?)\s.*', line)
        y = re.match(r'^[ATGCN-]+$', line)
        if x is not None:
            contig = x.group(1)
            contig_dict[contig] = 0
        elif y is not None:
            try:
                contig_dict[contig] += len(y[0].rstrip('\n'))
            except UnboundLocalError:
                pass
    for key, value in contig_dict.items():
        if value != 0:
            contig_dict_ok[key] = value
    if contig_dict_ok:
        return contig_dict_ok


# Appending contig length dictionary per assembly file to a list
def assembly_len_calc(in_dir=None, f_name_dic=None):
    assem_len_list = []
    for file_entry in f_name_dic:
        file = os.path.join(in_dir, file_entry)
        contig_dict_ok = {}
        if (f_name_dic[file_entry] == GZ_TRUE):
            with gzip.open(file, mode='rt') as open_file:
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
    for assem in assem_len_list:
        total_assembly_len = 0
        for key in assem:
            total_assembly_len += assem[key]
        l50_len = total_assembly_len/2
        rev_sorted_contig_dict = sorted(
            assem.items(),
            key=operator.itemgetter(1),
            reverse=True
            )
        temp_sum = 0
        for item in rev_sorted_contig_dict:
            temp_sum += item[1]
            if temp_sum >= l50_len:
                n50_array.append(item[1])
                n50_array_log.append(round(math.log(item[1], 10), 3))
                break
    return n50_array, n50_array_log


# Calculating n50 statistical summaries
def n50_stat_summary(n50_array=None, out_dir=None, genus=None, species=None):
    output_file = open(os.path.join(out_dir, 'summary_statistics.txt'), 'w')
    s = pd.Series(n50_array)
    summary_list = [
        s.count(),
        s.min(), s.max(),
        round(s.mean(), 3),
        s.median(),
        round(s.std(), 3)
    ]
    if s.count() == 0:
        for num in range(1, len(summary_list)):
            summary_list[num] = 'NA'
    if s.count() == 1:
        summary_list[-1] = 'NA'
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


# Plotting n50 histogram
def n50_fig_summary(
    n50_array_log=None, out_dir=None, genus=None, species=None
):
    fig = plt.figure()
    plt.hist(n50_array_log, bins=100, color='green')
    plt.xlabel('Log (base 10) N50 values')
    plt.ylabel('Counts')
    plt.title('Histogram of ' + genus + ' ' + species + ' assembly lengths')
    fig.savefig(os.path.join(out_dir, 'hist.pdf'))
    return fig


# Checking if the user input organism name is ok
def valid_organism_check(genus=None, species=None):
    genus_chk = re.match(r'^[A-Z]([a-z]+$|$|\.$)', genus)
    sp_chk = re.match(r'^[a-z]+$', species)
    if genus_chk is None or sp_chk is None:
        return False
    else:
        return True


def main():
    if len(sys.argv) != 5:
        sys.exit(
            """
            Command usage: contig_n50 INPUT_DIRECTORY OUTPUT_DIRECTORY GENUS
            SPECIES
            Need to pass 4 arguments corresponding to input directory
            containing fasta assembly files, custom output directory, genus and
            species name of the organism
            """
        )
    in_dir = sys.argv[1]
    out_dir = sys.argv[2]
    genus = sys.argv[3]
    species = sys.argv[4]
    if not os.path.exists(in_dir):
        sys.exit("Incorrect input directory..exiting..")
    file_dic = fasta_extn_check(in_dir)
    assem_len_list = assembly_len_calc(in_dir, file_dic)
    n50, n50_log = n50_calc(assem_len_list)
    if not valid_organism_check(genus, species):
        sys.exit(
            """
            User input organism name is incorrect.
            Either input complete scientific name (e.g. Homo sapiens) or
            closely related abbreviation (e.g. H sapiens, H. sapiens are
            fine too). For an approximately complete list of organism names,
            visit NCBI's genome list catalogue at
            https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/
            """
        )
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    n50_stat_summary(n50, out_dir, genus, species)
    print(n50_stat_summary(n50, out_dir, genus, species))
    n50_fig_summary(n50_log, out_dir, genus, species)
    print(
        """
        Hooray! you have calculated n50 summary statistics
        and plotted a histogram..
        """
    )


if __name__ == "__main__":
    main()
