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

# Importing functions from calculating_n50_assemblies.py
from n50.calculating_n50_assemblies import \
    fasta_extn_check, contigs_len_calc, assembly_len_calc, \
    n50_calc, n50_stat_summary, n50_fig_summary, \
    valid_organism_check

import sys
import os
import numpy as np
import filecmp
from pathlib import Path

START_DIR = Path(__file__).parent.absolute()
TEST_INPUT_DIR = '{}/data/'.format(str(START_DIR))
EMPTY_DIR = '{}/empty_dir/'.format(str(TEST_INPUT_DIR))

'''
Different examples of test fasta assembly files
have been created in folder data:
1. test_assem1.fasta - 10 contigs, valid file

2. test_assem2.fasta - 5 contigs, valid file

3. test_assem3.fasta - 6 valid contig headers:
                            2 have valid sequences,
                            1 has mix and match of valid sequence
                                and invalid sequence
                            3 are invalid contigs:
                                1 has no sequence
                                1 has blanks
                                1 has random text paragraph

4. test_assem4.fasta - 3 invalid contigs:
                            1 has no sequence
                            1 has blanks
                            1 has random text paragraph

5. test_assem5.fasta - Random text paragraph but will pass
        as a valid assembly due to the nature of embedded text

6. test_assem6.fasta - Random text paragraph which won't pass

7. test_assem7.fasta - Orphaned nucleotides (valid sequences
        without fasta headers). This won't pass

8. test_assem8.fasta - Empty file (with a valid fasta extension)
'''

# Defining return objects for successful test functions'''
 
'''Dictionary of filenames with value of 0/1 indicating unzipped/zipped'''
test_file_dic_good = {
    'test_assem1.fasta': 0, 'test_assem2.fasta': 0, 'test_assem3.fasta': 0,
    'test_assem4.fasta': 0, 'test_assem5.fasta': 0, 'test_assem6.fasta': 0,
    'test_assem7.fasta': 0, 'test_assem8.fasta': 0
}

'''List of dictionaries. Each dictionary corresponds to one assembly file
with contig names (keys) and lengths (values)'''
test_assem_list_good = [
    {'contig1': 8, 'contig2': 5, 'contig3': 4, 'contig4': 1, 'contig5': 7,
     'contig6': 2, 'contig7': 3, 'contig8': 10, 'contig9': 6, 'contig10': 9},
    {'contig1': 1, 'contig2': 2, 'contig3': 3, 'contig4': 4, 'contig5': 5},
    {'contig1': 5, 'contig2': 150, 'contig3': 12},
    {'I': 1}
]

'''Tuple of N50 values and their corresponding logarithms with base 10.
Each tuple element is a list of N50 values for all valid assembly files'''
test_n50_tuple_good = (
    [7, 4, 150, 1],
    [0.845, 0.602, 2.176, 0.0]
)

'''List of summary statistics of N50 values for all assembly files'''
test_n50_stat_summary_good = [4, 1, 150, 40.5, 5.5, 73.041]

'''Defining pre-existing histogram plot of log N50 values'''
expected_fig = os.path.join(TEST_INPUT_DIR, 'expected.png')

'''Creating good and bad combinations of organism names'''
genus_good = ['Homo', 'H.', 'H']
species_good = ['sapiens', 'sap']
genus_bad = ['homo', 'H,', 'Homo[Organism]']
species_bad = ['sap[Org]', 'sap.', 'Sapiens']


class TestCalcN50:
    @classmethod
    def setup_class(cls):
        sys.path.append(os.path.dirname(os.path.abspath(__file__)))


    # Test for checking good fasta file name extensions
    def test_fasta_extn_check_true(self):
        test_file_dic = {}
        
        '''Creating dictionary of 8 fasta file names with value 0'''
        for i in range(1, 9):
            filename = 'test_assem' + str(i) + '.fasta'
            test_file_dic[filename] = 0
        
        method_dic = fasta_extn_check(TEST_INPUT_DIR)
        assert test_file_dic == method_dic


    # Test for checking fasta files in empty directory
    def test_fasta_extn_check_empty(self):
        test_file_dic = {}
        method_dic = fasta_extn_check(EMPTY_DIR)
        assert test_file_dic == method_dic


    # Test for checking multiple contig lengths
    def test_contigs_len_calc_multiple(self):
        '''First element of list of contigs dictionary set to a variable'''  
        test_contig_dic = test_assem_list_good[0]
        
        empty_dicn = {}
        file = os.path.join(TEST_INPUT_DIR, 'test_assem1.fasta')
        
        '''Function reads first assembly file and outputs
        contig names and lengths dictionary, which is compared'''
        with open(file, mode='r') as open_file:
            method_contigs_dic = contigs_len_calc(open_file, empty_dicn)
            assert test_contig_dic == method_contigs_dic


    # Test for checking single contig lengths
    def test_contigs_len_calc_single(self):
        '''Last element of list of contigs dictionary set to a variable'''
        test_contig_dic = test_assem_list_good[-1]
        
        empty_dicn = {}
        file = os.path.join(TEST_INPUT_DIR, 'test_assem5.fasta')
        
        '''Function reads 5th assembly file and outputs
        contig names (just one) and length dictionary, which is compared'''
        with open(file, mode='r') as open_file:
            method_contigs_dic = contigs_len_calc(open_file, empty_dicn)
            assert test_contig_dic == method_contigs_dic


    # Test for creating dictionary of multiple contig lengths
    def test_assembly_len_calc_multiple(self):
        expected_assembly = test_assem_list_good
        
        '''Calling function to generate list of
        contig lengths dictionaries, which is compared'''
        method_assembly = assembly_len_calc(TEST_INPUT_DIR, test_file_dic_good)
        assert expected_assembly == method_assembly


    # Test for creating dictionary of single contig length
    def test_assembly_len_calc_single(self):
        '''Isolating a file with only one contig'''
        test_file_dic_single = {'test_assem5.fasta': 0}

        '''Resulting list of dictionaries will contain one element,
        which is compared'''
        expected_assembly = [test_assem_list_good[-1]]
        method_assembly = assembly_len_calc(TEST_INPUT_DIR, test_file_dic_single)
        assert expected_assembly == method_assembly


    # Test for creating dictionary of non-exiting contig lengths
    def test_assembly_len_calc_empty(self):
        '''Declaring empty (contig lengths)
        dictionary and empty list'''
        test_file_dic_empty = {}
        expected_assembly = []

        '''Function with empty contig lengths dictionary will return empty list,
        which is compared''' 
        method_assembly = assembly_len_calc(EMPTY_DIR, test_file_dic_empty)
        assert expected_assembly == method_assembly


    # Test for calculating n50 values of multiple contigs
    def test_n50_calc_multiple(self):
        expected_tuple = test_n50_tuple_good
        
        '''Calling function to generate tuple of
        N50 values, which is compared'''
        method_tuple = n50_calc(test_assem_list_good)
        assert expected_tuple == method_tuple


    # Test for calculating n50 values of single contig
    def test_n50_calc_single(self):
        '''Defining tuple for one N50 value and its corresponding log'''
        expected_tuple = ([test_n50_tuple_good[0][3]], [test_n50_tuple_good[1][3]])
        
        '''Calling function will generate tuple of single N50 value
        and its log, which is compared'''
        method_tuple = n50_calc([test_assem_list_good[-1]])
        assert expected_tuple == method_tuple


    # Test for calculating n50 values of non-existent contigs
    def test_n50_calc_empty(self):
        '''Defining empty tuples of N50 values and its log'''
        expected_tuple = ([], [])
        
        '''Calling function will generate empty tuples of lists'''
        method_tuple = n50_calc([])
        assert expected_tuple == method_tuple


    # Test for valid organism check
    def test_valid_organism_check_good(self):
        '''Iterating through passable genus and species names
        and pairing them'''
        for g in range(0, len(genus_good)):
            for s in range(0, len(species_good)):

                '''Calling function returns True for every acceptable pair
                of genus and species'''
                assert valid_organism_check(genus_good[g], species_good[s])


    # Test for invalid organism check
    def test_valid_organism_check_bad(self):
        '''Iterating through bad combinations of genus and species names
        and pairing them'''
        for g in range(0, len(genus_bad)):
            for s in range(0, len(species_bad)):

                '''Calling function returns False for every pair
                of genus and species'''
                assert not valid_organism_check(genus_bad[g], species_bad[s])


    # Test for checking n50 summary statistics
    def test_n50_stat_summary_multiple(self):
        expected_list = test_n50_stat_summary_good
        
        '''Calling function returns list of summary
        statistics of N50 values, which is compared'''  
        method_list = n50_stat_summary(test_n50_tuple_good[0],
                                    TEST_INPUT_DIR, 'desired', 'organism')
        assert expected_list == method_list


    # Test for checking n50 summary statistics for single file
    def test_n50_summary_single(self):
        '''Defining N50 summary statistics for only 5th file'''
        expected_list = [1, 1, 1, 1.0, 1.0, np.nan]
        
        '''Calling function will return N50 summary statistics for 5th
        file, which is compared'''
        method_list = n50_stat_summary(test_n50_tuple_good[0][3],
                                    TEST_INPUT_DIR, 'desired', 'organism')
        assert expected_list == method_list


    # Test for checking n50 summary statistics for no file
    def test_n50_summary_empty(self):
        '''Defining N50 summary statistics for no file, hence first
        element (count) is 0, rest are NaN's'''
        expected_list = [0]
        expected_list += 5 * [np.nan]

        '''Calling function will return list of 0 and NaN's which are
        compared'''
        test_empty_assembly_list = []
        method_list = n50_stat_summary(test_empty_assembly_list,
                                    TEST_INPUT_DIR, 'desired', 'organism')
        assert expected_list == method_list


    # Test for checking histogram plot
    def test_n50_fig_summary(self):
        '''Calling function will generate matplotlib histogram for a
        given list of log N50 values, which is compared by filecmp'''
        method_fig = n50_fig_summary(test_n50_tuple_good[1],
                                    TEST_INPUT_DIR, 'desired', 'organism')
        assert filecmp.cmp(expected_fig, method_fig)
