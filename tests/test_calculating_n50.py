from n50.calculating_n50_assemblies import \
    fasta_extn_check, contigs_len_calc, assembly_len_calc, \
    n50_calc, n50_stat_summary, n50_fig_summary, \
    valid_organism_check

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

test_file_dic_good = {
    'test_assem1.fasta': 0, 'test_assem2.fasta': 0, 'test_assem3.fasta': 0,
    'test_assem4.fasta': 0, 'test_assem5.fasta': 0, 'test_assem6.fasta': 0,
    'test_assem7.fasta': 0, 'test_assem8.fasta': 0
}

test_assem_list_good = [
    {'contig1': 8, 'contig2': 5, 'contig3': 4, 'contig4': 1, 'contig5': 7,
     'contig6': 2, 'contig7': 3, 'contig8': 10, 'contig9': 6, 'contig10': 9},
    {'contig1': 1, 'contig2': 2, 'contig3': 3, 'contig4': 4, 'contig5': 5},
    {'contig1': 5, 'contig2': 150, 'contig3': 12},
    {'I': 1}
]

test_n50_tuple_good = (
    [7, 4, 150, 1],
    [0.845, 0.602, 2.176, 0.0]
)

test_n50_stat_summary_good = [4, 1, 150, 40.5, 5.5, 73.041]

genus_good = ['Homo', 'H.', 'H']
species_good = ['sapiens', 'sap']
genus_bad = ['homo', 'H,', 'Homo[Organism]']
species_bad = ['sap[Org]', 'sap.', 'Sapiens']


# Tests for checking fasta file name extensions
def test_fasta_extn_check_true():
    test_file_dic = {}
    for i in range(1, 9):
        filename = 'test_assem' + str(i) + '.fasta'
        test_file_dic[filename] = 0
    method_dic = fasta_extn_check(TEST_INPUT_DIR)
    assert test_file_dic == method_dic


def test_fasta_extn_check_empty():
    test_file_dic = {}
    method_dic = fasta_extn_check(EMPTY_DIR)
    assert test_file_dic == method_dic


# Tests for checking multiple and single contig lengths
def test_contigs_len_calc_multiple():
    test_contig_dic = test_assem_list_good[0]
    empty_dicn = {}
    file = os.path.join(TEST_INPUT_DIR, 'test_assem1.fasta')
    with open(file, mode='r') as open_file:
        method_contigs_dic = contigs_len_calc(open_file, empty_dicn)
        assert test_contig_dic == method_contigs_dic


def test_contigs_len_calc_single():
    test_contig_dic = test_assem_list_good[-1]
    empty_dicn = {}
    file = os.path.join(TEST_INPUT_DIR, 'test_assem5.fasta')
    with open(file, mode='r') as open_file:
        method_contigs_dic = contigs_len_calc(open_file, empty_dicn)
        assert test_contig_dic == method_contigs_dic


# Tests for creating dictionary of contig lengths
def test_assembly_len_calc_multiple():
    expected_assembly = test_assem_list_good
    method_assembly = assembly_len_calc(TEST_INPUT_DIR, test_file_dic_good)
    assert expected_assembly == method_assembly


def test_assembly_len_calc_single():
    test_file_dic_single = {'test_assem5.fasta': 0}
    expected_assembly = [test_assem_list_good[-1]]
    method_assembly = assembly_len_calc(TEST_INPUT_DIR, test_file_dic_single)
    assert expected_assembly == method_assembly


def test_assembly_len_calc_empty():
    test_file_dic_empty = {}
    expected_assembly = []
    method_assembly = assembly_len_calc(EMPTY_DIR, test_file_dic_empty)
    assert expected_assembly == method_assembly


# Tests for calculating n50 values
def test_n50_calc_multiple():
    expected_tuple = test_n50_tuple_good
    method_tuple = n50_calc(test_assem_list_good)
    assert expected_tuple == method_tuple


def test_n50_calc_single():
    expected_tuple = ([test_n50_tuple_good[0][3]], [test_n50_tuple_good[1][3]])
    method_tuple = n50_calc([test_assem_list_good[-1]])
    assert expected_tuple == method_tuple


def test_n50_calc_empty():
    expected_tuple = ([], [])
    method_tuple = n50_calc([])
    assert expected_tuple == method_tuple


# Tests for organism check
def test_valid_organism_check_good():
    for g in range(0, len(genus_good)):
        for s in range(0, len(species_good)):
            assert valid_organism_check(genus_good[g], species_good[s])


def test_valid_organism_check_bad():
    for g in range(0, len(genus_bad)):
        for s in range(0, len(species_bad)):
            assert not valid_organism_check(genus_bad[g], species_bad[s])


# Tests for checking n50 summary statistics
def test_n50_stat_summary_multiple():
    expected_list = test_n50_stat_summary_good
    method_list = n50_stat_summary(test_n50_tuple_good[0],
                                   TEST_INPUT_DIR, 'desired', 'organism')
    assert expected_list == method_list


def test_n50_summary_single():
    expected_list = [1, 1, 1, 1.0, 1.0, np.nan]
    method_list = n50_stat_summary(test_n50_tuple_good[0][3],
                                   TEST_INPUT_DIR, 'desired', 'organism')
    assert expected_list == method_list


def test_n50_summary_empty():
    expected_list = [0]
    expected_list += 5 * [np.nan]
    test_empty_assembly_list = []
    method_list = n50_stat_summary(test_empty_assembly_list,
                                   TEST_INPUT_DIR, 'desired', 'organism')
    assert expected_list == method_list


# Test for checking histogram plot
def test_n50_fig_summary():
    expected_fig = os.path.join(TEST_INPUT_DIR, 'expected.png')
    method_fig = n50_fig_summary(test_n50_tuple_good[1],
                                 TEST_INPUT_DIR, 'desired', 'organism')
    assert filecmp.cmp(expected_fig, method_fig)
