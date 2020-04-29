from calculating_n50_assemblies import \
    input_dir_check, output_dir_check, input_dir_empty_check, \
    fasta_extn_check, contigs_len_calc, assembly_len_calc, \
    n50_calc, n50_summary
import sys
import os
import gzip
import numpy as np
import pytest

START_DIR = os.getcwd()
TEST_INPUT_DIR = START_DIR+'/data'
FAKE_INPUT_DIR = START_DIR+'/does_not_exist'
EMPTY_DIR = TEST_INPUT_DIR+'/empty_dir'

'''
Different examples of test fasta assembly files have been created in folder data:
1. test_assem1.fasta - 10 contigs, valid file
2. test_assem2.fasta - 5 contigs, valid file
3. test_assem3.fasta - 6 valid contig headers:
                            2 have valid sequences, 
                            1 has mix and match of valid sequence and invalid sequence
                            3 are invalid contigs:
                                1 has no sequence
                                1 has blanks
                                1 has random text paragraph
4. test_assem4.fasta - 3 invalid contigs:
                            1 has no sequence
                            1 has blanks
                            1 has random text paragraph
5. test_assem5.fasta - Random text paragraph but will pass as a valid assembly due to the nature of embedded text
6. test_assem6.fasta - Random text paragraph which won't pass
7. test_assem7.fasta - Orphaned nucleotides (valid sequences without fasta headers). This won't pass                            
'''

test_file_list = [
    'test_assem1.fasta','test_assem2.fasta','test_assem3.fasta',
    'test_assem4.fasta','test_assem5.fasta','test_assem6.fasta', 
    'test_assem7.fasta'
]

test_ass_list = [
    {'contig1': 8, 'contig2': 5, 'contig3': 4, 'contig4': 1, 'contig5': 7, 
    'contig6': 2, 'contig7': 3, 'contig8': 10, 'contig9': 6, 'contig10': 9}, 
    {'contig1': 1, 'contig2': 2, 'contig3': 3, 'contig4': 4, 'contig5': 5}, 
    {'contig1': 5, 'contig2': 150, 'contig3': 12}, 
    {'I': 29}
]

test_n50_tuple = (
    [7, 4, 150, 29],
    [0.845, 0.602, 2.176, 1.462]
)

test_summary_list = [4, 4, 150, 47.5, 18.0, 69.236, 'histogram desired organism generated']

def test_input_dir_check_true():
    expected_result = True
    method_result = input_dir_check(TEST_INPUT_DIR)
    assert expected_result == method_result

def test_input_dir_check_false():
    expected_result = False
    method_result = input_dir_check(FAKE_INPUT_DIR)
    assert expected_result == method_result

def test_input_dir_empty_check_true():
    expected_result = True
    method_result = input_dir_empty_check(TEST_INPUT_DIR)
    assert expected_result == method_result

def test_input_dir_empty_check_false():
    expected_result = False
    method_result = input_dir_empty_check(EMPTY_DIR)
    assert expected_result == method_result

def test_fasta_extn_check_true():
    true_file_list = []
    for i in range(1,8):
        filename = 'test_assem' + str(i) + '.fasta'
        true_file_list.append(filename)
    method_list = fasta_extn_check(TEST_INPUT_DIR)
    assert true_file_list == method_list

def test_fasta_extn_check_false():
    false_file_list = ['this','is','not','a','file','list']
    method_list = fasta_extn_check(TEST_INPUT_DIR)
    assert false_file_list != method_list

def test_assembly_len_calc_true():
    expected_assembly = test_ass_list
    method_assembly = assembly_len_calc(TEST_INPUT_DIR,test_ass_list)
    assert expected_assembly == method_assembly

def test_assembly_len_calc_false():
    test_ass_list_false = test_ass_list[1:]
    method_assembly = assembly_len_calc(TEST_INPUT_DIR,test_ass_list)
    assert test_ass_list_false != method_assembly

def test_n50_calc():
    expected_tuple = test_n50_tuple
    method_tuple = n50_calc(test_ass_list)
    assert expected_tuple == method_tuple
print(n50_summary(test_n50_tuple[0],test_n50_tuple[1],TEST_INPUT_DIR,'A','Histogram'))

def test_n50_summary():
    expected_list = test_summary_list
    method_list = n50_summary(test_n50_tuple[0],test_n50_tuple[1],TEST_INPUT_DIR,'desired','organism')
    assert expected_list == method_list

