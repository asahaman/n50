from calculating_n50_assemblies import \
    input_dir_check, output_dir_check, input_dir_empty_check, \
    fasta_extn_check, contigs_len_calc, assembly_len_calc, \
    n50_calc, n50_summary
import sys
import os
import gzip
import numpy as np
import pytest
import pathlib

START_DIR = os.getcwd()
TEST_INPUT_DIR = START_DIR+"/data"
FAKE_INPUT_DIR = START_DIR+"/does_not_exist"
EMPTY_DIR = TEST_INPUT_DIR+'/empty_dir'

list1 = ['test_assem1.fasta','test_assem2.fasta','test_assem3.fasta','test_assem4.fasta','test_assem5.fasta']
list2 = [
    {'contig1': 8, 'contig2': 5, 'contig3': 4, 'contig4': 1, 'contig5': 7, 'contig6': 2, 'contig7': 3, 'contig8': 10, 'contig9': 6, 'contig10': 9},
    {'contig1': 1, 'contig2': 2, 'contig3': 3, 'contig4': 4, 'contig5': 5},
    {'contig1': 8, 'contig2': 5, 'contig3': 150},
    {'contig2': 5, 'contig3': 150, 'contig4': 12},
    {'contig2': 5, 'contig3': 150}
]
tup1 = (
    [7, 4, 150, 150, 150], 
    [0.845, 0.602, 2.176, 2.176, 2.176]
)
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
    for i in range(1,6):
        filename = 'test_assem' + str(i) + '.fasta'
        true_file_list.append(filename)
    expected_list1 = true_file_list
    expected_list2 = ['this','is','not','a','file','list']
    method_list = fasta_extn_check(TEST_INPUT_DIR)
    assert expected_list1 == method_list
    assert expected_list2 != method_list

def test_assembly_len_calc():
    expected_assembly =  list2
    method_assembly = assembly_len_calc(TEST_INPUT_DIR,list1)
    assert expected_assembly == method_assembly

def test_n50_calc():
    expected_tuple = tup1
    method_tuple = n50_calc(list2)
    assert expected_tuple == method_tuple    

