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



