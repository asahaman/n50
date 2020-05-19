from n50.fetch_contigs import \
    entrez_api, urls_array, \
    download_needed, url_download

from pathlib import Path
import pytest
import tempfile
import urllib
import os

START_DIR = Path(__file__).parent.absolute()
test_genus_good = 'Escherichia'
test_species_good = 'coli'
test_num_good = 3
test_email = 'idontexist@gmail.com'
test_api_key = None

TEST_ID_LIST = ['7143111', '7140871', '7140861']
TEST_URL_LIST = [
    'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/106/615/GCA_013106615.1_PDT000734164.1',
    'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/106/435/GCA_013106435.1_PDT000733736.1',
    'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/106/395/GCA_013106395.1_PDT000734354.1'
]

TEST_FAKE_URL = 'ftp://fantastic_organism_genomic.fna.gz'


class MockHeader:
    headers = {"Content-Length": 0}

    def read(self):
        return ""


# Tests for good and bad API calls
def test_entrez_api_good():
    method_id_list = entrez_api(test_genus_good, test_species_good,
                                test_num_good, test_email, test_api_key)
    assert TEST_ID_LIST == method_id_list


def test_entrez_api_empty():
    test_num_bad = 0
    with pytest.raises(RuntimeError):
        assert entrez_api(test_genus_good, test_species_good,
                          test_num_bad, test_email, test_api_key)
    test_genus_bad = 'Ebola'
    with pytest.raises(RuntimeError):
        assert entrez_api(test_genus_bad, test_species_good,
                          test_num_good, test_email, test_api_key)
    test_species_bad = 'colus'
    with pytest.raises(RuntimeError):
        assert entrez_api(test_genus_good, test_species_bad,
                          test_num_good, test_email, test_api_key)


# Tests for list of url arrays
def test_urls_array_good():
    method_urls_array = urls_array(TEST_ID_LIST, test_num_good)
    assert TEST_URL_LIST == method_urls_array


def test_urls_array_empty():
    empty_url_list = []
    with pytest.raises(RuntimeError):
        assert urls_array(empty_url_list, test_num_good)


# Tests for download needed
def test_download_needed_is_needed(monkeypatch):
    def mock_open_header(req):
        return MockHeader
    monkeypatch.setattr(urllib.request, "urlopen", mock_open_header)
    with tempfile.TemporaryDirectory() as temp:
        assert download_needed(TEST_FAKE_URL, temp)


def test_download_needed_not_needed():
    test_file = os.path.basename(TEST_URL_LIST[0]) + '_genomic.fna.gz'
    test_url_file = os.path.join(TEST_URL_LIST[0], test_file)
    test_file_path = '{}/data/fetch_contigs_test/{}'.format(
        START_DIR, test_file
    )
    assert not download_needed(test_url_file, test_file_path)


# Tests for downloading contigs
def test_url_download():
    expected_string_list = [
        '1st file downloaded',
        '2nd file downloaded',
        '3rd file downloaded'
    ]
    assert expected_string_list == url_download(TEST_URL_LIST, START_DIR)


def test_url_download_exists():
    expected_string_list = [
        '1st file already exists',
        '2nd file already exists',
        '3rd file already exists'
    ]
    expected_path = '{}/data/fetch_contigs_test/'.format(START_DIR)
    assert expected_string_list == url_download(TEST_URL_LIST, expected_path)
