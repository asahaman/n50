from n50.fetch_contigs import \
    entrez_api, urls_array, \
    download_needed, url_download

from Bio import Entrez
from pathlib import Path
import pytest
import tempfile
import urllib

START_DIR = Path(__file__).parent.absolute()
TEST_OUTPUT_DIR = '{}/output/'.format(str(START_DIR))


test_genus = 'Escherichia'
test_species = 'coli'
test_num = 3
test_email = 'idontexist@gmail.com'
test_api_key = None

TEST_ID_LIST = ['7097491', '7097481', '7097471']
TEST_URL_LIST = [
    'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/041/555/GCA_013041555.1_PDT000409044.1',
    'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/041/545/GCA_013041545.1_PDT000409479.1',
    'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/041/525/GCA_013041525.1_PDT000409464.1'
]

TEST_FAKE_URL = 'ftp://ftp.nicegenome.com/fantastic_organism/fantastic_organism_genomic.fna.gz'

class MockHeader:
    headers = {"Content-Length": 0}
    def read(self):
        return ""

def test_entrez_api_good():
    method_id_list = entrez_api(test_genus, test_species, test_num, test_email, test_api_key)
    assert TEST_ID_LIST == method_id_list


def test_entrez_api_empty():
    expected_id_list = []
    method_id_list = entrez_api(test_genus, test_species, 0, test_email, test_api_key)
    assert expected_id_list == method_id_list


def test_urls_array_good():
    method_urls_array = urls_array(TEST_ID_LIST, test_num)
    assert TEST_URL_LIST == method_urls_array


def test_urls_array_empty():
    empty_url_list = []
    with pytest.raises(RuntimeError):
        assert urls_array(empty_url_list, test_num)

def test_download_needed(monkeypatch):
    def mock_open_header(req):
        return MockHeader
    monkeypatch.setattr(urllib.request, "urlopen", mock_open_header)
    with tempfile.TemporaryDirectory() as temp:
        assert download_needed(TEST_FAKE_URL,temp)