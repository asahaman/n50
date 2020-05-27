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

# Importing functions from fetch_contigs.py 
from n50.fetch_contigs import \
    entrez_api, urls_array, \
    download_needed, url_download

from pathlib import Path
import pytest
import tempfile
import urllib
import os
import sys
import mock
import Bio
import urllib

START_DIR = Path(__file__).parent.absolute()

# Defining parameters for valid NCBI api call
test_genus_good = 'Escherichia'
test_species_good = 'coli'
test_num_good = 3
test_email = 'idontexist@gmail.com'
test_api_key = None

# Defining biopython example return objects from Entrez API calls
TEST_RECORD1 = {'IdList':['7143261', '7143111', '7140871']}
TEST_RECORD2 = {
    'DocumentSummarySet':{
        'DocumentSummary': [ {
            'FtpPath_GenBank':
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/112/445/GCA_013112445.1_ASM1311244v1'
        },
        {
            'FtpPath_GenBank':
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/106/615/GCA_013106615.1_PDT000734164.1'
        },
        {
            'FtpPath_GenBank':
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/106/435/GCA_013106435.1_PDT000733736.1'
        }
        ]
    }
}

TEST_URL_LIST = []
for i in range(0, test_num_good):
    TEST_URL_LIST.append(TEST_RECORD2['DocumentSummarySet']
        ['DocumentSummary'][i]['FtpPath_GenBank']    
    )

TEST_FAKE_URL = 'ftp://fantastic_organism_genomic.fna.gz'

TEST_DOWNLOAD_STR1 = [
    '1st file downloaded',
    '2nd file downloaded',
    '3rd file downloaded'
]

TEST_DOWNLOAD_STR2 = [
    '1st file already exists',
    '2nd file already exists',
    '3rd file already exists'
]


# Defining mockheader for fake url request
class MockHeader:
    headers = {"Content-Length": 0}

    def read(self):
        return ""


class TestFetchContigs:
    @classmethod
    def setup_class(cls):
        sys.path.append(os.path.dirname(os.path.abspath(__file__)))


    # Test for valid API call using mock function
    @mock.patch('Bio.Entrez.esearch', mock.MagicMock(
        return_value='esearch_works_well'))
    @mock.patch('Bio.Entrez.read', mock.MagicMock(
        return_value=TEST_RECORD1))

    def test_entrez_api_good(self):
        method_id_list = entrez_api(test_genus_good, test_species_good,
                                    test_num_good, test_email, test_api_key)
        assert TEST_RECORD1['IdList'] == method_id_list


    # Test for bad API call by altering genus, species and number of contigs
    def test_entrez_api_empty(self):
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


    # Test for generating valid url arrays by mocking Entrez esummary and read 
    @mock.patch('Bio.Entrez.esummary', mock.MagicMock(
        return_value='esummary_works_well'))
    @mock.patch('Bio.Entrez.read', mock.MagicMock(
        return_value=TEST_RECORD2))
    def test_urls_array_good(self):
        method_urls_array = urls_array(TEST_RECORD1['IdList'], test_num_good)
        assert TEST_URL_LIST == method_urls_array


    # Error test for generating empty url arrays
    def test_urls_array_empty(self):
        empty_url_list = []
        with pytest.raises(RuntimeError):
            assert urls_array(empty_url_list, test_num_good)


    # Creating fake url request using monkeypatch
    def test_download_needed_is_needed(self, monkeypatch):
        def mock_open_header(req):
            return MockHeader
        monkeypatch.setattr(urllib.request, "urlopen", mock_open_header)

        # Creating need to download fake url
        with tempfile.TemporaryDirectory() as temp:
            assert download_needed(TEST_FAKE_URL, temp)


    # Test for 'download not needed' since files exist
    def test_download_needed_not_needed(self):
        '''Picking up one of the url links, generating file
        and path objects and evaluating condition where file
        already exists'''

        test_file = os.path.basename(TEST_URL_LIST[0]) + '_genomic.fna.gz'
        test_url_file = os.path.join(TEST_URL_LIST[0], test_file)
        test_file_path = '{}/data/fetch_contigs_test/{}'.format(
            START_DIR, test_file
        )
        assert not download_needed(test_url_file, test_file_path)


    # Test for downloading contigs returning list of strings
    @mock.patch('urllib.request.urlretrieve', mock.MagicMock(
        return_value=TEST_DOWNLOAD_STR1))

    def test_url_download(self):

        '''Downloading contigs to temporary directory returning string object'''
        with tempfile.TemporaryDirectory() as temp:
            assert TEST_DOWNLOAD_STR1 == url_download(TEST_URL_LIST, temp)


    # Test for downloading contigs where files already exist
    @mock.patch('urllib.request.urlretrieve', mock.MagicMock(
        return_value=TEST_DOWNLOAD_STR2))
    def test_url_download_exists(self):

        '''Return string object of existing files'''
        expected_path = '{}/data/fetch_contigs_test/'.format(START_DIR)
        assert TEST_DOWNLOAD_STR2 == url_download(TEST_URL_LIST, expected_path)
