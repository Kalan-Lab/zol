import os
import sys
import argparse
from Bio import SeqIO
from collections import defaultdict
from time import sleep
from zol import util
import subprocess
import traceback
import multiprocessing

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
	Program: setup_bigscape.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

    Downloads 
	""", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--install_path', help='Path to where BiG-SCAPE should be installed.', required=False, default=lsaBGC_main_directory + 'external_tools/')
    parser.add_argument('-p', '--pfam_download_path', help='Path to where Pfam should be downloaded and setup.', required=False, default=lsaBGC_main_directory + 'db/')

    args = parser.parse_args()
    return args

def setup_annot_dbs():
    myargs = create_parser()

    bigscape_install_path = os.path.abspath(myargs.install_path) + '/'
    pfam_download_path = os.path.abspath(myargs.pfam_download_path) + '/'

    try:
        assert(os.path.isdir(bigscape_install_path))
    except:
        sys.stderr.write('Error: Provided directory for installing BiG-SCAPE does not exist!')
        sys.exit(1)

    try:
        assert(os.path.isdir(pfam_download_path))
    except:
        sys.stderr.write('Error: Provided directory for downloading Pfam files does not exist!\n')
        sys.exit(1)

    # download files in requested directory
    try:
        os.chdir(pfam_download_path)

        pfam_hmm_file = pfam_download_path + 'Pfam-A.hmm'
        if not os.path.isfile(pfam_hmm_file):
            # Download Pfam HMM
            print('Setting up Pfam database!')
            os.system('wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz')
            os.system('gunzip Pfam-A.hmm.gz')

            pfam_hmm = pfam_download_path + 'Pfam-A.hmm'
            assert(os.path.isfile(pfam_hmm_file))
            os.system('hmmpress ' + pfam_hmm)

    except:
        sys.stderr.write('Error: issues with downloading or seting up Pfam database files! Please post to Github issues if unable to figure out!\n')
        sys.exit(1)

    try:
        bigscape_prog = bigscape_install_path + 'BiG-SCAPE-1.1.4/bigscape.py'
        if not os.path.isfile(bigscape_prog):
            os.chdir(bigscape_install_path)
            os.system('wget https://github.com/medema-group/BiG-SCAPE/archive/refs/tags/v1.1.4.tar.gz')
            os.system('tar -zxvf v1.1.4.tar.gz')
            assert(os.path.isfile(bigscape_prog))
    except:
        sys.stderr.write('Error: issues with installing BiG-SCAPE. Please post to Github issues if unable to figure out what the issue is!\n')
        sys.exit(1)

    if os.path.isfile(bigscape_prog) and os.path.isfile(pfam_hmm_file):
        bigscape_info_file = open(lsaBGC_main_directory + 'external_tools/bigscape_location.txt', 'w')
        bigscape_info_file.write(bigscape_prog + '\t' + pfam_download_path + '\n')
        bigscape_info_file.close()

setup_annot_dbs()