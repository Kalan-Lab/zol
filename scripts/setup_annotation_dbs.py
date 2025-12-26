#!/usr/bin/env python3

"""
Program: setup_annotation_dbs.py
Author: Rauf Salamzade
Kalan Lab
UW Madison, Department of Medical Microbiology and Immunology
"""

# BSD 3-Clause License
#
# Copyright (c) 2023-2025, Kalan-Lab
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
import os
import sys
import argparse
from Bio import SeqIO
import shutil
import traceback

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: setup_annotation_dbs.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
		
	Downloads annotation databases for annotations and lateral transfer inference: KOfam, 
	PGAP, PaperBlast, MIBiG, CARD, Conserved Ribosomal Proteins, VOGs, TnCentral, VFDB,
	MOB-suite MPF and MOB Proteins, and Pfam.
								  
	Location of where to download databases is controlled by setting the environmental
	variable ZOL_DATA_PATH e.g.:
								  							   
	$ export ZOL_DATA_PATH=/path/to/database_location/zol_db/

	CAUTION: This directory will be deleted and recreated when this script is run - it 
	should therefore be unique and specific to zol.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-c', '--threads', type=int, help="Number of threads to use [Default is 4].", required=False, default=4)
	parser.add_argument('-m', '--minimal', action='store_true', help="Minimal mode - will only download Pfam and PGAP HMMs.", required=False, default=False)
	parser.add_argument('-ld', '--lsabgc-minimal', action='store_true', help="Minimal mode for lsaBGC - will only download Pfam HMMs, PGAP HMMs, CARD proteins, & MIBiG proteins.", required=False, default=False)	
	parser.add_argument('-y', '--yes', action='store_true', help="Assume 'yes' for prompts (non-interactive).", required=False, default=False)

	args = parser.parse_args()
	return args

def setup_annot_dbs():
	myargs = create_parser()

	download_path = None
	if str(os.getenv("ZOL_DATA_PATH")) != 'None':
		download_path = str(os.getenv("ZOL_DATA_PATH"))
	if download_path == None or not os.path.isdir(download_path):
		sys.stderr.write('Issues validing database download directory exists.\n')
		sys.exit(1)
	else:
		download_path = os.path.abspath(download_path) + "/"
	
	threads = myargs.threads
	minimal_mode = myargs.minimal
	lsabgc_minimal_mode = myargs.lsabgc_minimal
	non_interactive_yes = myargs.yes or str(os.getenv('ZOL_NONINTERACTIVE', '0')).strip() in ['1', 'true', 'True']

	try:
		assert(os.path.isdir(download_path))
		if not non_interactive_yes:
			response = input(f"The directory {download_path}\nalready exists, will delete it and recreate it. (This directory\nshould be specific to zol not a general directory for\ndatabases) Proceed with deleting? (yes/no): ")
			if response.lower() != 'yes':
				os.system('Deletion not requested! Exiting ...')
				sys.exit(1)
	except Exception as e:
		sys.stderr.write('Error: Provided directory for downloading annotation files does not exist or user did not accept deleting the directory and recreating it!\n')

	if lsabgc_minimal_mode:
		sys.stdout.write('lsaBGC minimal mode requested, will only be downloading the Pfam, MIBiG, CARD, and PGAP databases.\n')
	elif minimal_mode:
		sys.stdout.write('Minimal mode requested, will only be downloading the Pfam and PGAP databases.\n')
	
	try:
		shutil.rmtree(download_path)
		os.mkdir(download_path)
		readme_outf = open(download_path + 'README.txt', 'w')
		readme_outf.write('Default space for downloading KOFam annotation databases.\n')
		readme_outf.close()
	except Exception as e:
		sys.stderr.write('Issues clearing contents of db/ to re-try downloads.\n')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

	# download vScore information and GECCO scores

	sys.stdout.write('Downloading vScore information and GECCO weights ...\n')
	download_links = [
		'https://anantharamanlab.github.io/V-Score-Search/VScoreDataNormalized.csv',
		'https://raw.githubusercontent.com/raufs/gtdb_gca_to_taxa_mappings/refs/heads/main/GECCO_Weights.txt'
	]

	# Download
	print('Starting download of files!')
	os.chdir(download_path)
	try:
		for dl in download_links:
			bf = dl.split('/')[-1].split('?')[0]
			download_output = download_path + bf
			curl_download_dbs_cmd = ['curl', '-L', dl, '-o', download_output]
			os.system(' '.join(curl_download_dbs_cmd))
			assert(os.path.isfile(download_output))
	except Exception as e:
		sys.stderr.write('Error occurred during downloading of weight files - please let us know on GitHub issues!\n')
		sys.stderr.write(str(e) + '\n')

	listing_file = download_path + 'database_location_paths.txt'
	issues_file = download_path + 'issues.txt'
	listing_handle = open(listing_file, 'w')
	issues_handle = open(issues_file, 'w')

	if lsabgc_minimal_mode:
		# Final annotation files
		pgap_info_file = download_path + 'hmm_PGAP.tsv'
		pgap_phmm_file = download_path + 'PGAP.hmm'
		mb_faa_file = download_path + 'mibig.dmnd'
		pfam_phmm_file = download_path + 'Pfam-A.hmm'
		card_faa_file = download_path + 'card.dmnd'

		download_links = ['https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz',
				  'https://dl.secondarymetabolites.org/mibig/mibig_prot_seqs_4.0.fasta',
				  'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv',
				  'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz',
				  'https://card.mcmaster.ca/latest/data']

		# Download
		print('Starting download of files!')
		os.chdir(download_path)
		try:
			for dl in download_links:
				bf = dl.split('/')[-1].split('?')[0]
				download_output = download_path + bf
				if bf == 'data':
					bf = 'card-data.tar.bz2'
					download_output = download_path + bf
				curl_download_dbs_cmd = ['curl', '-L', dl, '-o', download_output]
				os.system(' '.join(curl_download_dbs_cmd))
		except Exception as e:
			sys.stderr.write('Error occurred during downloading of databases!\n')
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up Pfam database!')
			os.system(' '.join(['gunzip', 'Pfam-A.hmm.gz']))
			assert(os.path.isfile(pfam_phmm_file))
			name = None
			desc = None
			pfam_descriptions_file = download_path + 'pfam_descriptions.txt'
			pdf_handle = open(pfam_descriptions_file, 'w')
			with open(pfam_phmm_file) as oppf:
				for line in oppf:
					line = line.strip()
					ls = line.split()
					if ls[0].strip() == 'NAME':
						name = ' '.join(ls[1:]).strip()
					elif ls[0].strip() == 'DESC':
						desc = ' '.join(ls[1:]).strip()
						pdf_handle.write(name + '\t' + desc + '\n') # type: ignore
			pdf_handle.close()
			os.system(' '.join(['hmmpress', pfam_phmm_file]))
			z = 0
			with open(pfam_phmm_file) as oppf:
				for line in oppf:
					if line.startswith('NAME'): z += 1
			listing_handle.write('pfam\t' + pfam_descriptions_file + '\t' + pfam_phmm_file + '\t' + str(z) + '\n')
		except Exception as e:
			sys.stderr.write('Issues setting up Pfam database.\n')
			issues_handle.write('Issues setting up Pfam database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up PGAP database!')
			extract_hmm_dir = 'hmm_PGAP.HMM/'
			os.mkdir(extract_hmm_dir)
			os.system(' '.join(['tar', '-zxf', 'hmm_PGAP.HMM.tgz', '-C', 'hmm_PGAP.HMM/']))
			assert (os.path.isfile(pgap_info_file))
			assert (os.path.isdir(download_path + 'hmm_PGAP.HMM/'))
			for folder, subs, files in os.walk(extract_hmm_dir):
				for filename in files:
					if filename.endswith('.HMM') or filename.endswith('.hmm'):
						hmm_file_path = os.path.abspath(folder) + '/' + filename
						os.system(' '.join(['cat', hmm_file_path, '>>', pgap_phmm_file]))
			pgap_descriptions_file = download_path + 'pgap_descriptions.txt'
			pdf_handle = open(pgap_descriptions_file, 'w')
			with open(pgap_info_file) as opil:
				for i, line in enumerate(opil):
					if i == 0: continue
					line = line.strip()
					ls = line.split('\t')
					label = ls[2]
					description = ls[10]
					pdf_handle.write(label + '\t' + description + '\n')
			pdf_handle.close()
			assert (os.path.isfile(pgap_phmm_file))
			os.system(' '.join(['hmmpress', pgap_phmm_file]))
			z = 0
			with open(pgap_phmm_file) as oppf:
				for line in oppf:
					if line.startswith('NAME'): z += 1

			listing_handle.write('pgap\t' + pgap_descriptions_file + '\t' + pgap_phmm_file + '\t' + str(z) + '\n')
			os.system(' '.join(['rm', '-rf', download_path + 'hmm_PGAP.HMM/', download_path + 'hmm_PGAP.HMM.tgz', pgap_info_file]))
		except Exception as e:
			sys.stderr.write('Issues setting up PGAP database.\n')
			issues_handle.write('Issues setting up PGAP database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up MIBiG database!')
			os.system(' '.join(['diamond', 'makedb', '--in', 'mibig_prot_seqs_4.0.fasta', '-d', mb_faa_file, '--threads', str(threads)]))
			assert(os.path.isfile(mb_faa_file))
			mibig_descriptions_file = download_path + 'mibig_descriptions.txt'
			mdf_handle = open(mibig_descriptions_file, 'w')
			with open(download_path + 'mibig_prot_seqs_4.0.fasta') as omf:
				for rec in SeqIO.parse(omf, 'fasta'):
					mdf_handle.write(rec.id + '\t' + rec.description + '\n')
			mdf_handle.close()
			listing_handle.write('mibig\t' + mibig_descriptions_file + '\t' + mb_faa_file + '\tNA\n')
			os.system(' '.join(['rm', '-rf', download_path + 'mibig_prot_seqs_4.0.fasta']))
		except Exception as e:
			sys.stderr.write('Issues setting up MIBiG database.\n')
			issues_handle.write('Issues setting up MIBiG database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')
			
		try:
			print('Setting up CARD database!')
			print('Check out the ARTS webserver for more detailed analysis in finding BGCs for synthesizing antibiotics.')
			os.mkdir(download_path + 'CARD_DB_Files/')
			os.system(' '.join(['tar', '-xf', download_path + 'card-data.tar.bz2', '-C', download_path + 'CARD_DB_Files/']))
			os.system(' '.join(['mv', download_path + 'CARD_DB_Files/protein_fasta_protein_homolog_model.fasta', download_path]))
			os.system(' '.join(['diamond', 'makedb', '--in', download_path + 'protein_fasta_protein_homolog_model.fasta', '-d', card_faa_file, '--threads', str(threads)]))
			card_descriptions_file = download_path + 'card_descriptions.txt'
			cdf_handle = open(card_descriptions_file, 'w')
			with open(download_path + 'protein_fasta_protein_homolog_model.fasta') as ocf:
				for rec in SeqIO.parse(ocf, 'fasta'):
					cdf_handle.write(rec.id + '\t' + rec.description + '\n')
			cdf_handle.close()
			os.system(' '.join(['rm', '-rf', download_path + 'CARD_DB_Files/', download_path + 'protein_fasta_protein_homolog_model.fasta', download_path + 'card-data.tar.bz2']))
			assert(os.path.isfile(card_faa_file))
			listing_handle.write('card\t' + card_descriptions_file + '\t' + card_faa_file + '\tNA\n')
		except Exception as e:
			sys.stderr.write('Issues setting up CARD database.\n')
			issues_handle.write('Issues setting up CARD database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')


	
	elif minimal_mode:
		# Final annotation files
		pgap_info_file = download_path + 'hmm_PGAP.tsv'
		pgap_phmm_file = download_path + 'PGAP.hmm'
		pfam_phmm_file = download_path + 'Pfam-A.hmm'

		download_links = ['https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz',
						  'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv',
						  'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz']

		# Download
		print('Starting download of files!')
		os.chdir(download_path)
		try:
			for dl in download_links:
				bf = dl.split('/')[-1].split('?')[0]
				download_output = download_path + bf
				curl_download_dbs_cmd = ['curl', '-L', dl, '-o', download_output]
				os.system(' '.join(curl_download_dbs_cmd))
		except Exception as e:
			sys.stderr.write('Error occurred during downloading of databases!\n')
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up Pfam database!')
			os.system(' '.join(['gunzip', 'Pfam-A.hmm.gz']))
			assert(os.path.isfile(pfam_phmm_file))
			name = None
			desc = None
			pfam_descriptions_file = download_path + 'pfam_descriptions.txt'
			pdf_handle = open(pfam_descriptions_file, 'w')
			with open(pfam_phmm_file) as oppf:
				for line in oppf:
					line = line.strip()
					ls = line.split()
					if ls[0].strip() == 'NAME':
						name = ' '.join(ls[1:]).strip()
					elif ls[0].strip() == 'DESC':
						desc = ' '.join(ls[1:]).strip()
						pdf_handle.write(name + '\t' + desc + '\n') # type: ignore
			pdf_handle.close()
			os.system(' '.join(['hmmpress', pfam_phmm_file]))
			z = 0
			with open(pfam_phmm_file) as oppf:
				for line in oppf:
					if line.startswith('NAME'): z += 1
			listing_handle.write('pfam\t' + pfam_descriptions_file + '\t' + pfam_phmm_file + '\t' + str(z) + '\n')
		except Exception as e:
			sys.stderr.write('Issues setting up Pfam database.\n')
			issues_handle.write('Issues setting up Pfam database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up PGAP database!')
			extract_hmm_dir = 'hmm_PGAP.HMM/'
			os.mkdir(extract_hmm_dir)
			os.system(' '.join(['tar', '-zxf', 'hmm_PGAP.HMM.tgz', '-C', 'hmm_PGAP.HMM/']))
			assert (os.path.isfile(pgap_info_file))
			assert (os.path.isdir(download_path + 'hmm_PGAP.HMM/'))
			for folder, subs, files in os.walk(extract_hmm_dir):
				for filename in files:
					if filename.endswith('.HMM') or filename.endswith('.hmm'):
						hmm_file_path = os.path.abspath(folder) + '/' + filename
						os.system(' '.join(['cat', hmm_file_path, '>>', pgap_phmm_file]))
			pgap_descriptions_file = download_path + 'pgap_descriptions.txt'
			pdf_handle = open(pgap_descriptions_file, 'w')
			with open(pgap_info_file) as opil:
				for i, line in enumerate(opil):
					if i == 0: continue
					line = line.strip()
					ls = line.split('\t')
					label = ls[2]
					description = ls[10]
					pdf_handle.write(label + '\t' + description + '\n')
			pdf_handle.close()
			assert (os.path.isfile(pgap_phmm_file))
			os.system(' '.join(['hmmpress', pgap_phmm_file]))
			z = 0
			with open(pgap_phmm_file) as oppf:
				for line in oppf:
					if line.startswith('NAME'): z += 1

			listing_handle.write('pgap\t' + pgap_descriptions_file + '\t' + pgap_phmm_file + '\t' + str(z) + '\n')
			os.system(' '.join(['rm', '-rf', download_path + 'hmm_PGAP.HMM/', download_path + 'hmm_PGAP.HMM.tgz', pgap_info_file]))
		except Exception as e:
			sys.stderr.write('Issues setting up PGAP database.\n')
			issues_handle.write('Issues setting up PGAP database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

	else:
		# Final annotation files
		pfam_phmm_file = download_path + 'Pfam-A.hmm'
		ko_annot_info_file = download_path + 'ko_list'
		ko_phmm_file = download_path + 'profile.hmm'
		pgap_info_file = download_path + 'hmm_PGAP.tsv'
		pgap_phmm_file = download_path + 'PGAP.hmm'
		pb_faa_file = download_path + 'paperblast.dmnd'
		vog_phmm_file = download_path + 'vog.hmm'
		vog_info_file = download_path + 'vog.annotations.tsv'
		is_faa_file = download_path + 'tn_is.dmnd'
		mb_faa_file = download_path + 'mibig.dmnd'
		card_faa_file = download_path + 'card.dmnd'
		vfdb_faa_file = download_path + 'vfdb.dmnd'
		mobs_faa_file = download_path + 'mobsuite_mpf_and_mob_proteins.dmnd'
		unip_hmm_file = download_path + 'Universal-Hug-et-al.hmm'

		download_links = ['https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz',
						  'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz',
						  'ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz',
						  'http://papers.genomics.lbl.gov/data/uniq.faa',
						  'http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz',
						  'https://dl.secondarymetabolites.org/mibig/mibig_prot_seqs_4.0.fasta',
						  'http://fileshare.csb.univie.ac.at/vog/latest/vog.hmm.tar.gz',
						  'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv',
						  'http://fileshare.csb.univie.ac.at/vog/latest/vog.annotations.tsv.gz',
						  'ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz',
						  'https://tncentral.ncc.unesp.br/api/download_blast/prot/tn_is',
						  'https://card.mcmaster.ca/latest/data',
						  'https://zenodo.org/records/10304948/files/data.tar.gz?download=1',
                    	  'https://zenodo.org/records/13858489/files/Universal-Hug-et-al.hmm?download=1'
						  ]

		# Download
		print('Starting download of files!')
		os.chdir(download_path)
		try:
			for dl in download_links:
				bf = dl.split('/')[-1].split('?')[0]
				download_output = download_path + bf
				if bf == 'data':
					bf = 'card-data.tar.bz2'
					download_output = download_path + bf
				curl_download_dbs_cmd = ['curl', '-L', dl, '-o', download_output]
				os.system(' '.join(curl_download_dbs_cmd))
		except Exception as e:
			sys.stderr.write('Error occurred during downloading!\n')
			issues_handle.write('Error occurred during downloading with curl.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up Pfam database!')
			os.system(' '.join(['gunzip', 'Pfam-A.hmm.gz']))
			assert(os.path.isfile(pfam_phmm_file))
			name = None
			desc = None
			pfam_descriptions_file = download_path + 'pfam_descriptions.txt'
			pdf_handle = open(pfam_descriptions_file, 'w')
			with open(pfam_phmm_file) as oppf:
				for line in oppf:
					line = line.strip()
					ls = line.split()
					if ls[0].strip() == 'NAME':
						name = ' '.join(ls[1:]).strip()
					elif ls[0].strip() == 'DESC':
						desc = ' '.join(ls[1:]).strip()
						pdf_handle.write(name + '\t' + desc + '\n') # type: ignore
			pdf_handle.close()
			os.system(' '.join(['hmmpress', pfam_phmm_file]))
			z = 0
			with open(pfam_phmm_file) as oppf:
				for line in oppf:
					if line.startswith('NAME'): z += 1
			listing_handle.write('pfam\t' + pfam_descriptions_file + '\t' + pfam_phmm_file + '\t' + str(z) + '\n')
		except Exception as e:
			sys.stderr.write('Issues setting up Pfam database.\n')
			issues_handle.write('Issues setting up Pfam database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up KO database!')
			os.system(' '.join(['gunzip', 'ko_list.gz']))
			extract_hmm_dir = 'profiles/'
			os.mkdir(extract_hmm_dir)
			os.system(' '.join(['tar', '-zxf', 'profiles.tar.gz', '-C', extract_hmm_dir]))
			assert(os.path.isfile(ko_annot_info_file))
			assert(os.path.isdir(download_path + 'profiles/'))
			for folder, subs, files in os.walk(extract_hmm_dir):
				for filename in files:
					if filename.endswith('.HMM') or filename.endswith('.hmm'):
						hmm_file_path = os.path.abspath(folder) + '/' + filename
						os.system(' '.join(['cat', hmm_file_path, '>>', ko_phmm_file]))

			assert(os.path.isfile(ko_phmm_file))
			os.system(' '.join(['hmmpress', ko_phmm_file]))

			ko_descriptions_file = download_path + 'ko_descriptions.txt'
			kdf_handle = open(ko_descriptions_file, 'w')
			with open(ko_annot_info_file) as oaif:
				for i, line in enumerate(oaif):
					if i == 0: continue
					line = line.strip()
					ls = line.split('\t')
					ko = ls[0]
					if ls[1] == '-': continue
					description = ls[-1]
					kdf_handle.write(ko + '\t' + description + '\n')
			kdf_handle.close()
			z = 0
			with open(ko_phmm_file) as oppf:
				for line in oppf:
					if line.startswith('NAME'): z += 1
			listing_handle.write('ko\t' + ko_descriptions_file + '\t' + ko_phmm_file + '\t' + str(z) + '\n')
			os.system(' '.join(['rm', '-rf', download_path + 'profiles/', download_path + 'profiles.tar.gz', ko_annot_info_file]))
		except Exception as e:
			sys.stderr.write('Issues setting up KO database.\n')
			issues_handle.write('Issues setting up KO database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up PGAP database!')
			extract_hmm_dir = 'hmm_PGAP.HMM/'
			os.mkdir(extract_hmm_dir)
			os.system(' '.join(['tar', '-zxf', 'hmm_PGAP.HMM.tgz', '-C', 'hmm_PGAP.HMM/']))
			assert (os.path.isfile(pgap_info_file))
			assert (os.path.isdir(download_path + 'hmm_PGAP.HMM/'))
			for folder, subs, files in os.walk(extract_hmm_dir):
				for filename in files:
					if filename.endswith('.HMM') or filename.endswith('.hmm'):
						hmm_file_path = os.path.abspath(folder) + '/' + filename
						os.system(' '.join(['cat', hmm_file_path, '>>', pgap_phmm_file]))
			pgap_descriptions_file = download_path + 'pgap_descriptions.txt'
			pdf_handle = open(pgap_descriptions_file, 'w')
			with open(pgap_info_file) as opil:
				for i, line in enumerate(opil):
					if i == 0: continue
					line = line.strip()
					ls = line.split('\t')
					label = ls[2]
					description = ls[10]
					pdf_handle.write(label + '\t' + description + '\n')
			pdf_handle.close()
			assert(os.path.isfile(pgap_phmm_file))
			os.system(' '.join(['hmmpress', pgap_phmm_file]))
			z = 0
			with open(pgap_phmm_file) as oppf:
				for line in oppf:
					if line.startswith('NAME'): z += 1
			listing_handle.write('pgap\t' + pgap_descriptions_file + '\t' + pgap_phmm_file + '\t' + str(z) + '\n')
			os.system(' '.join(['rm', '-rf', download_path + 'hmm_PGAP.HMM/', download_path + 'hmm_PGAP.HMM.tgz', pgap_info_file]))
		except Exception as e:
			sys.stderr.write('Issues setting up PGAP database.\n')
			issues_handle.write('Issues setting up PGAP database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up MIBiG database!')
			os.system(' '.join(['diamond', 'makedb', '--in', 'mibig_prot_seqs_4.0.fasta', '-d', mb_faa_file, '--threads', str(threads)]))
			assert(os.path.isfile(mb_faa_file))
			mibig_descriptions_file = download_path + 'mibig_descriptions.txt'
			mdf_handle = open(mibig_descriptions_file, 'w')
			with open(download_path + 'mibig_prot_seqs_4.0.fasta') as omf:
				for rec in SeqIO.parse(omf, 'fasta'):
					mdf_handle.write(rec.id + '\t' + rec.description + '\n')
			mdf_handle.close()
			listing_handle.write('mibig\t' + mibig_descriptions_file + '\t' + mb_faa_file + '\tNA\n')
			os.system(' '.join(['rm', '-rf', download_path + 'mibig_prot_seqs_4.0.fasta']))
		except Exception as e:
			sys.stderr.write('Issues setting up MIBiG database.\n')
			issues_handle.write('Issues setting up MIBiG database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up VFDB database!')
			os.system(' '.join(['gunzip', 'VFDB_setB_pro.fas.gz']))
			os.system(' '.join(['diamond', 'makedb', '--in', 'VFDB_setB_pro.fas', '-d', vfdb_faa_file, '--threads', str(threads)]))
			assert(os.path.isfile(mb_faa_file))
			vfdb_descriptions_file = download_path + 'vfdb_descriptions.txt'
			vdf_handle = open(vfdb_descriptions_file, 'w')
			with open(download_path + 'VFDB_setB_pro.fas') as ovf:
				for rec in SeqIO.parse(ovf, 'fasta'):
					vdf_handle.write(rec.id + '\t' + rec.description + '\n')
			vdf_handle.close()
			listing_handle.write('vfdb\t' + vfdb_descriptions_file + '\t' + vfdb_faa_file + '\tNA\n')
			os.system(' '.join(['rm', '-rf', download_path + 'VFDB_setB_pro.fas']))
		except Exception as e:
			sys.stderr.write('Issues setting up VFDB database.\n')
			issues_handle.write('Issues setting up VFDB database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')
		
		try:
			print('Setting up CARD database!')
			os.mkdir(download_path + 'CARD_DB_Files/')
			os.system(' '.join(['tar', '-xf', download_path + 'card-data.tar.bz2', '-C', download_path + 'CARD_DB_Files/']))
			os.system(' '.join(['mv', download_path + 'CARD_DB_Files/protein_fasta_protein_homolog_model.fasta', download_path]))
			os.system(' '.join(['diamond', 'makedb', '--in', download_path + 'protein_fasta_protein_homolog_model.fasta', '-d', card_faa_file, '--threads', str(threads)]))
			card_descriptions_file = download_path + 'card_descriptions.txt'
			cdf_handle = open(card_descriptions_file, 'w')
			with open(download_path + 'protein_fasta_protein_homolog_model.fasta') as ocf:
				for rec in SeqIO.parse(ocf, 'fasta'):
					cdf_handle.write(rec.id + '\t' + rec.description + '\n')
			cdf_handle.close()
			os.system(' '.join(['rm', '-rf', download_path + 'CARD_DB_Files/', download_path + 'protein_fasta_protein_homolog_model.fasta', download_path + 'card-data.tar.bz2']))
			assert(os.path.isfile(card_faa_file))
			listing_handle.write('card\t' + card_descriptions_file + '\t' + card_faa_file + '\tNA\n')
		except Exception as e:
			sys.stderr.write('Issues setting up CARD database.\n')
			issues_handle.write('Issues setting up CARD database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up VOG database!')
			os.mkdir(download_path + 'VOG_HMM_Files/')
			os.system(' '.join(['tar', '-xzf', 'vog.hmm.tar.gz', '-C', download_path + 'VOG_HMM_Files/']))
			os.system(' '.join(['gunzip', download_path + 'vog.annotations.tsv.gz']))
			vog_db_dir = download_path + 'VOG_HMM_Files/'
			for folder, subs, files in os.walk(vog_db_dir):
				for filename in files:
					if filename.endswith('.HMM') or filename.endswith('.hmm'):
						hmm_file_path = os.path.abspath(folder) + '/' + filename
						os.system(' '.join(['cat', hmm_file_path, '>>', vog_phmm_file]))
			vog_descriptions_file = download_path + 'vog_descriptions.txt'
			vdf_handle = open(vog_descriptions_file, 'w')
			with open(download_path + 'vog.annotations.tsv') as ovf:
				for i, line in enumerate(ovf):
					line = line.rstrip('\n')
					ls = line.split('\t')
					if i == 0: continue
					vdf_handle.write(ls[0] + '\t' + ls[4] + '[' + ls[3] + ']\n')
			vdf_handle.close()
			os.system(' '.join(['rm', '-rf', download_path + 'VOG_HMM_Files/', vog_info_file]))
			assert(os.path.isfile(vog_phmm_file))
			assert(os.path.isfile(vog_descriptions_file))
			os.system(' '.join(['hmmpress', vog_phmm_file]))
			z = 0
			with open(vog_phmm_file) as oppf:
				for line in oppf:
					if line.startswith('NAME'): z += 1
			listing_handle.write('vog\t' + vog_descriptions_file + '\t' + vog_phmm_file + '\t' + str(z) + '\n')
		except Exception as e:
			sys.stderr.write('Issues setting up VOG database.\n')
			issues_handle.write('Issues setting up VOG database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up PaperBlast database!')
			pb_faa_path = download_path + 'uniq.faa'
			os.system(' '.join(['diamond', 'makedb', '--in', pb_faa_path, '-d', pb_faa_file, '--threads', str(threads)]))
			paperblast_descriptions_file = download_path + 'paperblast_descriptions.txt'
			pbdf_handle = open(paperblast_descriptions_file, 'w')
			with open(download_path + 'uniq.faa') as opdf:
				for rec in SeqIO.parse(opdf, 'fasta'):
					pbdf_handle.write(rec.id + '\t' + rec.description + '\n')
			pbdf_handle.close()
			os.system(' '.join(['rm', '-rf', pb_faa_path]))
			assert(os.path.isfile(paperblast_descriptions_file))
			assert(os.path.isfile(pb_faa_file))
			listing_handle.write('paperblast\t' + paperblast_descriptions_file + '\t' + pb_faa_file + '\tNA\n')
		except Exception as e:
			sys.stderr.write('Issues setting up PaperBlast database.\n')
			issues_handle.write('Issues setting up PaperBlast database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up TnCentral database!')
			is_zip_path = download_path + 'tncentral_isfinder.prot.fa.zip'
			is_faa_path = download_path + 'tncentral_isfinder.prot.fa'
			is_faa_proc_path = download_path + 'tncentral_isfinder.processed.prot.faa'
			is_descriptions_file = download_path + 'tn_is_descriptions.txt'

			# Extract the .fa file from the zip
			os.system(' '.join(['unzip', '-j', is_zip_path, '-d', download_path]))
			assert(os.path.isfile(is_faa_path))

			idf_handle = open(is_descriptions_file, 'w')
			ifpf_handle = open(is_faa_proc_path, 'w')
			with open(is_faa_path) as oif:
				for rec in SeqIO.parse(oif, 'fasta'):
					ifpf_handle.write('>' + rec.id + '\n' + str(rec.seq).rstrip('*').replace('?', 'X').replace(' ', '') + '\n')
					idf_handle.write(rec.id + '\t' + rec.description + '\n')
			idf_handle.close()
			ifpf_handle.close()

			os.system(' '.join(['diamond', 'makedb', '--in', is_faa_proc_path, '-d', is_faa_file, '--threads', str(threads)]))

			assert(os.path.isfile(is_faa_file))
			assert(os.path.isfile(is_descriptions_file))
			listing_handle.write('tn_is\t' + is_descriptions_file + '\t' + is_faa_file + '\tNA\n')
			os.system(' '.join(['rm', '-rf', is_faa_path, is_faa_proc_path, is_zip_path]))
		except Exception as e:
			sys.stderr.write('Issues setting up TnCentral database.\n')
			issues_handle.write('Issues setting up TnCentral database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up MOB-suite plasmid associated proteins database!')
			ms_tar_path = download_path + 'data.tar.gz'
			ms_faa_path = download_path + 'MOB-suite_proteins.faa'
			ms_descriptions_file = download_path + 'MOB-suite_proteins.txt'
			os.mkdir(download_path + 'MOB-suite_files/')
			os.system(' '.join(['tar', '-xzf', ms_tar_path, '-C', download_path + 'MOB-suite_files/']))
			os.system(' '.join(['cat', download_path + 'MOB-suite_files/data/*.faa', '>', ms_faa_path]))
			os.system(' '.join(['diamond', 'makedb', '--in', ms_faa_path, '-d', mobs_faa_file, '--threads', str(threads)]))
			mdf_handle = open(ms_descriptions_file, 'w')
			with open(ms_faa_path) as oif:
				for rec in SeqIO.parse(oif, 'fasta'):
					mdf_handle.write(rec.id + '\t' + rec.description + '\n')
			mdf_handle.close()
			assert(os.path.isfile(mobs_faa_file))
			assert(os.path.isfile(ms_descriptions_file))
			listing_handle.write('mobsuite\t' + ms_descriptions_file + '\t' + mobs_faa_file + '\tNA\n')
			os.system(' '.join(['rm', '-rf', ms_faa_path, ms_tar_path, download_path + 'MOB-suite_files/']))
		except Exception as e:
			sys.stderr.write('Issues setting up MOB-suite proteins database.\n')
			issues_handle.write('Issues setting MOB-suite proteins database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up Hug et al. 2016 universal ribosomal proteins database!')
			assert(os.path.isfile(unip_hmm_file))

			os.system(' '.join(['hmmpress', unip_hmm_file]))
			z = 0
			with open(unip_hmm_file) as oppf:
				for line in oppf:
					if line.startswith('NAME'): z += 1
			listing_handle.write('riboprots\tNA\t' + unip_hmm_file + '\t' + str(z) + '\n')
		except Exception as e:
			sys.stderr.write('Issues setting up Hug et al. 2016 database.\n')
			issues_handle.write('Issues setting up Hug et al. 2016 database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

	listing_handle.close()
	issues_handle.close()

	print("Done setting up annotation databases!")
	print(f"Information on final files used by zol can be found at:\n{listing_file}")
	sys.exit(0)
setup_annot_dbs()
