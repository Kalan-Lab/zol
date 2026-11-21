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
import datetime
from pathlib import Path
import glob
import re
import concurrent.futures
import email.utils
import threading
import time
import shlex
import subprocess

import codecs


def _latin1_fallback(err):
	"""
	Codec error handler: decode bytes that are not valid UTF-8 as Latin-1 instead of failing. Some databases
	(e.g. VFDB) have stray Latin-1 / Windows-1252 bytes in FASTA headers - non-breaking spaces (0xa0) become
	regular spaces.
	"""
	bad_bytes = err.object[err.start:err.end]
	return (''.join(' ' if b == 0xa0 else bytes([b]).decode('latin-1') for b in bad_bytes), err.end)


codecs.register_error('zol_latin1_fallback', _latin1_fallback)


def _open_fasta_text(fasta_file):
	"""Open a downloaded FASTA file for parsing, tolerating non-UTF-8 bytes in headers."""
	return open(fasta_file, encoding='utf-8', errors='zol_latin1_fallback')


def _run_cmd(cmd):
	"""Run a command (given as a list of arguments, no shell) and warn if it fails."""
	cmd = [str(x) for x in cmd]
	try:
		result = subprocess.run(cmd)
		if result.returncode == 0:
			return True
	except Exception as e:
		sys.stderr.write(f'Error: {e}\n')
	sys.stderr.write('Warning: issue running: ' + ' '.join(shlex.quote(x) for x in cmd) + '\n')
	return False


def _append_file(source_file, target_file):
	"""Append the contents of source_file to target_file (replacement for `cat source >> target`)."""
	with open(source_file, 'rb') as src, open(target_file, 'ab') as tgt:
		shutil.copyfileobj(src, tgt)


def _concat_files(source_files, target_file):
	"""Concatenate source_files into target_file (replacement for `cat sources > target`)."""
	with open(target_file, 'wb') as tgt:
		for source_file in source_files:
			with open(source_file, 'rb') as src:
				shutil.copyfileobj(src, tgt)


SEGMENT_MIN_SIZE = 64 * 1024 * 1024  # files at least this large are downloaded in parallel byte-range segments
MAX_PARALLEL_FILES = 6
_print_lock = threading.Lock()


def _log(msg):
	with _print_lock:
		print(msg, flush=True)


def _default_download_name(url):
	bf = url.split('/')[-1].split('?')[0]
	if bf == 'data':  # CARD
		bf = 'card-data.tar.bz2'
	return bf


def _probe_remote(url):
	"""Return (size in bytes or None, Last-Modified as epoch seconds or None, whether byte ranges are honored)."""
	size, mtime = None, None
	try:
		res = subprocess.run(['curl', '-s', '-I', '-L', '--max-time', '60', url], capture_output=True, text=True)
		# with redirects, curl prints a header block per hop - the last value seen is the final one
		for line in res.stdout.splitlines():
			if line.upper().startswith('HTTP/'):
				size, mtime = None, None  # new response (redirect hop) - only the final one is relevant
				continue
			key, _, val = line.partition(':')
			key, val = key.strip().lower(), val.strip()
			if key == 'content-length' and val.isdigit():
				size = int(val) if int(val) > 0 else None
			elif key == 'last-modified':
				try:
					mtime = email.utils.parsedate_to_datetime(val).timestamp()
				except Exception:
					pass
	except Exception:
		pass
	supports_ranges = False
	if size is not None and size >= SEGMENT_MIN_SIZE:
		# servers do not always advertise range support (e.g. Zenodo), so test with an actual 1-byte range request
		# (reading at most 2 bytes so that a server ignoring the range does not start a full download)
		try:
			proc = subprocess.Popen(['curl', '-s', '-L', '--max-time', '60', '-r', '0-0', url], stdout=subprocess.PIPE)
			data = proc.stdout.read(2)
			proc.kill()
			proc.wait()
			supports_ranges = (len(data) == 1)
		except Exception:
			supports_ranges = False
	return size, mtime, supports_ranges


def _curl_ok(cmd):
	try:
		return subprocess.run(cmd).returncode == 0
	except Exception as e:
		sys.stderr.write(f'Error: {e}\n')
		return False


def _download_file(url, output_file, segments, progress):
	"""Download url to output_file - in parallel byte-range segments when the file is large and the server allows."""
	size, mtime, supports_ranges = _probe_remote(url)
	if os.path.isfile(output_file):
		os.remove(output_file)
	curl_base = ['curl', '-s', '-S', '-L', '--retry', '5', '--retry-delay', '5', '--connect-timeout', '60']

	if segments > 1 and supports_ranges and size is not None:
		chunk = -(-size // segments)  # ceiling division
		bounds = [(start, min(start + chunk, size) - 1) for start in range(0, size, chunk)]
		parts = [f'{output_file}.part{i}' for i in range(len(bounds))]
		progress[output_file] = (size, parts)

		def fetch(i):
			start, end = bounds[i]
			for _ in range(3):
				_curl_ok(curl_base + ['-r', f'{start}-{end}', '-o', parts[i], url])
				if os.path.isfile(parts[i]) and os.path.getsize(parts[i]) == end - start + 1:
					return True
			return False

		with concurrent.futures.ThreadPoolExecutor(max_workers=len(bounds)) as ex:
			success = all(ex.map(fetch, range(len(bounds))))
		if success:
			_concat_files(parts, output_file)
		for part in parts:
			if os.path.isfile(part):
				os.remove(part)
		success = success and os.path.getsize(output_file) == size
	else:
		progress[output_file] = (size, [output_file])
		success = _curl_ok(curl_base + ['-R', url, '-o', output_file]) and os.path.isfile(output_file)
		if success and size is not None and os.path.getsize(output_file) != size:
			_log(f'Warning: size of {output_file} differs from the size reported by the server for {url}')

	progress.pop(output_file, None)
	if success and mtime is not None:
		os.utime(output_file, (mtime, mtime))  # keep the server's timestamp (like curl -R) for release info
	return success


def _download_all(urls, download_path, segments=8, output_names=None):
	"""
	Download several files concurrently (large ones additionally in parallel segments), printing periodic progress.
	Returns the list of URLs which failed to download.
	"""
	output_names = output_names or {}
	jobs = [(url, download_path + output_names.get(url, _default_download_name(url))) for url in urls]
	progress = {}
	done = threading.Event()

	def report():
		start = time.time()
		while not done.wait(20):
			got, total = 0, 0
			for size, files in list(progress.values()):
				got += sum(os.path.getsize(f) for f in files if os.path.isfile(f))
				total += size or 0
			if progress:
				elapsed = time.time() - start
				_log(f'  ... {len(progress)} download(s) in progress: {got / 1e6:,.0f} of {total / 1e6:,.0f} MB '
					 f'({elapsed / 60:.1f} min elapsed)')

	reporter = threading.Thread(target=report, daemon=True)
	reporter.start()
	failures = []

	def run(job):
		url, output_file = job
		_log(f'Downloading {url}')
		ok = _download_file(url, output_file, segments, progress)
		_log(('Finished ' if ok else 'FAILED to download ') + url)
		return url, ok

	with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_PARALLEL_FILES) as ex:
		for url, ok in ex.map(run, jobs):
			if not ok:
				failures.append(url)
	done.set()
	if failures:
		sys.stderr.write('Warning: the following downloads failed:\n' + '\n'.join(failures) + '\n')
	return failures


def _write_text(text, target_file):
	"""Write a single line of text to target_file (replacement for `echo text > target`)."""
	with open(target_file, 'w') as tgt:
		tgt.write(text + '\n')


DEFENSEFINDER_MODELS_VERSION = '3.1.0'
DEFENSEFINDER_MODELS_URL = f'https://github.com/mdmparis/defense-finder-models/archive/refs/tags/{DEFENSEFINDER_MODELS_VERSION}.tar.gz'
DEFENSEFINDER_TAR_NAME = 'defense-finder-models.tar.gz'


def _setup_defensefinder_db(download_path, download_release_info_path, listing_handle, issues_handle):
	"""
	Download the DefenseFinder models (https://github.com/mdmparis/defense-finder-models), concatenate the
	profile HMMs into a single pressed HMM database, and write a descriptions file mapping each HMM name to
	the gene and the defense system (or anti-defense factor category) it belongs to.
	"""
	try:
		print('Setting up DefenseFinder database!')
		df_tar_file = download_path + DEFENSEFINDER_TAR_NAME
		df_extract_dir = download_path + 'DefenseFinder_Files/'
		df_phmm_file = download_path + 'defensefinder.hmm'
		df_descriptions_file = download_path + 'defensefinder_descriptions.txt'
		if not os.path.isfile(df_tar_file):
			_run_cmd(['curl', '-R', '-L', DEFENSEFINDER_MODELS_URL, '-o', df_tar_file])
		if not (os.path.isfile(df_tar_file)):
			raise AssertionError()
		os.mkdir(df_extract_dir)
		_run_cmd(['tar', '-xzf', df_tar_file, '-C', df_extract_dir])
		model_dirs = glob.glob(df_extract_dir + 'defense-finder-models-*/')
		if len(model_dirs) != 1:
			raise AssertionError()
		model_dir = model_dirs[0]
		profiles = sorted(glob.glob(model_dir + 'profiles/*.hmm'))
		if len(profiles) == 0:
			raise AssertionError()

		# Profiles named System__Gene belong to defense systems; the remaining ones are anti-defense factors
		# (ADFs), whose category (e.g. Anti_CRISPR) is given by the definitions/ADF/<category>/ folder of the
		# MacSyFinder model that references the gene.
		adf_category = {}
		for xml_file in glob.glob(model_dir + 'definitions/ADF/**/*.xml', recursive=True):
			category = os.path.relpath(xml_file, model_dir + 'definitions/ADF/').split(os.sep)[0]
			with open(xml_file) as oxf:
				for gene in re.findall(r'<gene[^>]*\sname="([^"]+)"', oxf.read()):
					adf_category.setdefault(gene, category)

		ddf_handle = open(df_descriptions_file, 'w')
		for profile in profiles:
			profile_id = os.path.basename(profile)[:-len('.hmm')]
			if '__' in profile_id:
				system, gene = profile_id.split('__', 1)
				description = f'{gene} ({system} defense system)'
			else:
				category = adf_category.get(profile_id, 'Other')
				description = f'{profile_id} ({category} anti-defense factor)'
			# key on the NAME field in the HMM itself, which is what hmmsearch reports (not always the file name)
			hmm_name = None
			with open(profile) as opf:
				for line in opf:
					if line.startswith('NAME'):
						hmm_name = line.split()[1].strip()
						break
			if hmm_name is None:
				raise AssertionError()
			ddf_handle.write(hmm_name + '\t' + description + '\n')
			_append_file(profile, df_phmm_file)
		ddf_handle.close()

		if os.path.isfile(model_dir + 'LICENSE'):
			shutil.copy(model_dir + 'LICENSE', download_release_info_path + 'defensefinder_license.txt')
		_write_text(f"v{DEFENSEFINDER_MODELS_VERSION} (https://github.com/mdmparis/defense-finder-models/releases/tag/{DEFENSEFINDER_MODELS_VERSION})",
			  download_release_info_path + 'defensefinder_release_notes.txt')

		_run_cmd(['hmmpress', '-f', df_phmm_file])
		z = 0
		with open(df_phmm_file) as oppf:
			for line in oppf:
				if line.startswith('NAME'): z += 1
		listing_handle.write('defensefinder\t' + df_descriptions_file + '\t' + df_phmm_file + '\t' + str(z) + '\n')
		_run_cmd(['rm', '-rf', df_extract_dir, df_tar_file])
	except Exception as e:
		sys.stderr.write('Issues setting up DefenseFinder database.\n')
		issues_handle.write('Issues setting up DefenseFinder database.\n')
		sys.stderr.write(traceback.format_exc())
		sys.stderr.write(str(e) + '\n')


def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: setup_annotation_dbs.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
		
	Downloads annotation databases for annotations and lateral transfer inference: KOfam, 
	PGAP, PaperBlast, MIBiG, CARD, Conserved Ribosomal Proteins, VOGs, ISFinder, VFDB,
	MOB-suite MPF and MOB Proteins, DefenseFinder, and Pfam.
								  
	Location of where to download databases is controlled by setting the environmental
	variable ZOL_DATA_PATH e.g.:
								  							   
	$ export ZOL_DATA_PATH=/path/to/database_location/zol_db/

	CAUTION: This directory will be deleted and recreated when this script is run - it 
	should therefore be unique and specific to zol.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-c', '--threads', type=int, help="Number of threads to use [Default is 4].", required=False, default=4)
	parser.add_argument('-m', '--minimal', action='store_true', help="Minimal mode - will only download Pfam and PGAP HMMs.", required=False, default=False)
	parser.add_argument('-ld', '--lsabgc-minimal', action='store_true', help="Minimal mode for lsaBGC - will only download Pfam HMMs, PGAP HMMs, & MIBiG proteins.", required=False, default=False)	
	parser.add_argument('-ds', '--download-segments', type=int, help="Number of parallel connections used to download each large file\n(if the server supports it) [Default is 8; use 1 to disable].", required=False, default=8)
	parser.add_argument('-y', '--yes', action='store_true', help="Assume 'yes' for prompts (non-interactive).", required=False, default=False)

	args = parser.parse_args()
	return args

def get_file_timestamp(infile):
	infile = Path(infile)
	
	# Fetch the raw Unix timestamp (float)
	modification_timestamp = infile.stat().st_mtime
	# Convert to a readable datetime object
	readable_time = datetime.datetime.fromtimestamp(modification_timestamp)

	return readable_time


def setup_annot_dbs():
	myargs = create_parser()

	download_path = None
	if str(os.getenv("ZOL_DATA_PATH")) != 'None':
		download_path = str(os.getenv("ZOL_DATA_PATH"))
	if download_path is None or not os.path.isdir(download_path):
		sys.stderr.write('Issues validing database download directory exists.\n')
		sys.exit(1)
	else:
		download_path = os.path.abspath(download_path) + "/"
	
	threads = myargs.threads
	download_segments = max(1, myargs.download_segments)
	minimal_mode = myargs.minimal
	lsabgc_minimal_mode = myargs.lsabgc_minimal
	non_interactive_yes = myargs.yes or str(os.getenv('ZOL_NONINTERACTIVE', '0')).strip() in ['1', 'true', 'True']

	try:
		if not (os.path.isdir(download_path)):
			raise AssertionError()
		if not non_interactive_yes:
			response = input(f"The directory {download_path}\nalready exists, will delete it and recreate it. (This directory\nshould be specific to zol not a general directory for\ndatabases) Proceed with deleting? (yes/no): ")
			if response.lower() != 'yes':
				sys.stderr.write('Deletion not requested! Exiting ...\n')
				sys.exit(1)
	except Exception as e:
		sys.stderr.write('Error: Provided directory for downloading annotation files does not exist or user did not accept deleting the directory and recreating it!\n')

	download_release_info_path = download_path + 'release_info/'

	if lsabgc_minimal_mode:
		sys.stdout.write('lsaBGC minimal mode requested, will only be downloading the Pfam, MIBiG, and PGAP databases.\n')
	elif minimal_mode:
		sys.stdout.write('Minimal mode requested, will only be downloading the Pfam and PGAP databases.\n')
	
	try:
		shutil.rmtree(download_path)
		os.mkdir(download_path)
		os.mkdir(download_release_info_path)
		readme_outf = open(download_path + 'README.txt', 'w')
		readme_outf.write('Default space for downloading annotation databases.\n')
		readme_outf.write('Date of database setup: ' + str(datetime.datetime.now()) + '\n')
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
		if _download_all(download_links, download_path, segments=download_segments):
			raise AssertionError('download of weight files failed')
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

		download_links = ['https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz',
				  'https://dl.secondarymetabolites.org/mibig/mibig_prot_seqs_4.0.fasta',
				  'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv',
				  'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz', # versioning files after this 
				  'https://ftp.ncbi.nlm.nih.gov/hmm/current/RELEASE_NOTES.txt',
				  'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt']

		# Download
		print('Starting download of files!')
		os.chdir(download_path)
		try:
			_download_all(download_links, download_path, segments=download_segments)
		except Exception as e:
			sys.stderr.write('Error occurred during downloading of databases!\n')
			sys.stderr.write(str(e) + '\n')

		# rename record files 
		try:
			if os.path.isfile(download_path + 'RELEASE_NOTES.txt'):
				rename_pgap_release_notes_cmd = ['mv', download_path + 'RELEASE_NOTES.txt', download_release_info_path + 'pgap_release_notes.txt']
				_run_cmd(rename_pgap_release_notes_cmd)

			if os.path.isfile(download_path + 'relnotes.txt'):
				rename_pfam_release_notes_cmd = ['mv', download_path + 'relnotes.txt', download_release_info_path + 'pfam_release_notes.txt']
				_run_cmd(rename_pfam_release_notes_cmd)
				
			if os.path.isfile(download_path + 'mibig_prot_seqs_4.0.fasta'):
				_write_text("v4.0 (November 15, 2024)", download_release_info_path + 'mibig_release_notes.txt')

		except Exception as e:
			sys.stderr.write('Issues with renaming record/provenance files.\n')
			issues_handle.write('Issues with renaming record/provenance files.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up Pfam database!')
			_run_cmd(['gunzip', 'Pfam-A.hmm.gz'])
			if not (os.path.isfile(pfam_phmm_file)):
				raise AssertionError()
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
			_run_cmd(['hmmpress', '-f', pfam_phmm_file])
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
			_run_cmd(['tar', '-zxf', 'hmm_PGAP.HMM.tgz', '-C', 'hmm_PGAP.HMM/'])
			if not (os.path.isfile(pgap_info_file)):
				raise AssertionError()
			if not (os.path.isdir(download_path + 'hmm_PGAP.HMM/')):
				raise AssertionError()
			for folder, subs, files in os.walk(extract_hmm_dir):
				for filename in files:
					if filename.endswith('.HMM') or filename.endswith('.hmm'):
						hmm_file_path = os.path.abspath(folder) + '/' + filename
						_append_file(hmm_file_path, pgap_phmm_file)
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
			if not (os.path.isfile(pgap_phmm_file)):
				raise AssertionError()
			_run_cmd(['hmmpress', '-f', pgap_phmm_file])
			z = 0
			with open(pgap_phmm_file) as oppf:
				for line in oppf:
					if line.startswith('NAME'): z += 1

			listing_handle.write('pgap\t' + pgap_descriptions_file + '\t' + pgap_phmm_file + '\t' + str(z) + '\n')
			_run_cmd(['rm', '-rf', download_path + 'hmm_PGAP.HMM/', download_path + 'hmm_PGAP.HMM.tgz', pgap_info_file])
		except Exception as e:
			sys.stderr.write('Issues setting up PGAP database.\n')
			issues_handle.write('Issues setting up PGAP database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up MIBiG database!')
			_run_cmd(['diamond', 'makedb', '--in', 'mibig_prot_seqs_4.0.fasta', '-d', mb_faa_file, '--threads', str(threads)])
			if not (os.path.isfile(mb_faa_file)):
				raise AssertionError()
			mibig_descriptions_file = download_path + 'mibig_descriptions.txt'
			mdf_handle = open(mibig_descriptions_file, 'w')
			with _open_fasta_text(download_path + 'mibig_prot_seqs_4.0.fasta') as omf:
				for rec in SeqIO.parse(omf, 'fasta'):
					mdf_handle.write(rec.id + '\t' + rec.description + '\n')
			mdf_handle.close()
			listing_handle.write('mibig\t' + mibig_descriptions_file + '\t' + mb_faa_file + '\tNA\n')
			_run_cmd(['rm', '-rf', download_path + 'mibig_prot_seqs_4.0.fasta'])
		except Exception as e:
			sys.stderr.write('Issues setting up MIBiG database.\n')
			issues_handle.write('Issues setting up MIBiG database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')
			


	elif minimal_mode:
		# Final annotation files
		pgap_info_file = download_path + 'hmm_PGAP.tsv'
		pgap_phmm_file = download_path + 'PGAP.hmm'
		pfam_phmm_file = download_path + 'Pfam-A.hmm'

		download_links = ['https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz',
						  'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv',
						  'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz', # versioning files after this 
						  'https://ftp.ncbi.nlm.nih.gov/hmm/current/RELEASE_NOTES.txt',
     			  		  'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt']

		# Download
		print('Starting download of files!')
		os.chdir(download_path)
		try:
			_download_all(download_links, download_path, segments=download_segments)
		except Exception as e:
			sys.stderr.write('Error occurred during downloading of databases!\n')
			sys.stderr.write(str(e) + '\n')

		# rename record files 
		try:
			if os.path.isfile(download_path + 'RELEASE_NOTES.txt'):
				rename_pgap_release_notes_cmd = ['mv', download_path + 'RELEASE_NOTES.txt', download_release_info_path + 'pgap_release_notes.txt']
				_run_cmd(rename_pgap_release_notes_cmd)

			if os.path.isfile(download_path + 'relnotes.txt'):
				rename_pfam_release_notes_cmd = ['mv', download_path + 'relnotes.txt', download_release_info_path + 'pfam_release_notes.txt']
				_run_cmd(rename_pfam_release_notes_cmd)

		except Exception as e:
			sys.stderr.write('Issues with renaming record/provenance files.\n')
			issues_handle.write('Issues with renaming record/provenance files.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up Pfam database!')
			_run_cmd(['gunzip', 'Pfam-A.hmm.gz'])
			if not (os.path.isfile(pfam_phmm_file)):
				raise AssertionError()
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
			_run_cmd(['hmmpress', '-f', pfam_phmm_file])
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
			_run_cmd(['tar', '-zxf', 'hmm_PGAP.HMM.tgz', '-C', 'hmm_PGAP.HMM/'])
			if not (os.path.isfile(pgap_info_file)):
				raise AssertionError()
			if not (os.path.isdir(download_path + 'hmm_PGAP.HMM/')):
				raise AssertionError()
			for folder, subs, files in os.walk(extract_hmm_dir):
				for filename in files:
					if filename.endswith('.HMM') or filename.endswith('.hmm'):
						hmm_file_path = os.path.abspath(folder) + '/' + filename
						_append_file(hmm_file_path, pgap_phmm_file)
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
			if not (os.path.isfile(pgap_phmm_file)):
				raise AssertionError()
			_run_cmd(['hmmpress', '-f', pgap_phmm_file])
			z = 0
			with open(pgap_phmm_file) as oppf:
				for line in oppf:
					if line.startswith('NAME'): z += 1

			listing_handle.write('pgap\t' + pgap_descriptions_file + '\t' + pgap_phmm_file + '\t' + str(z) + '\n')
			_run_cmd(['rm', '-rf', download_path + 'hmm_PGAP.HMM/', download_path + 'hmm_PGAP.HMM.tgz', pgap_info_file])
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
		is_faa_file = download_path + 'isfinder.dmnd'
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
						  'https://github.com/raufs/ISfinder-sequences/raw/refs/heads/master/IS.faa',
						  'https://card.mcmaster.ca/latest/data',
						  'https://zenodo.org/records/10304948/files/data.tar.gz?download=1',
                    	  'https://zenodo.org/records/13858489/files/Universal-Hug-et-al.hmm?download=1', # versioning files after this 
						  'https://ftp.ncbi.nlm.nih.gov/hmm/current/RELEASE_NOTES.txt',
     			  		  'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt', 
						  'https://fileshare.lisc.univie.ac.at//vog/latest/release.txt']

		# Download
		print('Starting download of files!')
		os.chdir(download_path)
		try:
			download_links.append(DEFENSEFINDER_MODELS_URL)
			failed_downloads = _download_all(download_links, download_path, segments=download_segments,
								output_names={DEFENSEFINDER_MODELS_URL: DEFENSEFINDER_TAR_NAME})
			for dl in failed_downloads:
				issues_handle.write(f'Failed to download: {dl}\n')
		except Exception as e:
			sys.stderr.write('Error occurred during downloading!\n')
			issues_handle.write('Error occurred during downloading with curl.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		# rename record files 
		try:
			if os.path.isfile(download_path + 'RELEASE_NOTES.txt'):
				rename_pgap_release_notes_cmd = ['mv', download_path + 'RELEASE_NOTES.txt', download_release_info_path + 'pgap_release_notes.txt']
				_run_cmd(rename_pgap_release_notes_cmd)

			if os.path.isfile(download_path + 'relnotes.txt'):
				rename_pfam_release_notes_cmd = ['mv', download_path + 'relnotes.txt', download_release_info_path + 'pfam_release_notes.txt']
				_run_cmd(rename_pfam_release_notes_cmd)
				
			if os.path.isfile(download_path + 'mibig_prot_seqs_4.0.fasta'):
				_write_text(f"v4.0 (November 15, 2024)", download_release_info_path + "mibig_release_notes.txt")

			if os.path.isfile(download_path + 'profiles.tar.gz'):
				timestamp = get_file_timestamp(download_path + 'profiles.tar.gz')
				_write_text(f"Timestamp: {str(timestamp)}", download_release_info_path + "kofam_release_notes.txt")
			
			if os.path.isfile(download_path + 'uniq.faa'):
				timestamp = get_file_timestamp(download_path + 'uniq.faa')
				_write_text(f"Timestamp: {str(timestamp)}", download_release_info_path + "paperblast_release_notes.txt")
			
			if os.path.isfile(download_path + 'VFDB_setB_pro.fas.gz'):
				timestamp = get_file_timestamp(download_path + 'VFDB_setB_pro.fas.gz')
				_write_text(f"Timestamp: {str(timestamp)}", download_release_info_path + "vfdb_release_notes.txt")
			
			if os.path.isfile(download_path + 'release.txt'):
				rename_vog_release_notes_cmd = ['mv', download_path + 'release.txt', download_release_info_path + 'vog_release_notes.txt']
				_run_cmd(rename_vog_release_notes_cmd)
			
			if os.path.isfile(download_path + 'card-data.tar.bz2'):
				timestamp = get_file_timestamp(download_path + 'card-data.tar.bz2')
				_write_text(f"Timestamp: {str(timestamp)}", download_release_info_path + "card_release_notes.txt")
			
			if os.path.isfile(download_path + 'data.tar.gz'):
				_write_text(f"v3.1.8 (December 8, 2023)", download_release_info_path + "mobsuite_release_notes.txt")

		except Exception as e:
			sys.stderr.write('Issues with renaming record/provenance files.\n')
			issues_handle.write('Issues with renaming record/provenance files.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up Pfam database!')
			_run_cmd(['gunzip', 'Pfam-A.hmm.gz'])
			if not (os.path.isfile(pfam_phmm_file)):
				raise AssertionError()
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
			_run_cmd(['hmmpress', '-f', pfam_phmm_file])
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
			_run_cmd(['gunzip', 'ko_list.gz'])
			extract_hmm_dir = 'profiles/'
			os.mkdir(extract_hmm_dir)
			_run_cmd(['tar', '-zxf', 'profiles.tar.gz', '-C', extract_hmm_dir])
			if not (os.path.isfile(ko_annot_info_file)):
				raise AssertionError()
			if not (os.path.isdir(download_path + 'profiles/')):
				raise AssertionError()
			for folder, subs, files in os.walk(extract_hmm_dir):
				for filename in files:
					if filename.endswith('.HMM') or filename.endswith('.hmm'):
						hmm_file_path = os.path.abspath(folder) + '/' + filename
						_append_file(hmm_file_path, ko_phmm_file)

			if not (os.path.isfile(ko_phmm_file)):
				raise AssertionError()
			_run_cmd(['hmmpress', '-f', ko_phmm_file])

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
			_run_cmd(['rm', '-rf', download_path + 'profiles/', download_path + 'profiles.tar.gz', ko_annot_info_file])
		except Exception as e:
			sys.stderr.write('Issues setting up KO database.\n')
			issues_handle.write('Issues setting up KO database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up PGAP database!')
			extract_hmm_dir = 'hmm_PGAP.HMM/'
			os.mkdir(extract_hmm_dir)
			_run_cmd(['tar', '-zxf', 'hmm_PGAP.HMM.tgz', '-C', 'hmm_PGAP.HMM/'])
			if not (os.path.isfile(pgap_info_file)):
				raise AssertionError()
			if not (os.path.isdir(download_path + 'hmm_PGAP.HMM/')):
				raise AssertionError()
			for folder, subs, files in os.walk(extract_hmm_dir):
				for filename in files:
					if filename.endswith('.HMM') or filename.endswith('.hmm'):
						hmm_file_path = os.path.abspath(folder) + '/' + filename
						_append_file(hmm_file_path, pgap_phmm_file)
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
			if not (os.path.isfile(pgap_phmm_file)):
				raise AssertionError()
			_run_cmd(['hmmpress', '-f', pgap_phmm_file])
			z = 0
			with open(pgap_phmm_file) as oppf:
				for line in oppf:
					if line.startswith('NAME'): z += 1
			listing_handle.write('pgap\t' + pgap_descriptions_file + '\t' + pgap_phmm_file + '\t' + str(z) + '\n')
			_run_cmd(['rm', '-rf', download_path + 'hmm_PGAP.HMM/', download_path + 'hmm_PGAP.HMM.tgz', pgap_info_file])
		except Exception as e:
			sys.stderr.write('Issues setting up PGAP database.\n')
			issues_handle.write('Issues setting up PGAP database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up MIBiG database!')
			_run_cmd(['diamond', 'makedb', '--in', 'mibig_prot_seqs_4.0.fasta', '-d', mb_faa_file, '--threads', str(threads)])
			if not (os.path.isfile(mb_faa_file)):
				raise AssertionError()
			mibig_descriptions_file = download_path + 'mibig_descriptions.txt'
			mdf_handle = open(mibig_descriptions_file, 'w')
			with _open_fasta_text(download_path + 'mibig_prot_seqs_4.0.fasta') as omf:
				for rec in SeqIO.parse(omf, 'fasta'):
					mdf_handle.write(rec.id + '\t' + rec.description + '\n')
			mdf_handle.close()
			listing_handle.write('mibig\t' + mibig_descriptions_file + '\t' + mb_faa_file + '\tNA\n')
			_run_cmd(['rm', '-rf', download_path + 'mibig_prot_seqs_4.0.fasta'])
		except Exception as e:
			sys.stderr.write('Issues setting up MIBiG database.\n')
			issues_handle.write('Issues setting up MIBiG database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up VFDB database!')
			_run_cmd(['gunzip', 'VFDB_setB_pro.fas.gz'])
			_run_cmd(['diamond', 'makedb', '--in', 'VFDB_setB_pro.fas', '-d', vfdb_faa_file, '--threads', str(threads)])
			if not (os.path.isfile(vfdb_faa_file)):
				raise AssertionError()
			vfdb_descriptions_file = download_path + 'vfdb_descriptions.txt'
			vdf_handle = open(vfdb_descriptions_file, 'w')
			with _open_fasta_text(download_path + 'VFDB_setB_pro.fas') as ovf:
				for rec in SeqIO.parse(ovf, 'fasta'):
					vdf_handle.write(rec.id + '\t' + rec.description + '\n')
			vdf_handle.close()
			listing_handle.write('vfdb\t' + vfdb_descriptions_file + '\t' + vfdb_faa_file + '\tNA\n')
			_run_cmd(['rm', '-rf', download_path + 'VFDB_setB_pro.fas'])
		except Exception as e:
			sys.stderr.write('Issues setting up VFDB database.\n')
			issues_handle.write('Issues setting up VFDB database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up CARD database!')
			print('Check out the ARTS webserver for more detailed analysis in finding BGCs for synthesizing antibiotics.')
			os.mkdir(download_path + 'CARD_DB_Files/')
			_run_cmd(['tar', '-xf', download_path + 'card-data.tar.bz2', '-C', download_path + 'CARD_DB_Files/'])
			_run_cmd(['mv', download_path + 'CARD_DB_Files/protein_fasta_protein_homolog_model.fasta', download_path])
			_run_cmd(['mv', download_path + 'CARD_DB_Files/CARD-Download-README.txt', download_release_info_path + '/card_license.txt'])

			_run_cmd(['diamond', 'makedb', '--in', download_path + 'protein_fasta_protein_homolog_model.fasta', '-d', card_faa_file, '--threads', str(threads)])
			card_descriptions_file = download_path + 'card_descriptions.txt'
			cdf_handle = open(card_descriptions_file, 'w')
			with _open_fasta_text(download_path + 'protein_fasta_protein_homolog_model.fasta') as ocf:
				for rec in SeqIO.parse(ocf, 'fasta'):
					cdf_handle.write(rec.id + '\t' + rec.description + '\n')
			cdf_handle.close()
			_run_cmd(['rm', '-rf', download_path + 'CARD_DB_Files/', download_path + 'protein_fasta_protein_homolog_model.fasta', download_path + 'card-data.tar.bz2'])
			if not (os.path.isfile(card_faa_file)):
				raise AssertionError()
			listing_handle.write('card\t' + card_descriptions_file + '\t' + card_faa_file + '\tNA\n')
		except Exception as e:
			sys.stderr.write('Issues setting up CARD database.\n')
			issues_handle.write('Issues setting up CARD database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up VOG database!')
			os.mkdir(download_path + 'VOG_HMM_Files/')
			_run_cmd(['tar', '-xzf', 'vog.hmm.tar.gz', '-C', download_path + 'VOG_HMM_Files/'])
			_run_cmd(['gunzip', download_path + 'vog.annotations.tsv.gz'])
			vog_db_dir = download_path + 'VOG_HMM_Files/'
			for folder, subs, files in os.walk(vog_db_dir):
				for filename in files:
					if filename.endswith('.HMM') or filename.endswith('.hmm'):
						hmm_file_path = os.path.abspath(folder) + '/' + filename
						_append_file(hmm_file_path, vog_phmm_file)
			vog_descriptions_file = download_path + 'vog_descriptions.txt'
			vdf_handle = open(vog_descriptions_file, 'w')
			with open(download_path + 'vog.annotations.tsv') as ovf:
				for i, line in enumerate(ovf):
					line = line.rstrip('\n')
					ls = line.split('\t')
					if i == 0: continue
					vdf_handle.write(ls[0] + '\t' + ls[4] + '[' + ls[3] + ']\n')
			vdf_handle.close()
			_run_cmd(['rm', '-rf', download_path + 'VOG_HMM_Files/', vog_info_file])
			if not (os.path.isfile(vog_phmm_file)):
				raise AssertionError()
			if not (os.path.isfile(vog_descriptions_file)):
				raise AssertionError()
			_run_cmd(['hmmpress', '-f', vog_phmm_file])
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
			_run_cmd(['diamond', 'makedb', '--in', pb_faa_path, '-d', pb_faa_file, '--threads', str(threads)])
			paperblast_descriptions_file = download_path + 'paperblast_descriptions.txt'
			pbdf_handle = open(paperblast_descriptions_file, 'w')
			with _open_fasta_text(download_path + 'uniq.faa') as opdf:
				for rec in SeqIO.parse(opdf, 'fasta'):
					pbdf_handle.write(rec.id + '\t' + rec.description + '\n')
			pbdf_handle.close()
			_run_cmd(['rm', '-rf', pb_faa_path])
			if not (os.path.isfile(paperblast_descriptions_file)):
				raise AssertionError()
			if not (os.path.isfile(pb_faa_file)):
				raise AssertionError()
			listing_handle.write('paperblast\t' + paperblast_descriptions_file + '\t' + pb_faa_file + '\tNA\n')
		except Exception as e:
			sys.stderr.write('Issues setting up PaperBlast database.\n')
			issues_handle.write('Issues setting up PaperBlast database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up ISFinder database!')
			is_faa_path = download_path + 'IS.faa'
			is_faa_proc_path = download_path + 'IS.processed.faa'
			is_descriptions_file = download_path + 'is_descriptions.txt'

			idf_handle = open(is_descriptions_file, 'w')
			ifpf_handle = open(is_faa_proc_path, 'w')
			with _open_fasta_text(is_faa_path) as oif:
				for rec in SeqIO.parse(oif, 'fasta'):
					if '~~~Passenger Gene~~~' in rec.description: continue
					ifpf_handle.write('>' + rec.description + '\n' + str(rec.seq) + '\n')
					idf_handle.write(rec.id + '\t' + rec.description + '\n')
			ifpf_handle.close()
			idf_handle.close()

			_run_cmd(['diamond', 'makedb', '--in', is_faa_proc_path, '-d', is_faa_file, '--threads', str(threads)])

			if not (os.path.isfile(is_faa_path)):
				raise AssertionError()
			if not (os.path.isfile(is_descriptions_file)):
				raise AssertionError()
			listing_handle.write('isfinder\t' + is_descriptions_file + '\t' + is_faa_file + '\tNA\n')
			_run_cmd(['rm', '-rf', is_faa_path, is_faa_proc_path])
		except Exception as e:
			sys.stderr.write('Issues setting up ISFinder database.\n')
			issues_handle.write('Issues setting up ISFinder database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up MOB-suite plasmid associated proteins database!')
			ms_tar_path = download_path + 'data.tar.gz'
			ms_faa_path = download_path + 'MOB-suite_proteins.faa'
			ms_descriptions_file = download_path + 'MOB-suite_proteins.txt'
			os.mkdir(download_path + 'MOB-suite_files/')
			_run_cmd(['tar', '-xzf', ms_tar_path, '-C', download_path + 'MOB-suite_files/'])
			_concat_files(sorted(glob.glob(download_path + 'MOB-suite_files/data/*.faa')), ms_faa_path)
			_run_cmd(['diamond', 'makedb', '--in', ms_faa_path, '-d', mobs_faa_file, '--threads', str(threads)])
			mdf_handle = open(ms_descriptions_file, 'w')
			with _open_fasta_text(ms_faa_path) as oif:
				for rec in SeqIO.parse(oif, 'fasta'):
					mdf_handle.write(rec.id + '\t' + rec.description + '\n')
			mdf_handle.close()
			if not (os.path.isfile(mobs_faa_file)):
				raise AssertionError()
			if not (os.path.isfile(ms_descriptions_file)):
				raise AssertionError()
			listing_handle.write('mobsuite\t' + ms_descriptions_file + '\t' + mobs_faa_file + '\tNA\n')
			_run_cmd(['rm', '-rf', ms_faa_path, ms_tar_path, download_path + 'MOB-suite_files/'])
		except Exception as e:
			sys.stderr.write('Issues setting up MOB-suite proteins database.\n')
			issues_handle.write('Issues setting MOB-suite proteins database.\n')
			sys.stderr.write(traceback.format_exc())
			sys.stderr.write(str(e) + '\n')

		try:
			print('Setting up Hug et al. 2016 universal ribosomal proteins database!')
			if not (os.path.isfile(unip_hmm_file)):
				raise AssertionError()

			_run_cmd(['hmmpress', '-f', unip_hmm_file])
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

		_setup_defensefinder_db(download_path, download_release_info_path, listing_handle, issues_handle)

	listing_handle.close()
	issues_handle.close()

	print("Done setting up annotation databases!")
	print(f"Information on final files used by zol can be found at:\n{listing_file}\n--------------------")
	print(f"Information on database version/timestamp information can be found at:\n{download_release_info_path}\n--------------------")

	sys.exit(0)

if __name__ == '__main__': 
	setup_annot_dbs()