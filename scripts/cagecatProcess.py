#!/usr/bin/env python3

### Program: cagecatProcess.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

# BSD 3-Clause License
#
# Copyright (c) 2022, Kalan-Lab
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
import subprocess
import concurrent.futures
from Bio import SeqIO
from collections import defaultdict
from operator import attrgetter, itemgetter
from zol import util
import shutil

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""               
	Program: cagecatProcess.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	
    CAGECAT (https://cagecat.bioinformatics.nl/) allows identification of homologous gene clusters to 
	a query set of co-located proteins. While the set of clusters should be directly compatible with 
	processing using zol, CAGECAT primarily uses protein_id features in gene cluster GenBank files as 
	protein identifiers whereas zol primarily expects locus_tag. This can be overcome using the "-r"
	option in zol which will create arbitrary locus tags with GenBank files to allow for 
	smooth processing in zol. However, if you wish to retain the protein_id information, this 
	script will take in as input the zip folder following gene cluster extraction from CAGECAT
	and create a directory of similar GenBank files with locus_tag features harboring the same value 
	as the protein_id features in the original GenBank files.
								  
	This is not designed for fungi/eukaryotes currently because CAGECAT CDS feature coordinates 
	do not appear to contain exon info needed for zol.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--extract_clusters_zip', help='Pat.', required=True)
	parser.add_argument('-o', '--output_dir', help='Path to output directory.', required=True)
	args = parser.parse_args()
	return args

def cagecatProcess():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE INPUTS
	"""
	myargs = create_parser()

	extract_clusters_zip_file = os.path.abspath(myargs.extract_clusters_zip) 
	outdir = os.path.abspath(myargs.output_dir) + '/'

	if not os.path.isdir(outdir):
		os.system('mkdir %s' % outdir)
	else:
		sys.stderr.write('Note, output directory exists already! Exiting ...\n')
		sys.exit(1)
	
	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)

	version = util.getVersion()
	sys.stdout.write('Running version: %s\n' % version)
	logObject.info("Running version: %s" % version)

	parameters_file = outdir + 'Command_Issued.txt'
	sys.stdout.write("Appending command issued for future records to: %s\n" % parameters_file)
	sys.stdout.write("Logging more details at: %s\n" % log_file)
	logObject.info("\nNEW RUN!!!\n**************************************")
	logObject.info('Running version %s' % version)
	logObject.info("Appending command issued for future records to: %s" % parameters_file)

	parameters_handle = open(parameters_file, 'a+')
	parameters_handle.write(' '.join(sys.argv) + '\n')
	parameters_handle.close()

	"""
	START WORKFLOW
	"""

	# Step 1: Uncompress zip folder into output directory
	msg = "--------------------\nStep 1\n--------------------\nUncompressing zipped cluster extraction from CAGECAT."
	sys.stdout.write(msg + "\n")
	logObject.info(msg)
	
	cagecat_dir = outdir + 'CAGECAT_Results/'
	cagecat_res_dir = cagecat_dir + 'results/'
	extraction_cmd = ['unzip', extract_clusters_zip_file, '-d', cagecat_dir]
	try:
		print(' '.join(extraction_cmd))
		subprocess.call(' '.join(extraction_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
		assert(os.path.isdir(cagecat_res_dir))
	except Exception as e:
		logObject.error("Issue with running: %s" % ' '.join(extraction_cmd))
		logObject.error(e)
		raise RuntimeError(e)
	
	# Step 2: Process GenBank files and produce final versions with locus_tags
	msg = "--------------------\nStep 1\n--------------------\nProcessing GenBank files and producing final versions."
	sys.stdout.write(msg + "\n")
	logObject.info(msg)

	final_genbank_dir = outdir + 'Processed_GenBank_Files/'
	if not os.path.isdir(final_genbank_dir):
		os.mkdir(final_genbank_dir)

	for f in os.listdir(cagecat_res_dir):
		if f.endswith('.gbk'):
			mod_gbk = final_genbank_dir + f
			mod_gbk_handle = open(mod_gbk, 'w')
			cds_count = 1
			try:
				with open(cagecat_res_dir + f) as ocrf:
					for rec in SeqIO.parse(ocrf, 'genbank'):
						for feat in rec.features:
							if feat.type == 'CDS':
								protein_id = 'CDS_' + str(cds_count)
								try:
									protein_id = feat.qualifiers.get('protein_id')[0]
								except:
									pass
								feat.qualifiers['locus_tag'] = protein_id
								cds_count += 1
						SeqIO.write(rec, mod_gbk_handle, 'genbank')
			except Exception as e:
				logObject.error("Issue with processing GenBank file: %s" % ' '.join(cagecat_res_dir + f))
				logObject.error(e)
				raise RuntimeError(e)
			mod_gbk_handle.close()
			
	# DONE!
	sys.stdout.write("--------------------\nDONE!\n--------------------\nDirectory of processed GenBank files from CAGECAT can be found at: %s\n" % final_genbank_dir)
	logObject.info("--------------------\nDONE!\n--------------------\nDirectory of processed GenBank files from CAGECAT can be found at: %s" % final_genbank_dir)

if __name__ == '__main__':
	cagecatProcess()

