#!/usr/bin/env python

### Program: extractBiG-SCAPEclusters.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

# BSD 3-Clause License
#
# Copyright (c) 2023, Kalan-Lab
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
from time import sleep
import shutil
from zol import util

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: extractBiG-SCAPEclusters.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--antismash_dir',
						help='antiSMASH results directory.',
						required=True)
	parser.add_argument('-b', '--bigscape_dir',
						help='BiG-SCAPE results directory. Should be generated using directory from "--antismash_dir" as input.',
						required=True)
	parser.add_argument('-o', '--output_dir', help='Parent output/workspace directory.', required=True)
	parser.add_argument('-g', '--gcf_id', help='GCF ID of interest.', required=True)
	parser.add_argument('-t', '--type',
						help='Annotation type of BiG-SCAPE, to make sure GCFs with same IDs across different annotation classes are not used (not required but recommended).',
						required=False, default=None)
	parser.add_argument('-s', '--use_symlink', action='store_true', help='Use symlink instead of performing a deep copy of the file.', default=False, required=False)
	args = parser.parse_args()
	return args


def extractBiGSCAPEclusters():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE INPUTS
	"""
	myargs = create_parser()

	antismash_dir = os.path.abspath(myargs.antismash_dir) + '/'
	bigscape_dir = os.path.abspath(myargs.bigscape_dir) + '/'
	outdir = os.path.abspath(myargs.output_dir) + '/'

	gcf_id = myargs.gcf_id
	type = myargs.type
	use_symlink = myargs.use_symlink

	try:
		assert (os.path.isdir(antismash_dir) and os.path.isdir(bigscape_dir))
	except:
		sys.stderr.write('Issue with validating that the antismash and bigscape directories provided are valid.\n')
		sys.exit(1)

	try:
		assert (util.is_integer(gcf_id))
	except:
		raise RuntimeError('Issue with path to BGC predictions Genbanks listing file.')

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n ")
		sleep(5)

	gbk_dir = outdir + "GeneCluster_GenBanks/"
	util.setupReadyDirectory([outdir, gbk_dir])

	"""
	START WORKFLOW
	"""

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)
	version_string = util.parseVersionFromSetupPy()
	logObject.info('Running zol version %s' % version_string)
	sys.stderr.write('Running zol version %s' % version_string)

	logObject.info("Saving parameters for future records.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [antismash_dir, bigscape_dir, outdir, gcf_id, type, use_symlink]
	parameter_names = ["AntiSMASH Results Directory", "BiG-SCAPE Results Directory", "Output Directory of Extracted GenBanks",
					   "GCF_ID", "Type of BGC for Focal GCF ID", "Use Symlink instead of Copy"]
	util.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Step 1: Parse BiG-SCAPE Results - Get Relevant BGCs and Samples
	gcf_bgcs = set([])
	for dirpath, dirnames, files in os.walk(bigscape_dir):
		for filename in files:
			if filename.endswith(".tsv") and "clustering" in filename:
				cluster_tsv = os.path.join(dirpath, filename)
				cluster_annotation_type = cluster_tsv.split('/')[-2]
				if type != None and type != cluster_annotation_type: continue
				with open(cluster_tsv) as oct:
					for i, line in enumerate(oct):
						if i == 0: continue
						line = line.strip()
						bgc_id, clust_gcf_id = line.split('\t')
						if not clust_gcf_id == gcf_id: continue
						gcf_bgcs.add(bgc_id)

	try:
		assert (len(gcf_bgcs) >= 1)
	except:
		logObject.error('Less than one BGC found to belong to the GCF.')
		sys.stderr.write('Error: Less than one BGC found to belong to the GCF.\n')
		sys.exit(1)

	# Step 2: Parse through AntiSMASH Results and Copy over to Output Directory
	for dirpath, dirnames, files in os.walk(antismash_dir):
		for filename in files:
			if filename.endswith('.gbk'):
				gbk_file = os.path.join(dirpath, filename)
				if '.region' in gbk_file:
					bgc = filename.split('.gbk')[0]
					if bgc in gcf_bgcs:
						if use_symlink:
							os.system('ln -s %s %s' % (gbk_file, gbk_dir))
						else:
							shutil.copyfile(gbk_file, gbk_dir)

if __name__ == '__main__':
	extractBiGSCAPEclusters()