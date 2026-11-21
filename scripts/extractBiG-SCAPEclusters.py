#!/usr/bin/env python3

"""
Program: extractBiG-SCAPEclusters.py
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
from time import sleep
import shutil
from zol import util
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

	if not (os.path.isdir(antismash_dir) and os.path.isdir(bigscape_dir)):
		sys.stderr.write('Issue with validating that the antismash and bigscape directories provided are valid.\n')
		sys.exit(1)

	if not (util.is_integer(gcf_id)):
		raise RuntimeError('Issue with path to BGC predictions Genbanks listing file.')

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n ")
		sleep(5)

	gbk_dir = outdir + "GeneCluster_GenBanks/"
	util.setup_ready_directory([outdir, gbk_dir])

	"""
	START WORKFLOW
	"""

	# create logging object
	log_file = outdir + 'Progress.log'
	log_object = util.create_logger_object(log_file)
	version_string = util.get_version()
	log_object.info(f'Running zol version {version_string}')
	sys.stderr.write(f'Running zol version {version_string}')

	log_object.info("Saving parameters for future records.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [antismash_dir, bigscape_dir, outdir, gcf_id, type, use_symlink]
	parameter_names = ["AntiSMASH Results Directory", "BiG-SCAPE Results Directory", "Output Directory of Extracted GenBanks",
					   "GCF_ID", "Type of BGC for Focal GCF ID", "Use Symlink instead of Copy"]
	util.log_parameters_to_file(parameters_file, parameter_names, parameter_values)
	log_object.info("Done saving parameters!")

	# Step 1: Parse BiG-SCAPE Results - Get Relevant BGCs and Samples
	# Supports both BiG-SCAPE 1.x and 2.x clustering file formats. The class (v1) / bin (v2) is the name of the
	# folder containing the clustering file.
	gcf_bgcs = set([])
	for cluster_tsv in util.find_bigscape_clustering_files(bigscape_dir):
		cluster_annotation_type = os.path.basename(os.path.dirname(cluster_tsv))
		if type is not None and type != cluster_annotation_type: continue
		_, gcf_to_bgcs = util.parse_bigscape_clustering_file(cluster_tsv)
		gcf_bgcs.update(gcf_to_bgcs.get(gcf_id, []))

	if not (len(gcf_bgcs) >= 1):
		log_object.error('Less than one BGC found to belong to the GCF.')
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
							os.symlink(gbk_file, os.path.join(gbk_dir, filename))
						else:
							shutil.copy2(gbk_file, os.path.join(gbk_dir, filename))

if __name__ == '__main__':
	extractBiGSCAPEclusters()