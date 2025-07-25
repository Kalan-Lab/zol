#!/usr/bin/env python3

"""
Program: selectSpecificGeneClusters.py
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
from zol import util
import shutil
import time

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: selectSpecificGeneClusters.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	
	Script to automatically subset GenBank's from a fai results directory.

	Users can manually assess and filter the XLSX spreadsheet for instances that meet certain requirements beyond 
	standard filters in fai. Afterwards, they can select a list of sample names (first column in the sheet 
	"Genome Wide - Results") or individual gene-cluster instances (second column in the sheet "Gene Cluster - Results"). 

	Copying and pasting such IDs into a simple text file (one ID per line) and providing the path to the general fai 
	results directory to this script, will allow subsampling only for gene clusters which meet the desired stringency 
	thresholds and copy/symlink them to an output folder. 
								  
	This output folder, with multiple homologous gene cluster representations in GenBank format, can then be used as an 
	input to tools such as zol or clinker.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--fai_results_dir', help='Path to the fai results directory.', required=True)
	parser.add_argument('-s', '--selections_file', help='File lisitng samples/genomes or gene-cluster instance GenBanks to keep. One per line.', required=True)
	parser.add_argument('-o', '--output_dir', help='Path to output directory where subset of selected gene cluster GenBanks will be written.', required=True)
	parser.add_argument('-t', '--type', help='Whether the selection listing is being made at the sample level ("sample") or for individual gene cluster instances ("instance"). [Default is "sample"].', required=True)
	parser.add_argument('-l', '--symlink', action='store_true', help='Soft symlink gene cluster GenBanks to the output directory instead of fully copying them over.', default=False, required=False)
	args = parser.parse_args()

	return args

def selectGbks():
	"""
	PARSE INPUTS
	"""
	myargs = create_parser()

	fai_results_dir = os.path.abspath(myargs.fai_results_dir) + '/'
	selections_file = myargs.selections_file
	outdir = os.path.abspath(myargs.output_dir) + '/'
	type = myargs.type
	symlink = myargs.symlink

	try:
		assert(os.path.isdir(fai_results_dir) and os.path.isfile(selections_file))
	except Exception as e:
		sys.stderr.write('Error validating input directory of homologous gene-clusters or the query fasta exists!\n')

	try:
		assert(type in set(["instance", "sample"]))
	except Exception as e:
		sys.stderr.write('Error, the type of listing provided is neither "sample" or "instance"\n')
	
	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n")
		time.sleep(5)
	util.setup_ready_directory([outdir])
	
	homologous_gbk_dir = fai_results_dir + 'Final_Results/Homologous_Gene_Cluster_GenBanks/'
	gc_gbk_dir = fai_results_dir + 'GC_Segment_Processing/GeneCluster_Genbanks/'

	if type == 'sample':
		selected = set([])
		with open(selections_file) as osf:
			for line in osf:
				line = line.strip()
				selected.add(line)
				
		for gbk in os.listdir(homologous_gbk_dir):
			sample = '.gbk'.join(gbk.split('.gbk')[:-1]).split('_fai-gene-cluster')[0]
			if sample in selected:
				if symlink:
					symlink_file = outdir + sample + '.gbk'
					os.symlink(homologous_gbk_dir + gbk, symlink_file)
				else:
					shutil.copy2(homologous_gbk_dir + gbk, outdir)
					
	elif type == 'instance':
		selected = set([])
		with open(selections_file) as osf:
			for line in osf:
				line = line.strip().split('/')[-1]
				selected.add(line)
				
		for gbk in os.listdir(homologous_gbk_dir):
			if gbk in selected:
				if symlink:
					symlink_file = outdir + gbk
					os.symlink(homologous_gbk_dir + gbk, symlink_file)
				else:
					shutil.copy2(homologous_gbk_dir + gbk, outdir)

		for gbk in os.listdir(gc_gbk_dir):
			if gbk in selected:
				if symlink:
					symlink_file = outdir + gbk
					os.symlink(gc_gbk_dir + gbk, symlink_file)
				else:
					shutil.copy2(gc_gbk_dir + gbk, outdir)

if __name__ == '__main__':
	selectGbks()
