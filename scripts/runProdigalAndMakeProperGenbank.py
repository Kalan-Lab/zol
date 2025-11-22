#!/usr/bin/env python3

"""
Program: runProdigalAndMakeProperGenbank.py
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

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: runProdigalAndMakeProperGenbank.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Program to run prodigal gene calling and then create a genbank file with CDS features. Final genome-wide predicted 
	proteomes + genbank will be edited to have locus tags provided with the "-l" option. 
	
	11/17/2022 - can now use pyrodigal instead of progial.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--input_genomic_fasta', help='Path to genomic assembly in FASTA or GenBank (needed for Eukaryotes) format.', required=True)
	parser.add_argument('-o', '--output_directory', help='Path to output directory. Should already be created!', required=True)
	parser.add_argument('-s', '--sample_name', help='Sample name', default='Sample', required=False)
	parser.add_argument('-l', '--locus_tag', help='Locus tag', default='AAA', required=False)
	parser.add_argument('-gcm', '--gene_calling_method', help='Method to use for gene calling. Options are: pyrodigal, prodigal,\nor prodigal-gv. [Default is pyrodigal].', required=False, default='pyrodigal')
	parser.add_argument('-m', '--meta_mode', action='store_true', help='Use meta mode instead of single for pyrodigal/prodigal. Automatically turned on if prodigal-gv is requested.', default=False, required=False)

	args = parser.parse_args()
	return args

def prodigalAndReformat():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE REQUIRED INPUTS
	"""
	myargs = create_parser()

	input_genomic_fasta_file = os.path.abspath(myargs.input_genomic_fasta)
	outdir = os.path.abspath(myargs.output_directory) + '/'
	gene_calling_method = (myargs.gene_calling_method).lower()
	meta_mode = myargs.meta_mode

	try:
		assert(os.path.isfile(input_genomic_fasta_file))
	except Exception as e:
		raise RuntimeError('Issue with input genomic assembly file path.')

	if not os.path.isdir(outdir):
		sys.stderr.write("Output directory does not exist! Please create and retry program.")

	"""
	PARSE OPTIONAL INPUTS
	"""

	sample_name = myargs.sample_name
	locus_tag = myargs.locus_tag

	"""
	START WORKFLOW
	"""

	# Use the core implementation from util.py to avoid code duplication
	try:
		util.prodigal_and_reformat_core(
			input_genomic_fasta_file,
			outdir,
			sample_name,
			locus_tag,
			gene_calling_method=gene_calling_method,
			meta_mode=meta_mode,
			log_object=None,
		)
		sys.exit(0)
	except Exception as e:
		sys.stderr.write(f'Error during gene calling and reformatting: {str(e)}\n')
		sys.exit(1)
	
if __name__ == '__main__':
	prodigalAndReformat()
