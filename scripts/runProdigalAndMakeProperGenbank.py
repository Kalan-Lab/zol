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
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import defaultdict
from zol import util
from operator import itemgetter

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

	# check if uncompression needed
	try:
		if input_genomic_fasta_file.endswith('.gz'):
			os.system(f'cp {input_genomic_fasta_file} {outdir}')
			updated_genomic_fasta_file = outdir + input_genomic_fasta_file.split('/')[-1].split('.gz')[0]
			os.system(f'gunzip {updated_genomic_fasta_file}.gz')
			input_genomic_fasta_file = updated_genomic_fasta_file
		assert(util.is_fasta(input_genomic_fasta_file))
	except Exception as e:
		raise RuntimeError('Input genomic assembly does not appear to be in FASTA format.')

	# Step 1: Run Prodigal (if needed)
	og_prod_pred_prot_file = outdir + sample_name + '.original_predicted_proteome'
	prodigal_cmd = []
	if gene_calling_method == 'pyrodigal':
		prodigal_cmd = ['pyrodigal', '-i', input_genomic_fasta_file, '-a', og_prod_pred_prot_file]
	elif gene_calling_method == 'prodigal':
		prodigal_cmd = ['prodigal', '-i', input_genomic_fasta_file, '-a', og_prod_pred_prot_file]
	elif gene_calling_method == 'prodigal-gv':
		prodigal_cmd = ['prodigal-gv', '-p', 'meta', '-i', input_genomic_fasta_file, '-a', og_prod_pred_prot_file]
	else:
		sys.stderr.write('The gene-calling method selected is not a valid option. Has to be either: prodigal, pyrodigal, or prodigal-gv.\n')
		sys.exit(1)

	if meta_mode and not gene_calling_method == 'prodigal-gv':
		prodigal_cmd += ['-p', 'meta']

	util.run_cmd_via_subprocess(prodigal_cmd, check_files = [og_prod_pred_prot_file])

	# Step 2: Process Prodigal predicted proteome and create polished version with locus tags
	pc_prod_pred_prot_file = outdir + sample_name + '.faa'
	pc_prod_pred_prot_handle = open(pc_prod_pred_prot_file, 'w')
	scaffold_prots = defaultdict(list)
	prot_sequences = {}
	prot_locations = {}
	with open(og_prod_pred_prot_file) as ooppf:
		for i, rec in enumerate(SeqIO.parse(ooppf, 'fasta')):
			if (i+1) < 10:
				pid = '00000'+str(i+1)
			elif (i+1) < 100:
				pid = '0000'+str(i+1)
			elif (i+1) < 1000:
				pid = '000'+str(i+1)
			elif (i+1) < 10000:
				pid = '00'+str(i+1)
			elif (i+1) < 100000:
				pid = '0'+str(i+1)
			else:
				pid = str(i+1)
			new_prot_id = locus_tag + '_' + pid
			scaffold = '_'.join(rec.id.split('_')[:-1])
			start = int(rec.description.split(' # ')[1])
			end = int(rec.description.split(' # ')[2])
			direction = int(rec.description.split(' # ')[3])
			scaffold_prots[scaffold].append([new_prot_id, start])
			prot_locations[new_prot_id] = [scaffold, start, end, direction]
			prot_sequences[new_prot_id] = str(rec.seq)
			dir_str = '+'
			if direction == -1: dir_str = '-'
			pc_prod_pred_prot_handle.write('>' + new_prot_id + ' ' + scaffold + ' ' + str(start) + ' ' + str(end) + ' ' + dir_str + '\n' + str(rec.seq) + '\n')
	pc_prod_pred_prot_handle.close()

	print(scaffold_prots.keys())

	# Step 3: Process Prodigal predicted genbank and create polished version with locus tags and protein + nucleotide
	# sequences
	pc_prod_genbank_file = outdir + sample_name + '.gbk'
	pc_prod_genbank_handle = open(pc_prod_genbank_file, 'w')
	with open(input_genomic_fasta_file) as oigff:
		for fasta_rec in SeqIO.parse(oigff, 'fasta'):
			seq = fasta_rec.seq
			record = SeqRecord(seq, id=fasta_rec.id, name=fasta_rec.id, description='')
			record.annotations['molecule_type'] = 'DNA'
			feature_list = []
			for prot in sorted(scaffold_prots[fasta_rec.id], key=itemgetter(1)):
				pstart = prot_locations[prot[0]][1]
				pend = prot_locations[prot[0]][2]
				pstrand = prot_locations[prot[0]][3]
				feature = SeqFeature(FeatureLocation(start=pstart-1, end=pend, strand=pstrand), type='CDS')
				feature.qualifiers['locus_tag'] = prot[0]
				feature.qualifiers['translation'] = prot_sequences[prot[0]].rstrip('*')
				feature_list.append(feature)
			record.features = feature_list
			SeqIO.write(record, pc_prod_genbank_handle, 'genbank')
	pc_prod_genbank_handle.close()

	sys.exit(0)
	
if __name__ == '__main__':
	prodigalAndReformat()
