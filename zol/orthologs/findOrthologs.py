#!/usr/bin/env python

### Program: findOrthologs.py
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
from operator import attrgetter
from zol import util
import shutil

split_diamond_results_prog = os.path.abspath(os.path.dirname(__file__) + '/') + '/splitDiamondResults'
rbh_prog = os.path.abspath(os.path.dirname(__file__) + '/') + '/runRBH'

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""               
	Program: findOrthologs.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Performs reflexive alignment of proteins across different protein-sets/genomes and determines orthologs, 
	co-orthologs, and in-paralogs.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-p', '--proteome_dir',
						help='Path to directory with genomes. Should end with .faa, .faa, or .fasta', required=True)
	parser.add_argument('-o', '--output_dir', help='Path to output directory.', required=True)
	parser.add_argument('-e', '--evalue', type=float, help='E-value cutoff for determining orthologs.', required=False, default=0.001)
	parser.add_argument('-i', '--identity', type=float, help='Percent identity cutoff for determining orthologs.', required=False, default=30.0)
	parser.add_argument('-q', '--coverage', type=float, help='Bi-directional coverage cutoff for determining orthologs.', required=False, default=50.0)
	parser.add_argument('-c', '--cpus', type=int, help='Maximum number of CPUs to use. Default is 1.', default=1,
						required=False)
	args = parser.parse_args()
	return args


def refactorProteomes(inputs):
	sample_name, prot_file, updated_prot_file, original_naming_file, logObject = inputs
	try:
		updated_prot_handle = open(updated_prot_file, 'w')
		original_naming_handle = open(original_naming_file, 'w')
		with open(prot_file) as opf:
			for i, rec in enumerate(SeqIO.parse(opf, 'fasta')):
				original_naming_handle.write(
					sample_name + '\t' + prot_file + '\t' + 'Protein_' + str(i) + '\t' + rec.id + '\n')
				updated_prot_handle.write('>' + sample_name + '|' + 'Protein_' + str(i) + '\n' + str(rec.seq) + '\n')
		original_naming_handle.close()
		updated_prot_handle.close()
	except Exception as e:
		logObject.error('Error with refactoring proteome file.')
		logObject.error(e)
		sys.exit(1)


def oneVsAllParse(inputs):
	sample, samp_algn_res, samp_forth_res, identity_cutoff, coverage_cutoff, logObject = inputs

	rbh_cmd = [rbh_prog, samp_algn_res, str(identity_cutoff), str(coverage_cutoff), '>', samp_forth_res]
	try:
		subprocess.call(' '.join(rbh_cmd), shell=True, stdout=subprocess.DEVNULL,
						stderr=subprocess.DEVNULL, executable='/bin/bash')
		assert (os.path.isfile(samp_forth_res))
	except Exception as e:
		logObject.error("Issue with running: %s" % ' '.join(rbh_cmd))
		logObject.error(e)
		raise RuntimeError(e)

	os.system('rm -f %s %s' % (samp_algn_res))

def createFinalResults(concat_result_file, proteome_listing_file, original_naming_file, result_tab_file,
					   result_mat_file, logObject):
	try:
		samples = []
		with open(proteome_listing_file) as oplf:
			for line in oplf:
				line = line.strip()
				ls = line.split('\t')
				samples.append('.'.join(ls[1].split('/')[-1].split('.')[:-1]))

		mapping = {}
		for f in os.listdir(original_naming_file):
			name_file = original_naming_file + f
			with open(name_file) as onf:
				for line in onf:
					line = line.strip()
					rn_sample, og_sample, rn_protein, og_protein = line.split('\t')
					og_sample = '.'.join(og_sample.split('/')[-1].split('.')[:-1])
					mapping[rn_sample + '|' + rn_protein] = tuple([og_sample, og_protein])

		result_tab_handle = open(result_tab_file, 'w')
		result_mat_handle = open(result_mat_file, 'w')

		sorted_samples = sorted(samples)
		result_mat_handle.write('Sample\t' + '\t'.join(sorted_samples) + '\n')
		with open(concat_result_file) as ocrf:
			for line in ocrf:
				line = line.strip()
				ls = line.split('\t')
				og = ls[0]
				prots = ls[1:]
				samp_prots = defaultdict(set)
				for p in sorted(prots):
					s, op = mapping[p]
					samp_prots[s].add(op)
					result_tab_handle.write(og + '\t' + s + '\t' + op + '\n')
				printlist = [og]
				for s in sorted_samples:
					printlist.append(', '.join(sorted(samp_prots[s])))
				result_mat_handle.write('\t'.join(printlist) + '\n')
		result_tab_handle.close()
		result_mat_handle.close()
	except Exception as e:
		logObject.error(e)
		sys.stderr.write(e)
		sys.exit(1)

def findOrthologs():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE INPUTS
	"""
	myargs = create_parser()

	proteome_dir = os.path.abspath(myargs.proteome_dir) + '/'
	outdir = os.path.abspath(myargs.output_dir) + '/'
	cpus = myargs.cpus
	coverage_cutoff = myargs.coverage
	identity_cutoff = myargs.identity
	evalue_cutoff = myargs.evalue

	if not os.path.isdir(outdir):
		os.system('mkdir %s' % outdir)
	else:
		sys.stderr.write('Overwriting results file in 5 seconds! Stop now if you do not wish to proceed!\n')

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)

	version_string = util.parseVersionFromSetupPy()
	sys.stdout.write('Running version: %s\n' % version_string)
	logObject.info("Running version: %s" % version_string)

	parameters_file = outdir + 'Command_Issued.txt'
	sys.stdout.write("Appending command issued for future records to: %s\n" % parameters_file)
	sys.stdout.write("Logging more details at: %s\n" % log_file)
	logObject.info("\nNEW RUN!!!\n**************************************")
	logObject.info('Running version %s' % version_string)
	logObject.info("Appending command issued for future records to: %s" % parameters_file)

	parameters_handle = open(parameters_file, 'a+')
	parameters_handle.write(' '.join(sys.argv) + '\n')
	parameters_handle.close()

	local_proteome_dir = outdir + 'Proteomes/'
	if not os.path.isdir(local_proteome_dir):
		os.system('mkdir %s' % local_proteome_dir)

	proteome_name_dir = outdir + 'Original_Naming_of_Proteomes/'
	if not os.path.isdir(proteome_name_dir):
		os.system('mkdir %s' % proteome_name_dir)

	checkpoint_dir = outdir + 'Checkpoint_Files/'
	if not os.path.isdir(checkpoint_dir):
		os.system('mkdir %s' % checkpoint_dir)

	sys.stdout.write("--------------------\nStep 0\n--------------------\nProcessing proteomes\n")
	logObject.info("--------------------\nStep 0\n--------------------\nProcessing proteomes")

	proteome_listing_file = outdir + 'Proteome_Listing.txt'
	step0_checkpoint_file = checkpoint_dir + 'Step0_Checkpoint.txt'
	if not os.path.isfile(step0_checkpoint_file):
		try:
			proteome_listing_handle = open(proteome_listing_file, 'w')
			proteome_counter = 0
			assert (os.path.isdir(proteome_dir))
			refactor_proteomes = []
			for f in sorted(os.listdir(proteome_dir)):
				suffix = f.split('.')[-1]
				if not suffix in set(['fa', 'faa', 'fasta']): continue
				proteome_counter += 1
				prot_file = proteome_dir + f
				if not util.is_fasta(prot_file):
					logObject.warning('File %s not in valid FASTA format and ignored.')
					continue
				sample_name = 'Proteome_' + str(proteome_counter)
				updated_prot_file = local_proteome_dir + sample_name + '.faa'
				original_naming_file = proteome_name_dir + sample_name + '.txt'
				refactor_proteomes.append([sample_name, prot_file, updated_prot_file, original_naming_file, logObject])
				proteome_listing_handle.write(sample_name + '\t' + prot_file + '\t' + updated_prot_file + '\n')

			with concurrent.futures.ThreadPoolExecutor(max_workers=cpus) as executor:
				executor.map(refactorProteomes, refactor_proteomes)


		except Exception as e:
			logObject.error(e)
			logObject.error('Difficulties with parsing directory of proteomes.')
			sys.exit(1)
		proteome_listing_handle.close()
		os.system('touch %s' % step0_checkpoint_file)

	"""
	START WORKFLOW
	"""

	# Step 1: Concatenate all proteins into single multi-FASTA file
	sys.stdout.write("--------------------\nStep 1\n--------------------\nConcatenating proteins into multi-FASTA\n")
	logObject.info("--------------------\nStep 1\n--------------------\nConcatenating proteins into multi-FASTA.")
	concat_faa = outdir + 'All_Proteins.faa'
	if not os.path.isfile(concat_faa):
		cf_handle = open(concat_faa, 'w')
		for f in os.listdir(local_proteome_dir):
			with open(local_proteome_dir + f) as olf:
				for rec in SeqIO.parse(olf, 'fasta'):
					cf_handle.write('>' + rec.id + '\n' + str(rec.seq) + '\n')
		cf_handle.close()

	########### Perform Intra-Cluster RBH using reflexive DIAMOND alignment + SwiftOrtho's find_orth.py

	# Step 2: Perform reflexive alignment via Diamond
	sys.stdout.write("--------------------\nStep 2\n--------------------\nRunning reflexive alignment using Diamond\n")
	logObject.info("--------------------\nStep 2\n--------------------\nRunning reflexive alignment using Diamond")
	step2_checkpoint_file = checkpoint_dir + 'Step2_Checkpoint.txt'

	align_res_dir = outdir + 'Alignments/'
	forth_res_dir = outdir + 'Ortholog_Inference_Results/'

	sample_prots = defaultdict(set)
	sample_inputs = set([])
	sample_listing_file = outdir + 'sample_listing.txt'
	sample_listing_handle = open(sample_listing_file, 'w')
	visited = set([])
	try:
		with open(concat_faa) as ocf:
			for rec in SeqIO.parse(ocf, 'fasta'):
				sample = rec.id.split('|')[0]
				sample_prots[sample].add(rec.id)
				samp_algn_res = align_res_dir + sample + '.out'
				samp_forth_res = forth_res_dir + sample + '.abc-like'
				sample_inputs.add(tuple([sample, samp_algn_res, samp_forth_res]))
				if not sample in visited:
					sample_listing_handle.write(sample + '\t' + align_res_dir + sample + '.out\n')
				visited.add(sample)
	except Exception as e:
		raise RuntimeError(e)
	sample_listing_handle.close()

	if not os.path.isfile(step2_checkpoint_file):
		os.system('mkdir %s' % align_res_dir)
		db_path = outdir + 'All_Proteins.dmnd'
		diamond_makedb_cmd = ['diamond', 'makedb', '--ignore-warnings', '--in', concat_faa, '-d', db_path]

		try:
			subprocess.call(' '.join(diamond_makedb_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL, executable='/bin/bash')
			assert (os.path.isfile(db_path))
		except Exception as e:
			logObject.error("Issue with running: %s" % ' '.join(diamond_makedb_cmd))
			logObject.error(e)
			raise RuntimeError(e)

		alignment_result_file = outdir + 'Reflexive_Alignment.out'
		diamond_blastp_cmd = ['diamond', 'blastp', '--ignore-warnings', '-d', db_path, '-q', concat_faa, '-o',
							  alignment_result_file, '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length',
							  'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
							  'qcovhsp', '-k0', '--very-sensitive', '-e', str(evalue_cutoff), '--masking', '0', '-p',
							  str(cpus)]
		try:
			subprocess.call(' '.join(diamond_blastp_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL, executable='/bin/bash')
			assert (os.path.isfile(alignment_result_file))
		except Exception as e:
			logObject.error("Issue with running: %s" % ' '.join(diamond_blastp_cmd))
			logObject.error(e)
			raise RuntimeError(e)

		split_diamond_cmd = [split_diamond_results_prog, alignment_result_file, sample_listing_file]

		try:
			subprocess.call(' '.join(split_diamond_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL, executable='/bin/bash')
		except Exception as e:
			logObject.error("Issue with running: %s" % ' '.join(split_diamond_cmd))
			logObject.error(e)
			raise RuntimeError(e)

		os.system('touch %s' % step2_checkpoint_file)

	sys.stdout.write("--------------------\nStep 3\n--------------------\nRunning RBH.\n")
	logObject.info("--------------------\nStep 3\n--------------------\nRunning RBH.")
	scaled_find_orthos_result_file = outdir + 'All.normalized.abc-like'
	step3_checkpoint_file = checkpoint_dir + 'Step3_Checkpoint.txt'
	if not os.path.isfile(step3_checkpoint_file):
		os.system('mkdir %s' % forth_res_dir)

		find_orthos_result_file = outdir + 'All.abc-like'

		parallelize_inputs = []
		for inp_tup in sample_inputs:
			parallelize_inputs.append(list(inp_tup) + [identity_cutoff, coverage_cutoff, logObject])

		with concurrent.futures.ThreadPoolExecutor(max_workers=cpus) as executor:
			executor.map(oneVsAllParse, parallelize_inputs)

		os.system('cat %s* > %s' % (forth_res_dir, find_orthos_result_file))

		scaled_find_orthos_result_handle = open(scaled_find_orthos_result_file, 'w')
		qh_pairs_accounted = set([])
		qshs_rbh = defaultdict(int)
		qshs_rbh_sum = defaultdict(float)
		with open(find_orthos_result_file) as osforf:
			for line in osforf:
				line = line.strip()
				query, query_samp, hit, hit_samp, bs = line.split('\t')
				qh_pair = tuple(sorted([query, hit]))
				if qh_pair in qh_pairs_accounted: continue
				qh_pairs_accounted.add(qh_pair)
				sample_pair = tuple(sorted([query_samp, hit_samp]))
				qshs_rbh[sample_pair] += 1
				qshs_rbh_sum[sample_pair] += float(bs)

		qshs_rbh_avg = {}
		for sp in qshs_rbh:
			qshs_rbh_avg[sp] = qshs_rbh_sum[sp]/qshs_rbh[sp]

		qh_pairs_accounted = set([])
		with open(find_orthos_result_file) as osforf:
			for line in osforf:
				line = line.strip()
				query, query_samp, hit, hit_samp, bs = line.split('\t')
				qh_pair = tuple(sorted([query, hit]))
				if qh_pair in qh_pairs_accounted: continue
				qh_pairs_accounted.add(qh_pair)
				if query_samp == hit_samp:
					scaled_find_orthos_result_handle.write(query + '\t' + hit + '\t' + str(float(bs)/100.0) + '\n')
				else:
					sample_pair = tuple(sorted([query_samp, hit_samp]))
					scaled_find_orthos_result_handle.write(query + '\t' + hit + '\t' + str(float(bs)/qshs_rbh_avg[sample_pair]) + '\n')
		scaled_find_orthos_result_handle.close()

	sys.stdout.write("--------------------\nStep 4\n--------------------\nRun MCL to determine OrthoGroups.\n")
	logObject.info("--------------------\nStep 4\n--------------------\nRun MCL to determine OrthoGroups.")
	tmp_result_file = outdir + 'Tmp_Results.txt'
	step4_checkpoint_file = checkpoint_dir + 'Step4_Checkpoint.txt'
	if not os.path.isfile(step4_checkpoint_file):
		cluster_result_file = outdir + 'MCL_Cluster_Results.txt'
		find_clusters_cmd = ['mcl', scaled_find_orthos_result_file, '--abc', '-I', '1.5', '-te',
								str(cpus), '-o', cluster_result_file]
		try:
			subprocess.call(' '.join(find_clusters_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL, executable='/bin/bash')
			assert (os.path.isfile(cluster_result_file))
		except Exception as e:
			logObject.error("Issue with running: %s" % ' '.join(find_clusters_cmd))
			logObject.error(e)
			raise RuntimeError(e)

		result_handle = open(tmp_result_file, 'w')
		subclust_id = 1
		clustered = set([])
		with open(cluster_result_file) as omrf:
			for line in sorted([l for l in omrf.readlines()]):
				line = line.strip()
				ls = line.split()
				clustered = clustered.union(set(ls))
				result_handle.write('OG_' + str(subclust_id) + '\t' + '\t'.join(ls) + '\n')
				subclust_id += 1

		with open(concat_faa) as ocf:
			fasta_dict = SeqIO.to_dict(SeqIO.parse(ocf, "fasta"))
			for rec in sorted(fasta_dict.values(), key=attrgetter('id')):
				if not rec.id in clustered:
					result_handle.write('OG_' + str(subclust_id) + '\t' + rec.id + '\n')
					subclust_id += 1
		result_handle.close()

		#os.system('rm %s %s' % (alignment_result_file, cluster_result_file))
		#shutil.rmtree(forth_res_dir)
		#shutil.rmtree(align_res_dir)

		os.system('touch %s' % step4_checkpoint_file)

	step5_checkpoint_file = checkpoint_dir + 'Step5_Checkpoint.txt'
	sys.stdout.write("--------------------\nStep 5\n--------------------\nCreate final result files!\n")
	logObject.info("--------------------\nStep 5\n--------------------\nCreate final result files!")

	result_tab_file = outdir + 'Orthogroups.Listing.tsv'
	result_mat_file = outdir + 'Orthogroups.tsv'
	if not os.path.isfile(step5_checkpoint_file):
		os.system('touch %s' % step5_checkpoint_file)
		createFinalResults(tmp_result_file, proteome_listing_file, proteome_name_dir, result_tab_file,
								result_mat_file, logObject)
		os.system('touch %s' % step5_checkpoint_file)

	# DONE!
	sys.stdout.write(
		"--------------------\nDONE!\n--------------------\nOrthogroup by sample matrix can be found at: %s\n" % result_mat_file)
	logObject.info(
		"--------------------\nDONE!\n--------------------\nOrthogroup by sample matrix can be found at: %s" % result_mat_file)


if __name__ == '__main__':
	findOrthologs()

