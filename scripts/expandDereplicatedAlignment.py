#!/usr/bin/env python3

import os
import sys
import argparse
from zol import util
from Bio import SeqIO

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: expandDereplicatedAlignment.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	
	Script to use query protein sequences for some orthogroup to search a set of homologous-gene clusters and produce
	a multiple sequence alignment (MSA). Can optionally run FastTree for gene-phylogeny construction and HyPhy's
	GARD and FUBAR for site-specific selection analysis. 
	
	This is primarily intended for creating a comprehensive MSA after running zol in "dereplicated" mode and identifying
	distinct ortholog groups. It is much more computationally tractable than running zol with a full set of highly
	redundant gene-clusters.
	
	Note query files will not be included in resulting MSA/gene-tree.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--input_dir', help='Directory with orthologous/homologous locus-specific GenBanks.\nFiles must end with ".gbk", ".gbff", or ".genbank".', required=False, default=None)
	parser.add_argument('-q', '--query_prots', help='Query protein(s) provided in FASTA format.', required=True)
	parser.add_argument('-o', '--output_dir', help='Path to output directory.', required=True)
	parser.add_argument('-mi', '--min_identity', type=float, help='Minimum identity to a known instance (sequence in query). Default is 95.0.', required=False, default=95.0)
	parser.add_argument('-mc', '--min_coverage', type=float, help='Minimum query coverage of a known instance (sequence in query). Default is 95.0.', required=False, default=95.0)
	parser.add_argument('-c', '--threads', type=int, help='The number of threads to use. Default is 1.', required=False, default=1)
	parser.add_argument('-f', '--run_fasttree', action='store_true', help='Run FastTree for approximate maximum-likelihood phylogeny generation of the orthogroup.', required=False, default=False)
	args = parser.parse_args()

	return args

def expandOg():
	"""
	PARSE INPUTS
	"""
	myargs = create_parser()

	input_dir = os.path.abspath(myargs.input_dir) + '/'
	query_fasta = myargs.query_prots
	outdir = os.path.abspath(myargs.output_dir) + '/'
	min_identity = myargs.min_identity
	min_coverage = myargs.min_coverage
	threads = myargs.threads
	run_fasttree = myargs.run_fasttree

	try:
		assert(os.path.isdir(input_dir) and util.is_fasta(query_fasta))
	except:
		sys.stderr.write('Error validating input directory of homologous gene-clusters or the query fasta exists!\n')

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n")
		#sleep(5)
	else:
		os.mkdir(outdir)

	db_faa = outdir + 'Database.faa'
	db_dmnd = outdir + 'Database.dmnd'

	db_handle = open(db_faa, 'w')
	for f in os.listdir(input_dir):
		if util.is_genbank(input_dir + f):
			with open(input_dir + f) as of:
				for rec in SeqIO.parse(of, 'genbank'):
					for feature in rec.features:
						if not feature.type == 'CDS': continue
						lt = feature.qualifiers.get('locus_tag')[0]
						translation = feature.qualifiers.get('translation')[0]
						db_handle.write('>' + '.'.join(f.split('.')[:-1]) + '|' + lt + '\n' + translation + '\n')
	db_handle.close()

	align_result_file = outdir + 'DIAMOND_Results.txt'
	diamond_cmd = ['diamond', 'makedb', '--ignore-warnings', '--in', db_faa, '-d', db_dmnd, ';',
				   'diamond', 'blastp', '--ignore-warnings', '--threads', str(threads), '--very-sensitive', '--query',
				   query_fasta, '--db', db_dmnd, '--outfmt', '6', 'qseqid', 'sseqid',
				   'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
				   'bitscore', 'qcovhsp', 'scovhsp', '-k0', '--out', align_result_file, '--evalue', '1e-3']

	util.run_cmd_via_subprocess(diamond_cmd, check_files = [align_result_file])
	
	matching_lts = set([])
	with open(align_result_file) as oarf:
		for line in oarf:
			line = line.strip()
			qseqid, sseqid, pident, length, mismatch, gapone, qstart, qend, sstart, send, evalue, bitscore, qcovhsp, scovhsp  = line.split('\t')
			pident = float(pident)
			qcovhsp = float(qcovhsp)
			if qcovhsp >= min_coverage and pident >= min_identity:
				matching_lts.add(sseqid)

	orthogroup_seqs_faa = outdir + 'OrthoGroup.faa'
	orthogroup_seqs_handle = open(orthogroup_seqs_faa, 'w')
	with open(db_faa) as odf:
		for rec in SeqIO.parse(odf, 'fasta'):
			if rec.id in matching_lts:
				orthogroup_seqs_handle.write('>' + rec.id + '\n' + str(rec.seq) + '\n')
	orthogroup_seqs_handle.close()

	orthogroup_seqs_msa = outdir + 'OrthoGroup.msa.faa'
	muscle_cmd = ['muscle', '-super5', orthogroup_seqs_faa, '-output', orthogroup_seqs_msa, '-amino', '-threads', str(threads), '-perturb', '12345']
	util.run_cmd_via_subprocess(muscle_cmd, check_files = [orthogroup_seqs_msa])

	phylogeny_tre = outdir + 'OrthoGroup.tre'
	if run_fasttree:
		fasttree_cmd = ['FastTree', orthogroup_seqs_msa, '>', phylogeny_tre]
		util.run_cmd_via_subprocess(fasttree_cmd, check_files = [phylogeny_tre])

	sys.exit(0)

if __name__ == '__main__':
	expandOg()
