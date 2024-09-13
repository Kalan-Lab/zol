#!/usr/bin/env python3

import os
import sys
import argparse
from zol import util, zol
from Bio import SeqIO

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: annotateProtSeqs.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	
    	Simple script to annotate a protein FASTA using zol's standard approach taken for consensus sequences for ortholog groups.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--protein_faa', help='Path to FASTA of protein sequences.', required=True)
	parser.add_argument('-o', '--output_dir', help='Path to output directory.', required=True)
	parser.add_argument('-c', '--threads', type=int, help='The number of threads to use [Default is 1].', required=False, default=1)
	args = parser.parse_args()

	return args

def annotateProtSeqs():
	"""
	PARSE INPUTS
	"""
	myargs = create_parser()

	protein_faa = myargs.protein_faa
	outdir = os.path.abspath(myargs.output_dir) + '/'
	threads = myargs.threads

	try:
		assert(util.is_fasta(protein_faa))
		protein_faa = os.path.abspath(protein_faa)
	except:
		sys.stderr.write('Could not validate %s as a FASTA file' % protein_faa) 

	util.setupReadyDirectory([outdir])		

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = util.createLoggerObject(log_file)

	# perform annotate consensus sequences	
	annotations = zol.annotateConsensusSequences(protein_faa, outdir, logObject, threads=threads)
		
	header = ['Protein ID', 'KO Annotation (E-value)', 'PGAP Annotation (E-value)',
			  'PaperBLAST Annotation (E-value)', 'CARD Annotation (E-value)', 'IS Finder (E-value)',
			  'MI-BiG Annotation (E-value)', 'VOG Annotation (E-value)',  'VFDB Annotation (E-value)',
			  'Pfam Domains', 'Protein Sequence']

	seqs = {}
	with open(protein_faa) as opf:
		for rec in SeqIO.parse(opf, 'fasta'):
			seqs[rec.id] = str(rec.seq)

	annotation_tsv_file = outdir + 'Annotation_Results.tsv'
	frt_handle = open(annotation_tsv_file, 'w')
	frt_handle.write('\t'.join(header) + '\n')

	for prot in seqs:
		ko_annot = util.gatherAnnotationFromDictForHomoloGroup(prot, 'ko', annotations)
		pgap_annot = util.gatherAnnotationFromDictForHomoloGroup(prot, 'pgap', annotations)
		pb_annot = util.gatherAnnotationFromDictForHomoloGroup(prot, 'paperblast', annotations)
		card_annot = util.gatherAnnotationFromDictForHomoloGroup(prot, 'card', annotations)
		isf_annot = util.gatherAnnotationFromDictForHomoloGroup(prot, 'isfinder', annotations)
		mibig_annot = util.gatherAnnotationFromDictForHomoloGroup(prot, 'mibig', annotations)
		vog_annot = util.gatherAnnotationFromDictForHomoloGroup(prot, 'vog', annotations)
		vfdb_annot = util.gatherAnnotationFromDictForHomoloGroup(prot, 'vfdb', annotations)
		pfam_annots = 'NA'
		if 'pfam' in annotations and prot in annotations['pfam']:
			pfam_annots = '; '.join(annotations['pfam'][prot][0])	
		prot_seq = seqs[prot]
		row = [prot, ko_annot, pgap_annot, pb_annot, card_annot, isf_annot, mibig_annot, vog_annot,
				vfdb_annot, pfam_annots, prot_seq]
		row = [str(x) for x in row]
		frt_handle.write('\t'.join(row) + '\n')
	frt_handle.close()

	# close logging object and exit
	logObject.info('******************\nannotateProtSeqs.py finished!\n******************\nConsolidated TSV can be found at: %s' % annotation_tsv_file)
	util.closeLoggerObject(logObject)
	sys.exit(0)

if __name__ == '__main__':
	annotateProtSeqs()
