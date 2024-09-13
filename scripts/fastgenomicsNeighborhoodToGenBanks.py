#!/usr/bin/env python3

### Program: fastgenomicsNeighborhoodToGenBanks.py
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
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import defaultdict
from operator import itemgetter
from zol import util
import traceback

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""               
	Program: fastgenomicsNeighborhoodToGenBanks.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

    fastgenomicsNeighborhoodToGenBanks.py: A program to download genomes for a neighborhood
	identified via fast.genomics (https://fast.genomics.lbl.gov/cgi/search.cgi) and extract
	neighborhood regions into individual GenBank files. These can then be provided into zol
	or as input for gene cluster visualization software such as clinker or pyGenomeViz.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--table_of_genes', help='The table_of_genes.tsv file downloaded from fast.genomics for the neighborhood of interest.', required=True)
	parser.add_argument('-o', '--output_dir', help='Path to output directory.', required=True)
	args = parser.parse_args()
	return args

def fastgenomicsNeighborhoodToGenBanks():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE INPUTS
	"""
	myargs = create_parser()

	table_of_genes_file = myargs.table_of_genes
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

	# Step 1: Concatenate all proteins into single multi-FASTA file
	msg = "--------------------\nStep 1\n--------------------\nProcessing fast.genomics table of genes TSV and getting list of genomic accessions."
	sys.stdout.write(msg + "\n")
	logObject.info(msg)
	
	gcf_accessions = set([])
	gca_accessions = set([])
	acc_scaff_cds_info = defaultdict(lambda: defaultdict(list))
	with open(table_of_genes_file) as otogf:
		for i, line in enumerate(otogf): 
			if i == 0: continue
			line = line.strip('\n')
			ls = line.split('\t')
			acc = ls[9]
			if acc.startswith('GCF_'):
				gcf_accessions.add(acc)
			elif acc.startswith('GCA_'):
				gca_accessions.add(acc)
			lt, pid, scaff, start, end, strand, desc = ls[10:17]
			cds_info = [lt, pid, scaff, int(start), int(end), strand, desc]
			acc_scaff_cds_info[acc][scaff].append(cds_info)

	# Step 2: Download GCAs and GCFs using ncbi-genome-download
	msg = "--------------------\nStep 2\n--------------------\nDownloading relevant genomes using ncbi-genome-download."
	sys.stdout.write(msg + "\n")
	logObject.info(msg)
	
	genome_fasta_dir = outdir + 'Downloaded_Genome_FASTAs/'
	os.mkdir(genome_fasta_dir)

	gca_accessions_file = outdir + 'GCA_Accessions.txt'
	gcf_accessions_file = outdir + 'GCF_Accessions.txt'

	gaah = open(gca_accessions_file, 'w')
	gfah = open(gcf_accessions_file, 'w')
	gaah.write('\n'.join(gca_accessions) + '\n')
	gfah.write('\n'.join(gcf_accessions) + '\n')
	gaah.close()
	gfah.close()

	ngd_gca_cmd = ['ncbi-genome-download', '--formats', 'fasta', '--retries', '2', '--section', 'genbank', '-A', gca_accessions_file, '-o', genome_fasta_dir, '--flat-output',  'all']	
	try:
		subprocess.call(' '.join(ngd_gca_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
	except:
		sys.stderr.write('Had an issue running: %s\n' % ' '.join(ngd_gca_cmd))
		sys.stderr.write(traceback.format_exc())
	
	ngd_gcf_cmd = ['ncbi-genome-download', '--formats', 'fasta', '--retries', '2', '--section', 'refseq', '-A', gcf_accessions_file, '-o', genome_fasta_dir, '--flat-output',  'all']	
	try:
		subprocess.call(' '.join(ngd_gcf_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
	except:
		sys.stderr.write('Had an issue running: %s\n' % ' '.join(ngd_gcf_cmd))
		sys.stderr.write(traceback.format_exc())

	os.system('gunzip ' + genome_fasta_dir + '*.gz')

	# Step 3: Create GenBank files
	msg = "--------------------\nStep 3\n--------------------\nCreating GenBank files for neighborhoods."
	sys.stdout.write(msg + "\n")
	logObject.info(msg)

	acc_to_genome_fasta = {}
	for f in os.listdir(genome_fasta_dir): 
		acc = '_'.join(f.split('_')[:2])
		acc_to_genome_fasta[acc] = genome_fasta_dir + f

	final_genbank_dir = outdir + 'Gene_Cluster_GenBank_Files/'
	os.mkdir(final_genbank_dir)

	for acc in acc_scaff_cds_info:
		acc_fasta_file = None 
		try:
			acc_fasta_file = acc_to_genome_fasta[acc]
			assert(acc_fasta_file)
		except:
			msg = 'Could not download genome for accession: %s, skipping it!' % acc
			logObject.warning(msg)
			sys.stderr.write(msg + '\n')
			continue

		scaff_seqs = {}
		try:
			with open(acc_fasta_file) as oaff:
				for rec in SeqIO.parse(oaff, 'fasta'):
					scaff_seqs[rec.id] = str(rec.seq)
		except:
			msg = 'Issues with reading genome file: %s' % acc_fasta_file
			logObject.warning(msg)
			sys.stderr.write(msg + '\n')

		for scaff in acc_scaff_cds_info[acc]:
			all_gene_bounds = []	
			for cds in acc_scaff_cds_info[acc][scaff]:
				all_gene_bounds.append(cds[3])
				all_gene_bounds.append(cds[4])

			region_start = min(all_gene_bounds)
			region_end = max(all_gene_bounds)
			scaff_seq = scaff_seqs[scaff] 

			region_seq = scaff_seq[region_start-1:region_end]
			
			gbk_file = final_genbank_dir + acc + '_scaffold_' + scaff + '.gbk'
			gbk_handle = open(gbk_file, 'w')
			record = SeqRecord(Seq(region_seq), id=scaff, name=scaff, description='')
			record.annotations['molecule_type'] = 'DNA'
			feature_list = []
			for cds in acc_scaff_cds_info[acc][scaff]:
				lt, pid, _, gene_start, gene_end, gene_strand, desc = cds
				rel_start = gene_start - region_start + 1 
				rel_end = gene_end - region_start + 1
				gene_seq = Seq(scaff_seq[gene_start-1:gene_end])
				gene_strand_num = 1
				if gene_strand == '-':
					gene_strand_num = -1
					gene_seq = str(Seq(scaff_seq[gene_start:gene_end]).reverse_complement())
				prot_seq = str(Seq(gene_seq).translate())
				feature = SeqFeature(FeatureLocation(start=rel_start-1, end=rel_end, strand=gene_strand_num), type='CDS')
				feature.qualifiers['locus_tag'] = lt
				feature.qualifiers['product'] = desc
				feature.qualifiers['protein_id'] = pid
				feature.qualifiers['translation'] = prot_seq.rstrip('*')
				feature_list.append(feature)
			record.features = feature_list
			SeqIO.write(record, gbk_handle, 'genbank')
			gbk_handle.close()

	# DONE!
	sys.stdout.write("--------------------\nDONE!\n--------------------\nDirectory of GenBank files for gene neighborhood can be found at: %s\n" % final_genbank_dir)
	logObject.info("--------------------\nDONE!\n--------------------\nDirectory of GenBank files for gene neighborhood can be found at: %s" % final_genbank_dir)

if __name__ == '__main__':
	fastgenomicsNeighborhoodToGenBanks()

