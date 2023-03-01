#!/usr/bin/env python

### Program: convertMiniprotGffToGbkAndProt.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

# BSD 3-Clause License
#
# Copyright (c) 2021, Kalan-Lab
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
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import defaultdict
import subprocess
from operator import itemgetter

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: convertMiniprotGffToGbkAndProt.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-g', '--miniprot_gff3', help='Path to miniprot produced GFF3.', required=True)
	parser.add_argument('-f', '--genome_fasta', help='Path to target genome used for miniprot.', required=True)
	parser.add_argument('-og', '--output_genbank', help='Path to output GenBank.', required=True)
	parser.add_argument('-op', '--output_proteome', help='Path to output Proteome (FASTA format).', required=True)
	parser.add_argument('-l', '--locus_tag', help='Locus tag. Default is AAA.', required=False, default='AAA')

	args = parser.parse_args()
	return args


def convertMiniProtGFFtoGenbank():
	myargs = create_parser()

	genome_fasta = myargs.genome_fasta
	miniprot_gff3 = myargs.miniprot_gff3
	output_genbank = myargs.output_genbank
	output_proteome = myargs.output_proteome
	locus_tag = myargs.locus_tag

	try:
		assert (os.path.isfile(genome_fasta) and os.path.isfile(miniprot_gff3))
	except:
		raise RuntimeError('Issue with validating inputs for miniprot_gff3 and/or genome_fasta are existing files.')

	# Step 1: Parse GFF3 for PAF CIGAR Strings
	query_mrna_paf_info = defaultdict(list)
	paf = None
	with open(miniprot_gff3) as omg:
		for line in omg:
			line = line.strip()
			ls = line.split('\t')
			if line.startswith('##PAF'):
				paf = ls[1:]
			elif len(ls) >= 6:
				if ls[2] == 'mRNA' and paf != None:
					query = line.split('Target=')[1].split(';')[0].split()[0]
					mrna_start_coord = int(ls[3])
					mrna_end_coord = int(ls[4])
					scaffold = ls[0]
					paf_start_coord = int(paf[7])
					paf_end_coord = int(paf[8])
					mrna_coords = set(range(mrna_start_coord, mrna_end_coord+1))
					paf_coords = set(range(paf_start_coord, paf_end_coord+1))
					if len(mrna_coords.intersection(paf_coords))/len(paf_coords) >= 0.95 and len(mrna_coords.intersection(paf_coords))/len(paf_coords) <= 1.05:
						paf_cigar = paf[17]
						paf_matches = re.findall(r'(\d+)([A-Z]{1})', paf_cigar)
						paf_cigar_parsed = [{'type': m[1], 'length': int(m[0])} for m in paf_matches]
						query_mrna_paf_info[tuple([query, mrna_start_coord])] = [scaffold, paf_start_coord, paf_end_coord, paf_cigar_parsed, '\t'.join(paf)]
				else:
					paf = None

	# Step 1: Parse GFF3 for transcript coordinates
	query_mrna_coords = defaultdict(list)
	scaffold_queries = defaultdict(list)
	with open(miniprot_gff3) as omg:
		for line in omg:
			line = line.strip()
			if line.startswith('#'): continue
			ls = line.split('\t')
			if ls[2] != 'mRNA': continue
			query = line.split('Target=')[1].split(';')[0].split()[0]
			scaffold = ls[0]
			start_coord = int(ls[3])
			end_coord = int(ls[4])
			score = float(ls[5])
			identity = float(ls[8].split('Identity=')[1].split(';')[0].split()[0])
			if identity < 0.80: continue
			dire = 1
			if ls[6] == '-':
				dire = -1
			# in case there are paralogs
			query_mrna_coords[tuple([query, start_coord])] = [scaffold, start_coord, end_coord, dire, score]
			scaffold_queries[scaffold].append([query, start_coord])

	# Step 2: Resolve overlap between mRNA transcripts
	redundant = set([])
	for i, mrna1 in enumerate(sorted(query_mrna_coords.items())):
		m1_coords = set(range(mrna1[1][1], mrna1[1][2]))
		m1_score = mrna1[1][4]
		for j, mrna2 in enumerate(sorted(query_mrna_coords.items())):
			if i >= j: continue
			# check scaffolds are the same
			if mrna1[1][0] != mrna2[1][0]: continue
			m2_coords = set(range(mrna2[1][1], mrna2[1][2]))
			if len(m1_coords.intersection(m2_coords)) > 0:
				m2_score = mrna2[1][4]
				if m1_score >= m2_score:
					redundant.add(mrna2[0])
				else:
					redundant.add(mrna1[0])

	# Step 3: Parse GFF3 for CDS coordinates
	query_cds_coords = defaultdict(list)
	with open(miniprot_gff3) as omg:
		for line in omg:
			line = line.strip()
			if line.startswith('#'): continue
			ls = line.split('\t')
			if ls[2] != 'CDS': continue
			query = line.split('Target=')[1].split(';')[0].split()[0]
			scaffold = ls[0]
			start_coord = int(ls[3])
			end_coord = int(ls[4])
			dire = 1
			if ls[6] == '-':
				dire = -1
			query_cds_coords[query].append([scaffold, start_coord, end_coord, dire])

	# Step 4: Go through FASTA scaffold/contig by scaffold/contig and create output GenBank
	gbk_handle = open(output_genbank, 'w')
	faa_handle = open(output_proteome, 'w')

	lt_iter = 0
	with open(genome_fasta) as ogf:
		for rec in SeqIO.parse(ogf, 'fasta'):
			seq = rec.seq
			gbk_rec = SeqRecord(seq, id=rec.id, name=rec.id, description=rec.description)
			gbk_rec.annotations['molecule_type'] = 'DNA'
			feature_list = []
			for mrna in sorted(scaffold_queries[rec.id], key=itemgetter(1)):
				qid, start = mrna
				key = tuple([qid, start])
				if key in redundant: continue
				mrna_info = query_mrna_coords[key]
				mrna_coords = set(range(mrna_info[1], mrna_info[2]))

				mrna_exon_locs = []
				all_coords = []
				for cds_info in query_cds_coords[qid]:
					cds_coords = set(range(cds_info[1], cds_info[2]))
					if len(mrna_coords.intersection(cds_coords)) > 0 and cds_info[0] == mrna_info[0]:
						all_coords.append([cds_info[1], cds_info[2]])
						assert(mrna_info[3] == cds_info[3])
						mrna_exon_locs.append(FeatureLocation(cds_info[1]-1, cds_info[2], strand=cds_info[3]))
				feature_loc = sum(mrna_exon_locs)
				feature = SeqFeature(feature_loc, type='CDS')
				lt = None
				if (lt_iter + 1) < 10:
					lt = '00000' + str(lt_iter + 1)
				elif (lt_iter + 1) < 100:
					lt = '0000' + str(lt_iter + 1)
				elif (lt_iter + 1) < 1000:
					lt = '000' + str(lt_iter + 1)
				elif (lt_iter + 1) < 10000:
					lt = '00' + str(lt_iter + 1)
				elif (lt_iter + 1) < 100000:
					lt = '0' + str(lt_iter + 1)
				else:
					lt = str(lt_iter + 1)
				lt_iter += 1

				"""
				nucl_seq = ''
				prot_seq = None
				for sc, ec in sorted(all_coords, key=itemgetter(0)):
					if ec >= len(seq):
						nucl_seq += seq[sc - 1:]
					else:
						nucl_seq += seq[sc - 1:ec]
				if mrna_info[3] == -1:
					prot_seq = str(Seq(nucl_seq).reverse_complement().translate())
				else:
					prot_seq = str(Seq(nucl_seq).translate())
				assert(prot_seq != None)
				"""

				scaffold, paf_start_coord, paf_end_coord, paf_cigar_parsed, paf_string = query_mrna_paf_info[key]

				paf_nucl_seq = ''
				paf_coord = paf_start_coord
				if mrna_info[3] == -1:
					paf_cigar_parsed.reverse()
				for op in paf_cigar_parsed:
					length = op['length']
					if op['type'] == 'M':
						paf_nucl_seq += seq[paf_coord:(paf_coord+(length*3))]
						paf_coord += length*3
					elif op['type'] == 'D':
						paf_coord += length*3
					elif op['type'] == 'G':
						#paf_coord += seq[paf_coord:paf_coord+length]
						paf_coord += length
					elif op['type'] in set(['F', 'N', 'U', 'V']):
						paf_coord += length

				paf_prot_seq = None
				paf_upstream_nucl_seq = None
				#paf_upstream_check = None
				if mrna_info[3] == -1:
					paf_nucl_seq = str(Seq(paf_nucl_seq).reverse_complement())
					paf_prot_seq = str(Seq(paf_nucl_seq).translate())
					paf_upstream_nucl_seq = str(Seq(seq[paf_end_coord:paf_end_coord+100]).reverse_complement())
					#paf_upstream_check = str(Seq(seq[paf_end_coord-1:paf_end_coord+99]).reverse_complement())
				else:
					paf_prot_seq = str(Seq(paf_nucl_seq).translate())
					paf_upstream_nucl_seq = seq[paf_start_coord-100:paf_start_coord]

				"""
				print('----------------------')
				print(paf_string)
				print(locus_tag + '_' + lt)
				print(paf_nucl_seq)
				print('++++++++++++++++++++++')
				print(paf_upstream_nucl_seq)
				print('++++++++++++++++++++++')
				print(paf_upstream_check)
				#print(paf_prot_seq)
				"""

				faa_handle.write('>' + locus_tag + '_' + lt + '\n' + paf_prot_seq + '\n')
				feature.qualifiers['locus_tag'] = locus_tag + '_' + lt
				feature.qualifiers['translation'] = paf_prot_seq
				feature.qualifiers['open_reading_frame'] = paf_nucl_seq
				feature.qualifiers['orf_upstream'] = paf_upstream_nucl_seq
				feature_list.append(feature)
			gbk_rec.features = feature_list
			SeqIO.write(gbk_rec, gbk_handle, 'genbank')
	gbk_handle.close()
	faa_handle.close()

if __name__ == '__main__':
	convertMiniProtGFFtoGenbank()
