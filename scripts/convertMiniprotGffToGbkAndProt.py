#!/usr/bin/env python3

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
from operator import itemgetter
import gzip

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: convertMiniprotGffToGbkAndProt.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-g', '--miniprot_gff3', help='Path to miniprot produced GFF3.', required=True)
	parser.add_argument('-f', '--genome_fasta', help='Path to target genome used for searching with miniprot.', required=True)
	parser.add_argument('-og', '--output_genbank', help='Path to output GenBank.', required=True)
	parser.add_argument('-op', '--output_proteome', help='Path to output Proteome (FASTA format).', required=True)
	parser.add_argument('-l', '--locus_tag', help='Locus tag [Default is AAA].', required=False, default='AAA')

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
						paf_cigar = None
						for item in paf:
							if item.startswith('cg:Z:'):
								paf_cigar = item.split('cg:Z:')[1]
						if paf_cigar == None: continue
						paf_matches = re.findall(r'(\d+)([A-Z]{1})', paf_cigar)
						paf_cigar_parsed = [{'type': m[1].strip(), 'length': int(m[0])} for m in paf_matches]
						query_mrna_paf_info[tuple([query, mrna_start_coord])] = [scaffold, paf_start_coord, paf_end_coord, paf_cigar_parsed, paf_cigar]
				else:
					paf = None

	# Step 1: Parse GFF3 for transcript coordinates
	query_mrna_coords = defaultdict(list)
	scaffold_queries = defaultdict(list)
	map_counter = 0
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
			query_mrna_coords[tuple([query, map_counter])] = [scaffold, start_coord, end_coord, dire, score]
			scaffold_queries[scaffold].append([query, start_coord, map_counter])
			map_counter += 1

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

	"""
	# legacy for previous handling in versions 1.5.2 and prior.

	# Step 3: Parse GFF3 for CDS coordinates
	query_cds_coords = defaultdict(list)
	query_cds_coords_accounted = defaultdict(set)
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
			score = float(ls[5])
			curr_coord = tuple([scaffold, start_coord, end_coord, dire, score])	
			if not curr_coord in query_cds_coords_accounted[query]:
				query_cds_coords[query].append([scaffold, start_coord, end_coord, dire, score])
				query_cds_coords_accounted[query].add(curr_coord)

	# handle overlapping exons
	cds_redundant = defaultdict(set)
	for query in sorted(query_cds_coords):
		for i, cdsc1 in enumerate(sorted(query_cds_coords[query], key=itemgetter(1))):
			c1_coords = set(range(cdsc1[1], cdsc1[2]))
			c1_score = cdsc1[4]
			for j, cdsc2 in enumerate(sorted(query_cds_coords[query], key=itemgetter(1))):
				if i >= j: continue
				# check scaffolds are the same
				if cdsc1[0] != cdsc2[0]: continue
				c2_coords = set(range(cdsc2[1], cdsc2[2]))
				if len(c1_coords.intersection(c2_coords)) > 0:
					c2_score = cdsc2[4]
					if c1_score >= c2_score:
						cds_redundant[query].add(tuple(cdsc2))
					else:
						cds_redundant[query].add(tuple(cdsc1))
	"""
	
	# Step 4: Go through FASTA scaffold/contig by scaffold/contig and create output GenBank
	gbk_handle = open(output_genbank, 'w')
	faa_handle = open(output_proteome, 'w')

	lt_iter = 0
	
	ogf = None
	if not genome_fasta.endswith('.gz'):
		ogf = open(genome_fasta)
	else:
		ogf = gzip.open(genome_fasta, 'rt')

	for rec in SeqIO.parse(ogf, 'fasta'):
		seq = str(rec.seq)
		gbk_rec = SeqRecord(seq, id=rec.id, name=rec.id, description=rec.description)
		gbk_rec.annotations['molecule_type'] = 'DNA'
		feature_list = []
		for mrna in sorted(scaffold_queries[rec.id], key=itemgetter(1)):
			qid, start, map_counter = mrna
			key_for_redundancy_check = tuple([qid, map_counter])
			key_for_paf = tuple([qid, start])
			if key_for_redundancy_check in redundant: continue
			mrna_info = query_mrna_coords[key_for_redundancy_check]
			mrna_start = mrna_info[1]
			mrna_end = mrna_info[2]
			mrna_coords = set(range(mrna_start, mrna_end))
			
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

			scaffold, paf_start_coord, paf_end_coord, paf_cigar_parsed, paf_string = query_mrna_paf_info[key_for_paf]
			
			paf_nucl_seq = ''
			paf_coord = paf_start_coord
			cds_positions = []
			if mrna_info[3] == -1:
				paf_cigar_parsed.reverse()
			for op in paf_cigar_parsed:
				length = op['length']
				if op['type'] in set(['M', 'D']):
					paf_nucl_seq += seq[paf_coord:(paf_coord+(length*3))]
					for pos in range(paf_coord, (paf_coord+(length*3))):
						cds_positions.append(pos)
					paf_coord += length*3
				elif op['type'] in set(['F', 'N']):
					paf_coord += length
				elif op['type'] == 'U':
					tmp = seq[paf_coord:(paf_coord+length)]
					if mrna_info[3] == -1:
						paf_nucl_seq += tmp[:2] + tmp[-1]
						cds_positions.append(paf_coord)
						cds_positions.append(paf_coord+1)
						cds_positions.append(paf_coord+length-1)
					else:
						paf_nucl_seq += tmp[0] + tmp[-2:]
						cds_positions.append(paf_coord)
						cds_positions.append(paf_coord+length-1)
						cds_positions.append(paf_coord+length-2)
					paf_coord += length
				elif op['type'] == 'V':
					tmp = seq[paf_coord:(paf_coord+length)]
					if mrna_info[3] == -1:
						paf_nucl_seq += tmp[0] + tmp[-2:]
						cds_positions.append(paf_coord)
						cds_positions.append(paf_coord+length-1)
						cds_positions.append(paf_coord+length-2)
					else:
						paf_nucl_seq += tmp[:2] + tmp[-1]
						cds_positions.append(paf_coord)
						cds_positions.append(paf_coord+1)
						cds_positions.append(paf_coord+length-1)
					paf_coord += length
				elif op['type'] == 'G':
					paf_coord += length

			cds_locs = []
			all_coords = []
			for cds_part in ranges(sorted(cds_positions)):
				cds_locs.append(FeatureLocation(cds_part[0], cds_part[1]+1, strand=mrna_info[3]))
				all_coords.append([cds_part[0]+1, cds_part[1]+1])

			feature_loc = sum(cds_locs)
			feature = SeqFeature(feature_loc, type='CDS')
			
			paf_prot_seq = None
			paf_upstream_nucl_seq = None
			if mrna_info[3] == -1:
				paf_nucl_seq = str(Seq(paf_nucl_seq).reverse_complement())
				paf_prot_seq = str(Seq(paf_nucl_seq).translate())
				paf_upstream_nucl_seq = str(Seq(seq[paf_end_coord:paf_end_coord+100]).reverse_complement())
			else:
				paf_prot_seq = str(Seq(paf_nucl_seq).translate())
				paf_upstream_nucl_seq = seq[paf_start_coord-100:paf_start_coord]

			cds_nucl_seq = ''
			cds_prot_seq = None
			for sc, ec in sorted(all_coords, key=itemgetter(0)):
				if ec >= len(seq):
					cds_nucl_seq += seq[sc - 1:]
				else:
					cds_nucl_seq += seq[sc - 1:ec]
			if mrna_info[3] == -1:
				cds_prot_seq = str(Seq(cds_nucl_seq).reverse_complement().translate())
			else:
				cds_prot_seq = str(Seq(cds_nucl_seq).translate())

			assert(paf_prot_seq == cds_prot_seq)

			faa_handle.write('>' + locus_tag + '_' + lt + '\n' + paf_prot_seq + '\n')
			feature.qualifiers['prot_from_ref'] = qid
			feature.qualifiers['locus_tag'] = locus_tag + '_' + lt
			feature.qualifiers['translation'] = paf_prot_seq
			feature.qualifiers['paf_nucl_seq'] = paf_nucl_seq
			#feature.qualifiers['cds_prot_seq'] = cds_prot_seq
			feature.qualifiers['cigar_string'] = paf_string
			#feature.qualifiers['full_orf'] = seq[paf_start_coord:paf_end_coord]
			feature.qualifiers['paf_upstream'] = paf_upstream_nucl_seq
			feature_list.append(feature)
		gbk_rec.features = feature_list
		SeqIO.write(gbk_rec, gbk_handle, 'genbank')

	ogf.close()
	gbk_handle.close()
	faa_handle.close()

def ranges(seq):
	"""
	Solution by Gareth Latty on StackOverflow:
	https://stackoverflow.com/questions/10420464/group-list-of-ints-by-continuous-sequence
	"""
	start, end = seq[0], seq[0]
	count = start
	for item in seq:
		if not count == item:
			yield start, end
			start, end = item, item
			count = item
		end = item
		count += 1
	yield start, end

if __name__ == '__main__':
	convertMiniProtGFFtoGenbank()
