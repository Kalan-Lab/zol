import os
import sys
from zol import util, zol
from Bio import SeqIO
from collections import defaultdict
import multiprocessing
import subprocess
import decimal
from operator import itemgetter
import gzip
from pomegranate import *
import numpy
import traceback
from scipy.stats import pearsonr
import shutil
import pickle

zol_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
plot_prog = zol_main_directory + 'zol/plotSegments.R'
split_diamond_results_prog = os.path.abspath(os.path.dirname(__file__) + '/') + '/splitDiamondResultsForFai'

def subsetGenBankForQueryLocus(full_gw_genbank, locus_genbank, locus_proteins, reference_contig, reference_start, reference_end, logObject):
	try:
		util.createGenbank(full_gw_genbank, locus_genbank, reference_contig, reference_start, reference_end)

		locus_proteins_handle = open(locus_proteins, 'w')
		with open(locus_genbank) as olg:
			for rec in SeqIO.parse(olg, 'genbank'):
				for feature in rec.features:
					if feature.type != 'CDS': continue
					lt = feature.qualifiers.get('locus_tag')[0]
					prot_seq = feature.qualifiers.get('translation')[0]
					locus_proteins_handle.write('>' + lt + '\n' + str(prot_seq) + '\n')
		locus_proteins_handle.close()

	except Exception as e:
		sys.stderr.write('Issues subsetting locus-specific GenBank from genome-wide GenBank.\n')
		logObject.error('Issues subsetting locus-specific GenBank from genome-wide GenBank.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def parseHomologGroupMatrix(orthogroup, logObject):
	protein_to_hg = {}
	try:
		with open(orthogroup) as oog:
			for i, line in enumerate(oog):
				line = line.strip()
				ls = line.split('\t')
				if i == 0: continue
				hg = ls[0]
				for lts in ls[1:]:
					for lt in lts.split(', '):
						protein_to_hg[lt] = hg
	except Exception as e:
		sys.stderr.write('Issues mapping proteins to homolog groups.\n')
		logObject.error('Issues mapping proteins to homolog groups.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)
	return protein_to_hg

def collapseProteinsUsingCDHit(protein_fasta, nr_protein_fasta, logObject):
	protein_to_hg = {}
	try:
		cdhit_cluster_file = nr_protein_fasta + '.clstr'
		cdhit_cmd = ['cd-hit', '-i', protein_fasta, '-o', nr_protein_fasta, '-c', '0.95', '-aL', '0.90', '-aS', '0.90',
					 '-d', '0']

		try:
			subprocess.call(' '.join(cdhit_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL, executable='/bin/bash')
			assert(os.path.isfile(nr_protein_fasta))
			assert(os.path.isfile(cdhit_cluster_file))
		except Exception as e:
			logObject.error("Issue with running: %s" % ' '.join(cdhit_cmd))
			logObject.error(e)
			raise RuntimeError(e)

		prefix = '.'.join(protein_fasta.split('/')[-1].split('.')[:-1])
		cluster_rep = None
		tmp = []
		with open(cdhit_cluster_file) as occf:
			for line in occf:
				line = line.strip()
				if line.startswith('>'):
					if len(tmp) > 0 and cluster_rep != None:
						for plt in tmp:
							protein_to_hg[plt] = cluster_rep
					tmp = []
					cluster_rep = None
					continue
				ls = line.split()
				lt = ls[2][1:-3]
				if line.endswith(' *'):
					cluster_rep = lt
				tmp.append(prefix + '|' + lt)
		if len(tmp) > 0 and cluster_rep != None:
			for plt in tmp:
				protein_to_hg[plt] = cluster_rep
	except Exception as e:
		sys.stderr.write('Issues running CD-hit.\n')
		logObject.error('Issues running CD-hit.\n')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)
	return protein_to_hg

def createDummyProteinToHomologGroupMapping(protein_fasta, logObject):
	protein_to_hg = {}
	try:
		with open(protein_fasta, 'fasta') as opf:
			for rec in SeqIO.parse(opf, 'fasta'):
				protein_to_hg[rec.id] = rec.id
	except Exception as e:
		sys.stderr.write('Issues creating fake mapping between protein IDs to themselves.\n')
		logObject.error('Issues creating fake mapping between protein IDs to themselves.\n')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)
	return protein_to_hg

def parseCoordsFromGenbank(genbanks, logObject):
	""" Simple function to parse """
	# get coords for individual genes
	comp_gene_info = defaultdict(dict)
	try:
		for gbk in genbanks:
			prefix = '.'.join(gbk.split('/')[-1].split('.')[:-1])
			gene_locations = util.parseGbk(gbk, prefix, logObject)
			for g in gene_locations:
				comp_gene_info[prefix][g] = gene_locations[g]
	except Exception as e:
		sys.stderr.write('Issues with parsing coordinates of genes in GenBanks!\n')
		logObject.error('Issues with parsing coordinates of genes in GenBanks!')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)
	return dict(comp_gene_info)

def genConsensusSequences(genbanks, outdir, logObject, cpus=1, use_super5=False):
		# determine orthologs
		prot_dir = outdir + 'CDS_Protein/'
		nucl_dir = outdir + 'CDS_Nucleotide/'
		og_dir = outdir + 'Determine_Orthogroups/'
		ortho_matrix_file = og_dir + 'Orthogroups.tsv'
		util.setupReadyDirectory([prot_dir, nucl_dir])
		try:
			for gbk in genbanks:
				try:
					prefix = '.'.join(gbk.split('/')[-1].split('.')[:-1])
					proteins, nucleotides, upstream_regions = util.parseGenbankForCDSProteinsAndDNA(gbk, logObject)
					protein_outf = prot_dir + prefix + '.faa'
					protein_handle = open(protein_outf, 'w')
					for lt in proteins:
						protein_handle.write('>' + prefix + '|' + lt + '\n' + proteins[lt] + '\n')
					protein_handle.close()
				except Exception as e:
					sys.stderr.write('Issues with parsing the GenBank %s\n' % gbk)
					logObject.error('Issues with parsing the GenBank %s' % gbk)
					sys.stderr.write(str(e) + '\n')
					sys.exit(1)
			fo_cmd = ['findOrthologs.py', '-p', prot_dir, '-o', og_dir, '-c', str(cpus)]
			try:
				subprocess.call(' '.join(fo_cmd), shell=True, stdout=subprocess.DEVNULL,
								stderr=subprocess.DEVNULL, executable='/bin/bash')
				assert (os.path.isfile(ortho_matrix_file))
				assert (util.checkCoreHomologGroupsExist(ortho_matrix_file))
			except Exception as e:
				logObject.error("Issue with running: %s" % ' '.join(fo_cmd))
				logObject.error(e)
				raise RuntimeError(e)

		except Exception as e:
			sys.stderr.write('Issues with determining ortholog/homolog groups!\n')
			logObject.error('Issues with determining ortholog/homolog groups!')
			sys.stderr.write(str(e) + '\n')
			sys.exit(1)

		# create Alignments, phylogenies and consensus sequences
		proc_dir = outdir + 'Homolog_Group_Processing/'
		hg_prot_dir = proc_dir + 'HG_Protein_Sequences/'
		hg_nucl_dir = proc_dir + 'HG_Nucleotide_Sequences/'
		prot_algn_dir = proc_dir + 'HG_Protein_Alignments/'
		phmm_dir = proc_dir + 'HG_Profile_HMMs/'
		cons_dir = proc_dir + 'HG_Consensus_Sequences/'
		consensus_prot_seqs_faa = outdir + 'HG_Consensus_Seqs.faa'
		util.setupReadyDirectory([proc_dir, prot_algn_dir, phmm_dir, cons_dir, hg_prot_dir])
		zol.partitionSequencesByHomologGroups(ortho_matrix_file, prot_dir, nucl_dir, hg_prot_dir, hg_nucl_dir, logObject)
		zol.createProteinAlignments(hg_prot_dir, prot_algn_dir, logObject, use_super5=use_super5, cpus=cpus)
		zol.createProfileHMMsAndConsensusSeqs(prot_algn_dir, phmm_dir, cons_dir, logObject, cpus=1)
		consensus_prot_seqs_handle = open(consensus_prot_seqs_faa, 'w')
		for f in os.listdir(cons_dir):
			with open(cons_dir + f) as ocf:
				for rec in SeqIO.parse(ocf, 'fasta'):
					consensus_prot_seqs_handle.write('>' + f.split('.cons.faa')[0] + '\n' + str(rec.seq) + '\n')
		consensus_prot_seqs_handle.close()
		return([ortho_matrix_file, consensus_prot_seqs_faa])

def loadTargetGenomeInfo(target_annotation_information, target_genomes_pkl_dir, diamond_reuslts, valid_tg_samples, logObject, min_genes_per_scaffold=2, lowmem_mode=True):
	gene_locations = {}
	scaffold_genes = {}
	boundary_genes = {}
	gene_id_to_order = {}
	gene_order_to_id = {}

	genes_hit_in_diamond_search = set(diamond_reuslts.keys())
	try:
		for sample in target_annotation_information:
			if not sample in valid_tg_samples: continue
			sample_info_pkl_file = target_genomes_pkl_dir + sample + '.pkl'
			genbank_info = None
			try:
				assert(os.path.isfile(sample_info_pkl_file))
				with open(sample_info_pkl_file, 'rb') as handle:
					genbank_info = pickle.load(handle)
			except:
				sys.stderr.write('Could not find pickle file with CDS coordinate information for sample %s\n' % sample)
				logObject.error('Could not find pickle file with CDS coordinate information for sample %s' % sample)
				sys.exit(1)

			gene_locs = genbank_info[0]
			scaff_genes = genbank_info[1]
			bound_genes = genbank_info[2]
			gito = genbank_info[3]
			goti = genbank_info[4]

			if lowmem_mode:
				valid_scaffolds = set([])
				gene_locs_filt = {}
				scaff_genes_filt = defaultdict(set)
				bound_genes_filt = set([])
				gito_filt = defaultdict(dict)
				goti_filt = defaultdict(dict)
				for scaffold in scaff_genes:
					if len(scaff_genes[scaffold].intersection(genes_hit_in_diamond_search)) >= min_genes_per_scaffold:
						valid_scaffolds.add(scaffold)
						scaff_genes_filt[scaffold] = scaff_genes[scaffold]
						gito_filt = gito[scaffold]
						goti_filt = goti[scaffold]
						for lt in scaff_genes[scaffold]:
							gene_locs_filt[lt] = gene_locs[lt]
						bound_genes_filt = bound_genes_filt.union(scaff_genes[scaffold].intersection(bound_genes))
				gene_locations[sample] = gene_locs_filt
				scaffold_genes[sample] = scaff_genes_filt
				boundary_genes[sample] = bound_genes_filt
				gene_id_to_order[sample] = gito_filt
				gene_order_to_id[sample] = goti_filt
			else:
				gene_locations[sample] = gene_locs
				scaffold_genes[sample] = scaff_genes
				boundary_genes[sample] = bound_genes
				gene_id_to_order[sample] = gito
				gene_order_to_id[sample] = goti

	except Exception as e:
		logObject.error('Issues with parsing CDS location information from genes of target genomes.')
		sys.stderr.write('Issues with parsing CDS location information from genes of target genomes.\n')
		sys.stderr.write(str(e) + '\n')
		raise RuntimeError(traceback.format_exc())
		sys.exit(1)


	target_genome_gene_info = {'gene_locations': gene_locations, 'scaffold_genes': scaffold_genes,
							   'boundary_genes': boundary_genes, 'gene_id_to_order': gene_id_to_order,
							   'gene_order_to_id': gene_order_to_id}
	return(target_genome_gene_info)

def runDiamondBlastp(target_concat_genome_db, query_fasta, diamond_work_dir, logObject, diamond_sensitivity='very-sensitive', evalue_cutoff=1e-10, cpus=1):
	diamond_results_file = diamond_work_dir + 'DIAMOND_Results.txt'
	try:
		diamond_blastp_cmd = ['diamond', 'blastp', '--ignore-warnings', '--threads', str(cpus), '--' + diamond_sensitivity,
					   '--query', query_fasta, '--db', target_concat_genome_db, '--outfmt', '6', 'qseqid', 'sseqid',
					   'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
					   'bitscore', 'qlen', 'slen', '-k0', '--out', diamond_results_file, '--evalue', str(evalue_cutoff)]
		try:
			subprocess.call(' '.join(diamond_blastp_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			logObject.info('Successfully ran: %s' % ' '.join(diamond_blastp_cmd))
		except Exception as e:
			logObject.error('Had an issue running: %s' % ' '.join(diamond_blastp_cmd))
			sys.stderr.write('Had an issue running: %s' % ' '.join(diamond_blastp_cmd))
			logObject.error(e)
			sys.exit(1)
		assert(os.path.isfile(diamond_results_file) and os.path.getsize(diamond_results_file) > 0)

	except Exception as e:
		logObject.error('Issues with running DIAMOND blastp.')
		sys.stderr.write('Issues with running DIAMOND blastp.\n')
		sys.stderr.write(str(e) + '\n')
		raise RuntimeError(traceback.format_exc())
		sys.exit(1)
	return(diamond_results_file)

def processDiamondBlastp(target_annot_information, diamond_results_file, work_dir, logObject):
	"""
	Function to process DIAMOND results
	"""
	diamond_results = defaultdict(list)
	try:
		search_res_dir = work_dir + 'Alignment_Results/'
		util.setupReadyDirectory([search_res_dir])

		split_mapping_file = work_dir + 'Sample_to_Listing.txt'
		split_mapping_handle = open(split_mapping_file, 'w')
		for sample in target_annot_information:
			sample_result_file = search_res_dir + sample + '.txt'
			split_mapping_handle.write(sample + '\t' + sample_result_file + '\n')
		split_mapping_handle.close()

		split_diamond_cmd = [split_diamond_results_prog, diamond_results_file, split_mapping_file]

		try:
			subprocess.call(' '.join(split_diamond_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			logObject.info('Successfully ran: %s' % ' '.join(split_diamond_cmd))
		except Exception as e:
			logObject.error('Had an issue running: %s' % ' '.join(split_diamond_cmd))
			sys.stderr.write('Had an issue running: %s' % ' '.join(split_diamond_cmd))
			logObject.error(e)
			sys.exit(1)

		for sample in target_annot_information:
			result_file = search_res_dir + sample + '.txt'
			try:
				assert (os.path.isfile(result_file))
				best_hit_per_lt = defaultdict(lambda: [[], 0.0, [], [], []])
				with open(result_file) as orf:
					for line in orf:
						line = line.strip()
						ls = line.split()
						hg = ls[0]
						lt = ls[1].split('|')[1]
						identity = float(ls[2])
						eval = decimal.Decimal(ls[10])
						bitscore = float(ls[11])
						qlen = int(ls[12])
						slen = int(ls[13])
						sql_ratio = float(slen)/float(qlen)
						if bitscore > best_hit_per_lt[lt][1]:
							best_hit_per_lt[lt][0] = [hg]
							best_hit_per_lt[lt][1] = bitscore
							best_hit_per_lt[lt][2] = [eval]
							best_hit_per_lt[lt][3] = [identity]
							best_hit_per_lt[lt][4] = [sql_ratio]
						elif bitscore == best_hit_per_lt[lt][1]:
							best_hit_per_lt[lt][0].append(hg)
							best_hit_per_lt[lt][2].append(eval)
							best_hit_per_lt[lt][3].append(identity)
							best_hit_per_lt[lt][4].append(sql_ratio)

				for lt in best_hit_per_lt:
					for i, hg in enumerate(best_hit_per_lt[lt][0]):
						diamond_results[lt].append([hg, best_hit_per_lt[lt][1], best_hit_per_lt[lt][2][i], sample,
													best_hit_per_lt[lt][3][i], best_hit_per_lt[lt][4][i]])
			except:
				#raise RuntimeError(traceback.format_exc())
				sys.stderr.write('Warning: Did not detect homology whatsoever at requested e-value for sample %s\'s proteome.\n' % sample)
				logObject.warning('Did not detect homology whatsoever at requested e-value for sample %s\'s proteome.' % sample)

	except Exception as e:
		logObject.error('Issues with running DIAMOND blastp or processing of results.')
		sys.stderr.write('Issues with running DIAMOND blastp or processing of results.\n')
		sys.stderr.write(str(e) + '\n')
		raise RuntimeError(traceback.format_exc())
		sys.exit(1)

	return(dict(diamond_results))

def parseGenbankAndFindBoundaryGenes(inputs):
	"""
	Function to parse Genbanks from Prokka and return a dictionary of genes per scaffold, gene to scaffold, and a
	set of genes which lie on the boundary of scaffolds.
	:param sample_genbank: Prokka generated Genbank file.
	:param distance_to_scaffold_boundary: Distance to scaffold edge considered as boundary.
	:return gene_to_scaffold: Dictionary mapping each gene's locus tag to the scaffold it is found on.
	:return scaffold_genes: Dictionary with keys as scaffolds and values as a set of genes found on that scaffold.
	:return boundary_genes: Set of gene locus tag ids which are found within proximity to scaffold edges.
	"""

	distance_to_scaffold_boundary = 2000
	gene_location = {}
	scaffold_genes = defaultdict(set)
	boundary_genes = set([])
	gene_id_to_order = defaultdict(dict)
	gene_order_to_id = defaultdict(dict)

	sample, sample_genbank, sample_gbk_info = inputs
	osg = None
	if sample_genbank.endswith('.gz'):
		osg = gzip.open(sample_genbank, 'rt')
	else:
		osg = open(sample_genbank)
	for rec in SeqIO.parse(osg, 'genbank'):
		scaffold = rec.id
		scaffold_length = len(str(rec.seq))
		boundary_ranges = set(range(1, distance_to_scaffold_boundary + 1)).union(
			set(range(scaffold_length - distance_to_scaffold_boundary, scaffold_length + 1)))
		gene_starts = []
		for feature in rec.features:
			if not feature.type == 'CDS': continue
			locus_tag = feature.qualifiers.get('locus_tag')[0]

			start = None
			end = None
			direction = None
			if not 'join' in str(feature.location):
				start = min(
					[int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
				end = max([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
				direction = str(feature.location).split('(')[1].split(')')[0]
			else:
				all_starts = []
				all_ends = []
				all_directions = []
				for exon_coord in str(feature.location)[5:-1].split(', '):
					start = min([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
					end = max([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
					direction = exon_coord.split('(')[1].split(')')[0]
					all_starts.append(start)
					all_ends.append(end)
					all_directions.append(direction)
				start = min(all_starts)
				end = max(all_ends)
				direction = all_directions[0]

			gene_location[locus_tag] = {'scaffold': scaffold, 'start': start, 'end': end, 'direction': direction}
			scaffold_genes[scaffold].add(locus_tag)

			gene_range = set(range(start, end + 1))
			if len(gene_range.intersection(boundary_ranges)) > 0:
				boundary_genes.add(locus_tag)

			gene_starts.append([locus_tag, start])

		for i, g in enumerate(sorted(gene_starts, key=itemgetter(1))):
			gene_id_to_order[scaffold][g[0]] = i
			gene_order_to_id[scaffold][i] = g[0]
	osg.close()
	sample_gbk_info[sample] = [gene_location, dict(scaffold_genes), boundary_genes, dict(gene_id_to_order),
							   dict(gene_order_to_id)]


def mapKeyProteinsToHomologGroups(query_fasta, key_protein_queries_fasta, work_dir, logObject, cpus=1):
	"""
	Function to align key protein queries fasta to query hg fasta file and determine set of hgs which are key.
	"""
	key_hgs = set([])
	try:
		search_res_dir = work_dir + 'Alignment_Results/'
		util.setupReadyDirectory([search_res_dir])

		query_dmnd_db = query_fasta.split('.faa')[0] + '.dmnd'
		align_result_file = work_dir + 'Key_Proteins_to_General_Queries_Alignment.txt'
		diamond_cmd = ['diamond', 'makedb', '--ignore-warnings', '--in', query_fasta, '-d', query_dmnd_db, ';',
					   'diamond', 'blastp', '--ignore-warnings', '--threads', str(cpus), '--very-sensitive', '--query',
					   key_protein_queries_fasta, '--db', query_dmnd_db, '--outfmt', '6', '--out',
					   align_result_file, '--evalue', '1e-10']

		try:
			subprocess.call(' '.join(diamond_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL, executable='/bin/bash')
			assert(os.path.isfile(align_result_file))
		except Exception as e:
			logObject.error("Issue with running: %s" % ' '.join(diamond_cmd))
			logObject.error(e)
			raise RuntimeError(e)

		kq_top_hits = defaultdict(lambda: [[], 0.0])
		with open(align_result_file) as oarf:
			for line in oarf:
				line = line.strip()
				ls = line.split('\t')
				kq, hg = ls[:2]
				bs = float(ls[11])
				if bs > kq_top_hits[kq][1]:
					kq_top_hits[kq] = [[hg], bs]
				elif bs == kq_top_hits[kq][1]:
					kq_top_hits[kq][0].append(hg)

		with open(key_protein_queries_fasta) as okpqf:
			for rec in SeqIO.parse(okpqf, 'fasta'):
				try:
					assert(len(kq_top_hits[rec.id][0]) >= 1)
					for hg in kq_top_hits[rec.id][0]:
						key_hgs.add(hg)
				except Exception as e2:
					logObject.error('Found no mapping for key query protein %s to general proteins.' % rec.id)
					sys.stderr.write('Found no mapping for key query protein %s to general proteins.\n' % rec.id)
					sys.stderr.write(str(e2) + '\n')
					sys.exit(1)

	except Exception as e:
		logObject.error('Issues with mapping key proteins to general proteins or consensus sequence of homolog groups.')
		sys.stderr.write('Issues with mapping key proteins to general proteins or consensus sequence of homolog groups.\n')
		sys.stderr.write(str(e) + '\n')
		raise RuntimeError(traceback.format_exc())
		sys.exit(1)
	return key_hgs


gc_genbanks_dir, gc_info_dir, query_gene_info, lt_to_hg, model, target_annotation_info, boundary_genes = [None]*7
sample_lt_to_evalue, sample_lt_to_identity, sample_lt_to_sqlratio, sample_lt_to_bitscore = [None]*4
lts_ordered_dict, hgs_ordered_dict, gene_locations, gene_id_to_order, gene_order_to_id = [None]*5

def identifyGCInstances(query_information, target_information, diamond_results, work_dir, logObject, min_hits=5,
						 min_key_hits=3, draft_mode=False, gc_to_gc_transition_prob=0.9, bg_to_bg_transition_prob=0.9,
						 gc_emission_prob_with_hit=0.95,  gc_emission_prob_without_hit=0.2,
						 syntenic_correlation_threshold=0.8, max_int_genes_for_merge=0,  kq_evalue_threshold=1e-20,
						 flanking_context=1000, cpus=1, block_size=3000, gc_delineation_mode="GENE-CLUMPER"):
	"""
	Function to search for instances of Gene Cluster in samples using HMM based approach based on homolog groups as
	characters. This function utilizes the convenient Python library Pomegranate.

	:param query_information: Dictionary with 3 keys: protein_to_hg (dict), query_fasta (file path), comp_gene_info (dict)
	:param target_information: Dictionary with 2 keys: target_annotation_information (dict), target_genome_gene_info (dict)
	:param diamond_results: Dictionary of DIAMOND alignment results meeting e-value threshold
	"""

	try:
		global gc_genbanks_dir, gc_info_dir, query_gene_info, lt_to_hg, model, target_annotation_info, boundary_genes
		global sample_lt_to_evalue, sample_lt_to_identity, sample_lt_to_sqlratio, sample_lt_to_bitscore
		global lts_ordered_dict, hgs_ordered_dict, gene_locations, gene_id_to_order, gene_order_to_id

		gc_genbanks_dir = os.path.abspath(work_dir + 'GeneCluster_Genbanks') + '/'
		gc_info_dir = os.path.abspath(work_dir + 'GeneCluster_Info') + '/'
		util.setupReadyDirectory([gc_genbanks_dir, gc_info_dir])

		# unpack information in dictionaries

		query_gene_info = query_information['comp_gene_info']
		lt_to_hg = query_information['protein_to_hg']
		all_hgs = set(lt_to_hg.values())
		key_hgs = query_information['key_hgs']

		target_annotation_info = target_information['target_annotation_information']
		target_genome_gene_info = target_information['target_genome_gene_info']

		gene_locations = target_genome_gene_info['gene_locations']
		scaffold_genes = target_genome_gene_info['scaffold_genes']
		boundary_genes = target_genome_gene_info['boundary_genes']
		gene_id_to_order = target_genome_gene_info['gene_id_to_order']
		gene_order_to_id = target_genome_gene_info['gene_order_to_id']

		# Create HMM
		gc_hg_probabilities = {'background': gc_emission_prob_without_hit}
		bg_hg_probabilities = {'background': 1.0-gc_emission_prob_without_hit}
		for hg in all_hgs:
			gc_hg_probabilities[hg] = gc_emission_prob_with_hit
			bg_hg_probabilities[hg] = 1.0-gc_emission_prob_with_hit

		gc_distribution = DiscreteDistribution(dict(gc_hg_probabilities))
		bg_distribution = DiscreteDistribution(dict(bg_hg_probabilities))

		gc_state = State(gc_distribution, name='Gene-Cluster')
		bg_state = State(bg_distribution, name='Other')

		model = HiddenMarkovModel()
		model.add_states(gc_state, bg_state)

		gc_to_gc = gc_to_gc_transition_prob
		gc_to_bg = 1.0 - gc_to_gc_transition_prob
		bg_to_bg = bg_to_bg_transition_prob
		bg_to_gc = 1.0 - bg_to_bg_transition_prob

		start_to_gc = 0.5
		start_to_bg = 0.5
		gc_to_end = 0.5
		bg_to_end = 0.5

		model.add_transition(model.start, gc_state, start_to_gc)
		model.add_transition(model.start, bg_state, start_to_bg)
		model.add_transition(gc_state, model.end, gc_to_end)
		model.add_transition(bg_state, model.end, bg_to_end)
		model.add_transition(gc_state, gc_state, gc_to_gc)
		model.add_transition(gc_state, bg_state, gc_to_bg)
		model.add_transition(bg_state, gc_state, bg_to_gc)
		model.add_transition(bg_state, bg_state, bg_to_bg)

		model.bake()

		gc_hmm_evalues_file = work_dir + 'GeneCluster_NewInstances_HMMEvalues.txt'
		gc_list_file = work_dir + 'GeneCluster_Homolog_Listings.txt'

		# open handle to file where expanded GCF listings will be written
		total_samples = sorted(target_annotation_info.keys())
		for block_samp_list in util.chunks(total_samples, block_size):
			block_samp_set = set(block_samp_list)

			# start process of finding new BGCs
			sample_lt_to_hg = defaultdict(dict)
			sample_hgs = defaultdict(set)
			sample_lt_to_evalue = defaultdict(dict)
			sample_lt_to_identity = defaultdict(dict)
			sample_lt_to_sqlratio = defaultdict(dict)
			sample_lt_to_bitscore = defaultdict(dict)
			for lt in diamond_results:
				for i, hits in enumerate(sorted(diamond_results[lt], key=itemgetter(1))):
					if i == 0 and hits[3] in block_samp_set: # for now just selecting one - but really it should never be the case that this array is greater than one
						sample_lt_to_hg[hits[3]][lt] = hits[0]
						sample_hgs[hits[3]].add(hits[0])
						sample_lt_to_bitscore[hits[3]][lt] = hits[1]
						sample_lt_to_evalue[hits[3]][lt] = decimal.Decimal(hits[2])
						sample_lt_to_identity[hits[3]][lt] = float(hits[4])
						sample_lt_to_sqlratio[hits[3]][lt] = float(hits[5])
			identify_gc_segments_input = []
			hgs_ordered_dict = defaultdict(dict)
			lts_ordered_dict = defaultdict(dict)
			for sample in sample_hgs:
				if len(sample_hgs[sample]) < 3: continue
				for scaffold in scaffold_genes[sample]:
					lts_with_start = []
					for lt in scaffold_genes[sample][scaffold]:
						lts_with_start.append([lt, gene_locations[sample][lt]['start']])

					hgs_ordered = []
					lts_ordered = []
					for lt_start in sorted(lts_with_start, key=itemgetter(1)):
						lt, start = lt_start
						lts_ordered.append(lt)
						if lt in sample_lt_to_hg[sample]:
							hgs_ordered.append(sample_lt_to_hg[sample][lt])
						else:
							hgs_ordered.append('background')
					hgs_ordered_dict[sample][scaffold] = hgs_ordered
					lts_ordered_dict[sample][scaffold] = lts_ordered

				identify_gc_segments_input.append([sample, min_hits, min_key_hits, key_hgs, kq_evalue_threshold,
												   syntenic_correlation_threshold, max_int_genes_for_merge,
												   flanking_context, draft_mode, gc_delineation_mode])

			p = multiprocessing.Pool(cpus)
			p.map(identify_gc_instances, identify_gc_segments_input)
			p.close()

		os.system('find %s -type f -name "*.bgcs.txt" -exec cat {} + >> %s' % (gc_info_dir, gc_list_file))
		os.system('find %s -type f -name "*.hg_evalues.txt" -exec cat {} + >> %s' % (gc_info_dir, gc_hmm_evalues_file))

	except Exception as e:
		logObject.error('Issues with managing running of HMM to find gene-cluster homolog segments.')
		sys.stderr.write('Issues with managing running of HMM to find gene-cluster homolog segments.\n')
		sys.stderr.write(str(e) + '\n')
		raise RuntimeError(traceback.format_exc())
		sys.exit(1)

def identify_gc_instances(input_args):
	sample, min_hits, min_key_hits, key_hgs, kq_evalue_threshold, syntenic_correlation_threshold, max_int_genes_for_merge, flanking_context, draft_mode, gc_delineation_mode = input_args
	sample_gc_predictions = []
	if gc_delineation_mode == 'GENE-CLUMPER':
		gcs_id = 1
		for scaffold in hgs_ordered_dict[sample]:
			hgs_ordered = hgs_ordered_dict[sample][scaffold]
			lts_ordered = lts_ordered_dict[sample][scaffold]
			tmp = []
			dist_counter = 0
			hg_counter = 0
			last_hg_i = 0
			active_search_mode = False
			for i, hg in enumerate(hgs_ordered):
				if hg != 'background':
					tmp.append(i)
					dist_counter = 0
					hg_counter += 1
					last_hg_i = i
					active_search_mode = True
				elif dist_counter < max_int_genes_for_merge and active_search_mode:
					tmp.append(i)
					dist_counter += 1
				else:
					if hg_counter >= 3:
						gc_state_lts = []
						gc_state_hgs = []

						begin_i = min(tmp)
						if last_hg_i == (len(hgs_ordered)-1):
							gc_state_lts = lts_ordered[begin_i:]
							gc_state_hgs = hgs_ordered[begin_i:]
						else:
							gc_state_lts = lts_ordered[min(tmp):last_hg_i+1]
							gc_state_hgs = hgs_ordered[min(tmp):last_hg_i+1]
						boundary_lt_featured = False
						features_key_hg = False
						if len(key_hgs.intersection(set(gc_state_hgs).difference('background'))) > 0:
							for j, lt in enumerate(gc_state_lts):
								if hg in key_hgs and lt in sample_lt_to_evalue[sample] and sample_lt_to_evalue[sample][lt] <= kq_evalue_threshold:
									features_key_hg = True
						if len(boundary_genes[sample].intersection(set(gc_state_lts).difference('background'))) > 0:
							boundary_lt_featured = True
						sample_gc_predictions.append([gc_state_lts, gc_state_hgs, len(gc_state_lts),
													   len(set(gc_state_hgs).difference("background")),
													   len(set(gc_state_hgs).difference("background").intersection(key_hgs)),
													   scaffold, boundary_lt_featured, features_key_hg, gcs_id])
						gcs_id += 1
					tmp = []
					hg_counter = 0
					dist_counter = 0
					active_search_mode = False
			if hg_counter >= 3:
				gc_state_lts = []
				gc_state_hgs = []
				begin_i = min(tmp)
				if last_hg_i == (len(hgs_ordered)-1):
					gc_state_lts = lts_ordered[begin_i:]
					gc_state_hgs = hgs_ordered[begin_i:]
				else:
					gc_state_lts = lts_ordered[min(tmp):last_hg_i + 1]
					gc_state_hgs = hgs_ordered[min(tmp):last_hg_i + 1]
				boundary_lt_featured = False
				features_key_hg = False
				if len(key_hgs.intersection(set(gc_state_hgs).difference('background'))) > 0:
					for j, lt in enumerate(gc_state_lts):
						if hg in key_hgs and lt in sample_lt_to_evalue[sample] and sample_lt_to_evalue[sample][
							lt] <= kq_evalue_threshold:
							features_key_hg = True
				if len(boundary_genes[sample].intersection(set(gc_state_lts).difference('background'))) > 0:
					boundary_lt_featured = True
				sample_gc_predictions.append([gc_state_lts, gc_state_hgs, len(gc_state_lts),
											  len(set(gc_state_hgs).difference("background")),
											  len(set(gc_state_hgs).difference("background").intersection(key_hgs)),
											  scaffold, boundary_lt_featured, features_key_hg, gcs_id])
				gcs_id += 1
	else:
		gcs_id = 1
		for scaffold in hgs_ordered_dict[sample]:
			hgs_ordered = hgs_ordered_dict[sample][scaffold]
			lts_ordered = lts_ordered_dict[sample][scaffold]
			hg_seq = numpy.array(list(hgs_ordered))
			hmm_predictions = model.predict(hg_seq)

			int_bg_coords = set([])
			tmp = []
			for i, hg_state in enumerate(hmm_predictions):
				if hg_state == 1:
					tmp.append(i)
				else:
					if len(tmp) <= max_int_genes_for_merge:
						for lt_index in tmp:
							int_bg_coords.add(lt_index)
					tmp = []

			gc_state_lts = []
			gc_state_hgs = []
			for i, hg_state in enumerate(hmm_predictions):
				lt = lts_ordered[i]
				hg = hgs_ordered[i]
				if hg_state == 0:
					gc_state_lts.append(lt)
					gc_state_hgs.append(hg)
				if hg_state == 1 or i == (len(hmm_predictions) - 1):
					if len(set(gc_state_hgs).difference('background')) >= 3:
						boundary_lt_featured = False
						features_key_hg = False
						if len(key_hgs.intersection(set(gc_state_hgs).difference('background'))) > 0:
							for j, lt in enumerate(gc_state_lts):
								if hg in key_hgs and lt in sample_lt_to_evalue[sample] and sample_lt_to_evalue[sample][lt] <= kq_evalue_threshold:
									features_key_hg = True
						if len(boundary_genes[sample].intersection(set(gc_state_lts).difference('background'))) > 0:
							boundary_lt_featured = True
						sample_gc_predictions.append([gc_state_lts, gc_state_hgs, len(gc_state_lts),
													   len(set(gc_state_hgs).difference("background")),
													   len(set(gc_state_hgs).difference("background").intersection(key_hgs)),
													   scaffold, boundary_lt_featured, features_key_hg, gcs_id])
						gcs_id += 1
					gc_state_lts = []
					gc_state_hgs = []

			gc_state_lts = []
			gc_state_hgs = []
			gc_states = []
			for i, hg_state in enumerate(hmm_predictions):
				lt = lts_ordered[i]
				hg = hgs_ordered[i]
				if hg_state == 0 or i in int_bg_coords:
					gc_state_lts.append(lt)
					gc_state_hgs.append(hg)
					gc_states.append(hg_state)
				elif hg_state == 1:
					if len(set(gc_state_hgs).difference('background')) >= 3:
						if '1' in set(gc_states):
							gc_state_lts = []
							gc_state_hgs = []
							gc_states = []
							continue
						boundary_lt_featured = False
						features_key_hg = False
						if len(key_hgs.intersection(set(gc_state_hgs).difference('background'))) > 0:
							for j, lt in enumerate(gc_state_lts):
								if hg in key_hgs and lt in sample_lt_to_evalue[sample] and sample_lt_to_evalue[sample][lt] <= kq_evalue_threshold:
									features_key_hg = True
						if len(boundary_genes[sample].intersection(set(gc_state_lts).difference('background'))) > 0:
							boundary_lt_featured = True
						sample_gc_predictions.append([gc_state_lts, gc_state_hgs, len(gc_state_lts),
													   len(set(gc_state_hgs).difference("background")),
													   len(set(gc_state_hgs).difference("background").intersection(key_hgs)),
													   scaffold, boundary_lt_featured, features_key_hg, gcs_id])

						gcs_id += 1
					gc_state_lts = []
					gc_state_hgs = []
					gc_states = []
				if i == (len(hmm_predictions) - 1):
					if len(set(gc_state_hgs).difference('background')) >= 3:
						if '1' in set(gc_states):
							gc_state_lts = []
							gc_state_hgs = []
							gc_states = []
							continue
						boundary_lt_featured = False
						features_key_hg = False
						if len(key_hgs.intersection(set(gc_state_hgs).difference('background'))) > 0:
							for j, lt in enumerate(gc_state_lts):
								if hg in key_hgs and lt in sample_lt_to_evalue[sample] and sample_lt_to_evalue[sample][lt] <= kq_evalue_threshold:
									features_key_hg = True
						if len(boundary_genes[sample].intersection(set(gc_state_lts).difference('background'))) > 0:
							boundary_lt_featured = True
						sample_gc_predictions.append([gc_state_lts, gc_state_hgs, len(gc_state_lts),
													   len(set(gc_state_hgs).difference("background")),
													   len(set(gc_state_hgs).difference("background").intersection(key_hgs)),
													   scaffold, boundary_lt_featured, features_key_hg, gcs_id])
						gcs_id += 1
					gc_state_lts = []
					gc_state_hgs = []
					gc_states = []

	if len(sample_gc_predictions) == 0: return

	sorted_sample_gc_predictions = [x for x in sorted(sample_gc_predictions, key=itemgetter(3), reverse=True)]

	cumulative_edge_hgs = set([])
	visited_scaffolds_with_edge_gc_segment = set([])
	sample_gc_predictions_filtered = []
	sample_edge_gc_predictions_filtered = []

	for gc_segment in sorted_sample_gc_predictions:
		if (gc_segment[3] >= min_hits and gc_segment[4] >= min_key_hits) or (gc_segment[7]) or (gc_segment[3] >= 3 and gc_segment[6] and not gc_segment[5] in visited_scaffolds_with_edge_gc_segment):
			# code to determine whether syntenically, the considered segment aligns with what is expected.
			# (skipped if input mode was 3)
			input_mode_3 = False
			if query_gene_info == None:
				input_mode_3 = True

			best_corr = 'irrelevant'
			if not input_mode_3 and syntenic_correlation_threshold > 0.0:
				best_corr = None
				segment_hg_order = []
				segment_hg_direction = []
				gc_hg_orders = defaultdict(list)
				gc_hg_directions = defaultdict(list)

				copy_count_of_hgs_in_segment = defaultdict(int)
				for hg in gc_segment[1]:
					copy_count_of_hgs_in_segment[hg] += 1

				for gi, g in enumerate(gc_segment[0]):
					hg = gc_segment[1][gi]
					if copy_count_of_hgs_in_segment[hg] != 1: continue
					gene_midpoint = (gene_locations[sample][g]['start'] + gene_locations[sample][g]['end']) / 2.0
					segment_hg_order.append(gene_midpoint)
					segment_hg_direction.append(gene_locations[sample][g]['direction'])

					for gc in query_gene_info:
						g_matching = []
						for g in query_gene_info[gc]:
							if g in lt_to_hg:
								hg_of_g = lt_to_hg[g]
								if hg_of_g == hg: g_matching.append(g)
						if len(g_matching) == 1:
							gc_gene_midpoint = (query_gene_info[gc][g_matching[0]]['start'] + query_gene_info[gc][g_matching[0]]['end']) / 2.0
							gc_hg_orders[gc].append(gc_gene_midpoint)
							gc_hg_directions[gc].append(query_gene_info[gc][g_matching[0]]['direction'])
						else:
							gc_hg_orders[gc].append(None)
							gc_hg_directions[gc].append(None)

				best_corr = None
				for gc in query_gene_info:
					try:
						assert (len(segment_hg_order) == len(gc_hg_orders[gc]))
						list1_same_dir = []
						list2_same_dir = []
						list1_comp_dir = []
						list2_comp_dir = []
						for iter, hgval1 in enumerate(segment_hg_order):
							hgdir1 = segment_hg_direction[iter]
							hgval2 = gc_hg_orders[gc][iter]
							hgdir2 = gc_hg_directions[gc][iter]
							if hgval1 == None or hgval2 == None: continue
							if hgdir1 == None or hgdir2 == None: continue
							if hgdir1 == hgdir2:
								list1_same_dir.append(hgval1)
								list2_same_dir.append(hgval2)
							else:
								list1_comp_dir.append(hgval1)
								list2_comp_dir.append(hgval2)

						if len(list1_same_dir) >= 3:
							corr, pval = pearsonr(list1_same_dir, list2_same_dir)
							corr = abs(corr)
							if (pval < 0.1) and ((best_corr and best_corr < corr) or (not best_corr)):
								best_corr = corr
						if len(list1_comp_dir) >= 3:
							corr, pval = pearsonr(list1_comp_dir, list2_comp_dir)
							corr = abs(corr)
							if (pval < 0.1) and ((best_corr and best_corr < corr) or (not best_corr)):
								best_corr = corr
					except:
						raise RuntimeError(traceback.format_exc())
						pass
			if best_corr != 'irrelevant' and (best_corr == None or best_corr < syntenic_correlation_threshold): continue
			if (gc_segment[3] >= min_hits and gc_segment[4] >= min_key_hits):
				sample_gc_predictions_filtered.append(gc_segment)
				if gc_segment[6]:
					cumulative_edge_hgs = cumulative_edge_hgs.union(set(gc_segment[1]))
					visited_scaffolds_with_edge_gc_segment.add(gc_segment[5])
			elif gc_segment[3] >= 3 and gc_segment[6] and not gc_segment[5] in visited_scaffolds_with_edge_gc_segment:
				sample_edge_gc_predictions_filtered.append(gc_segment)
				visited_scaffolds_with_edge_gc_segment.add(gc_segment[5])
				cumulative_edge_hgs = cumulative_edge_hgs.union(set(gc_segment[1]))

	if len(sample_edge_gc_predictions_filtered) >= 1 and draft_mode:
		if len(cumulative_edge_hgs) >= min_hits and len(cumulative_edge_hgs.intersection(key_hgs)) >= min_key_hits:
			sample_gc_predictions_filtered += sample_edge_gc_predictions_filtered

	# dereplicate nested segments!
	redundant_gcs = set([])
	for i, gcs1 in enumerate(sorted(sample_gc_predictions_filtered, key=itemgetter(2), reverse=True)):
		for j, gcs2 in enumerate(sorted(sample_gc_predictions_filtered, key=itemgetter(2), reverse=True)):
			if i < j and len(set(gcs1[0]).intersection(set(gcs2[0]))) > 0:
					redundant_gcs.add(gcs2[8])

	dereplicated_sample_gc_predictions_filtered = []
	for gcs in sample_gc_predictions_filtered:
		if gcs[8] in redundant_gcs: continue
		dereplicated_sample_gc_predictions_filtered.append(gcs)

	gc_sample_listing_handle = open(gc_info_dir + sample + '.bgcs.txt', 'w')
	gc_hg_evalue_handle = open(gc_info_dir + sample + '.hg_evalues.txt', 'w')

	sample_gc_id = 1
	for gc_segment in dereplicated_sample_gc_predictions_filtered:
		gc_genbank_file = gc_genbanks_dir + sample + '_fai-gene-cluster-' + str(sample_gc_id) + '.gbk'
		sample_gc_id += 1

		gc_segment_scaff = gc_segment[5]
		min_gc_pos = min([gene_locations[sample][g]['start'] for g in gc_segment[0]])-flanking_context
		max_gc_pos = max([gene_locations[sample][g]['end'] for g in gc_segment[0]])+flanking_context

		#print(gc_segment)
		#print([gene_locations[g]['start'] for g in gc_segment[0]])
		#print([gene_locations[g]['end'] for g in gc_segment[0]])
		util.createGenbank(target_annotation_info[sample]['genbank'], gc_genbank_file, gc_segment_scaff,
							  min_gc_pos, max_gc_pos)
		gc_sample_listing_handle.write('\t'.join([sample, gc_genbank_file]) + '\n')

		for i, lt in enumerate(gc_segment[0]):
			hg = gc_segment[1][i]
			identity = 0.0
			sqlratio = 0.0
			bitscore = 0.0
			if lt in sample_lt_to_identity[sample]: identity = sample_lt_to_identity[sample][lt]
			if lt in sample_lt_to_sqlratio[sample]: sqlratio = sample_lt_to_sqlratio[sample][lt]
			if lt in sample_lt_to_bitscore[sample]: bitscore = sample_lt_to_bitscore[sample][lt]
			gc_hg_evalue_handle.write('\t'.join([gc_genbank_file, sample, lt, hg, str(bitscore), str(identity), str(sqlratio), str(hg in key_hgs)]) + '\n')

	gc_hg_evalue_handle.close()
	gc_sample_listing_handle.close()

def filterParalogousSegmentsAndConcatenateIntoMultiRecordGenBanks(hmm_work_dir, homologous_gbk_dir, logObject):
	try:
		gbk_info_dir = hmm_work_dir + 'GeneCluster_Info/'
		gbk_filt_dir = hmm_work_dir + 'GeneCluster_Filtered_Segments/'
		util.setupReadyDirectory([gbk_filt_dir])

		sample_gc_segs = defaultdict(set)
		sample_gcs_segs_to_filter = defaultdict(set)
		sample_gcs_hgs = defaultdict(lambda: defaultdict(set))
		sample_gcs_hg_bs = defaultdict(lambda: defaultdict(float))
		for gcs_info in os.listdir(gbk_info_dir):
			if not gcs_info.endswith('.hg_evalues.txt'): continue
			sample = gcs_info.split('.hg_evalues.txt')[0]
			gcs_info_file = gbk_info_dir + gcs_info
			with open(gcs_info_file) as ogif:
				for line in ogif:
					gc_seg_gbk, _, lt, hg, bits, iden, sqlr, is_key_hg = line.split('\t')
					sample_gc_segs[sample].add(gc_seg_gbk)
					if hg != 'background':
						sample_gcs_hgs[sample][gc_seg_gbk].add(hg)
						sample_gcs_hg_bs[sample][gc_seg_gbk] += float(bits)

		# greedy removal of paralogs by prioritization of higher aggregate bitscore gene-cluster segments
		for sample in sample_gcs_hgs:
			gbk_filt_file = gbk_filt_dir + sample + '.txt'
			gbk_filt_handle = open(gbk_filt_file, 'w')
			for i, gcs1 in enumerate(sorted(sample_gcs_hg_bs[sample], key=itemgetter(1), reverse=True)):
				for j, gcs2 in enumerate(sorted(sample_gcs_hg_bs[sample], key=itemgetter(1), reverse=True)):
					if i >= j: continue
					gcs1_hg = sample_gcs_hgs[sample][gcs1]
					gcs2_hg = sample_gcs_hgs[sample][gcs2]
					# consider segments paralogous if more than 2 reference proteins/homolog groups are overlapping
					# suggesting paralogy beyond fragmentation that might have split a gene in two.
					intersection_hgs = gcs1_hg.intersection(gcs2_hg)
					if len(intersection_hgs) >= 2:
						sample_gcs_segs_to_filter[sample].add(gcs2)
						gbk_filt_handle.write(gcs2 + '\n')
			gbk_filt_handle.close()

		for sample in sample_gc_segs:
			retained_gcs = set([])
			for gcs in sample_gc_segs[sample]:
				if not gcs in sample_gcs_segs_to_filter[sample]:
					retained_gcs.add(gcs)
			if len(retained_gcs) == 1:
				shutil.move(list(retained_gcs)[0], homologous_gbk_dir)
			elif len(retained_gcs) > 1:
				sample_gbk_handle = open(homologous_gbk_dir + sample + '.gbk', 'a+')
				for gcs_gbk in retained_gcs:
					with open(gcs_gbk) as oggbk:
						for line in oggbk:
							sample_gbk_handle.write(line)
				sample_gbk_handle.close()
	except Exception as e:
		logObject.error('Issues resolving paralogous gene segments.')
		sys.stderr.write('Issues resolving paralogous gene segments.\n')
		sys.stderr.write(str(e) + '\n')
		raise RuntimeError(traceback.format_exc())
		sys.exit(1)

def plotOverviews(target_annotation_info, hmm_work_dir, protein_to_hg, plot_work_dir, plot_result_pdf, plot_naming_pdf, logObject, height=10, width=15):
	try:
		que_info_table_file = plot_work_dir + 'Name_Mapping.txt'
		que_info_table_handle = open(que_info_table_file, 'w')
		que_info_table_handle.write('NR_ID\tCDS_Locus_Tags\n')
		hg_prots = defaultdict(set)
		for p in protein_to_hg:
			hg_prots[protein_to_hg[p]].add(p)
		for hg in hg_prots:
			que_info_table_handle.write(hg + '\t' + ', '.join(sorted(list([x.split('|')[-1].strip() for x in hg_prots[hg]]))) + '\n')
		que_info_table_handle.close()

		gbk_info_dir = hmm_work_dir + 'GeneCluster_Info/'
		gbk_filt_dir = hmm_work_dir + 'GeneCluster_Filtered_Segments/'

		relevant_samples = set([])
		relevant_lts = set([])
		if os.path.isdir(gbk_info_dir):
			for gcs_info in os.listdir(gbk_info_dir):
				if not gcs_info.endswith('.hg_evalues.txt') or os.path.getsize(gbk_info_dir + gcs_info) == 0: continue
				sample = gcs_info.split('.hg_evalues.txt')[0]
				relevant_samples.add(sample)
				with open(gbk_info_dir + gcs_info) as ogif:
					for line in ogif:
						line = line.strip()
						ls = line.split('\t')
						relevant_lts.add(ls[2])
		lt_dists_to_scaff_starts = defaultdict(dict)
		lt_dists_to_scaff_ends = defaultdict(dict)
		lt_starts = defaultdict(dict)
		for sample in target_annotation_info:
			if not sample in relevant_samples: continue
			with open(target_annotation_info[sample]['genbank']) as osgbk:
				for rec in SeqIO.parse(osgbk, 'genbank'):
					scaff_len = len(str(rec.seq))
					scaff_is_relevant = False
					for feature in rec.features:
						if feature.type != 'CDS': continue
						lt = feature.qualifiers.get('locus_tag')[0]
						if lt in relevant_lts:
							scaff_is_relevant = True
							break
					if not scaff_is_relevant: continue
					for feature in rec.features:
						lt = feature.qualifiers.get('locus_tag')[0]
						all_coords = []
						if not 'join' in str(feature.location):
							start = min([int(x.strip('>').strip('<')) for x in
										 str(feature.location)[1:].split(']')[0].split(':')]) + 1
							end = max([int(x.strip('>').strip('<')) for x in
									   str(feature.location)[1:].split(']')[0].split(':')])
							direction = str(feature.location).split('(')[1].split(')')[0]
							all_coords.append([start, end, direction])
						else:
							for exon_coord in str(feature.location)[5:-1].split(', '):
								start = min([int(x.strip('>').strip('<')) for x in
											 exon_coord[1:].split(']')[0].split(':')]) + 1
								end = max([int(x.strip('>').strip('<')) for x in
										   exon_coord[1:].split(']')[0].split(':')])
								direction = exon_coord.split('(')[1].split(')')[0]
								all_coords.append([start, end, direction])
						start = 1e16
						end = -1
						for sc, ec, dc in sorted(all_coords, key=itemgetter(0), reverse=False):
							if sc < start:
								start = sc
							if ec > end:
								end = ec
						lt_dists_to_scaff_starts[sample][lt] = start
						lt_dists_to_scaff_ends[sample][lt] = scaff_len - end
						lt_starts[sample][lt] = start

		sample_filt_gcs = defaultdict(set)
		if os.path.isdir(gbk_filt_dir):
			for gfi in os.listdir(gbk_filt_dir):
				sample = gfi.split('.txt')[0]
				with open(gbk_filt_dir + gfi) as ogfi:
					for gcs in ogfi:
						sample_filt_gcs[sample].add(gcs.strip())

		plot_input_file = plot_work_dir + 'Bar_Plot_Input.txt'
		plot_input_handle = open(plot_input_file, 'w')
		plot_input_handle.write('\t'.join(['sample', 'segment_title', 'gene', 'key', 'gene_order', 'identity', 'sql_ratio']) + '\n')
		for gcs_info in os.listdir(gbk_info_dir):
			if not gcs_info.endswith('.hg_evalues.txt'): continue
			sample = gcs_info.split('.hg_evalues.txt')[0]
			gcs_info_file = gbk_info_dir + gcs_info

			dists_to_scaff_start = defaultdict(list)
			dists_to_scaff_end = defaultdict(list)
			filt_status = 'Retained'
			with open(gcs_info_file) as ogif:
				for line in ogif:
					gc_seg_gbk, _, lt, hg, bits, iden, sqlr, is_key_hg = line.strip().split('\t')
					dists_to_scaff_start[gc_seg_gbk].append(lt_dists_to_scaff_starts[sample][lt])
					dists_to_scaff_end[gc_seg_gbk].append(lt_dists_to_scaff_ends[sample][lt])

			seg_to_name = {}
			for gc_seg_gbk in dists_to_scaff_start:
				gcs_id = gc_seg_gbk.split('/')[-1].split('fai-gene-cluster-')[1].split('.gbk')[0]
				filt_status = 'Retained'
				if gc_seg_gbk in sample_filt_gcs[sample]:
					filt_status = 'Filtered'
				gcs_name = sample + ' Segment_' + gcs_id + ' (' + filt_status + '); ' + str(min(dists_to_scaff_start[gc_seg_gbk])) + ' bp to Scaffold Start; ' + str(min(dists_to_scaff_end[gc_seg_gbk])) + ' bp to Scaffold End'
				seg_to_name[gc_seg_gbk] = gcs_name

			hg_iters = defaultdict(int)
			with open(gcs_info_file) as ogif:
				for line in ogif:
					gc_seg_gbk, _, lt, hg, bits, iden, sqlr, is_key_hg = line.strip().split('\t')
					if hg == 'background':
						hg = 'nov_aux'
					hg_iters[hg] += 1
					hg_id = hg + '|' + str(hg_iters[hg])
					plot_input_handle.write('\t'.join([sample, seg_to_name[gc_seg_gbk], hg_id, is_key_hg, str(lt_starts[sample][lt]), str(iden), str(sqlr)]) + '\n')
		plot_input_handle.close()

		plot_cmd = ['Rscript', plot_prog, que_info_table_file, plot_input_file, plot_result_pdf, plot_naming_pdf, str(height), str(width)]
		try:
			subprocess.call(' '.join(plot_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(plot_result_pdf))
			logObject.info('Successfully ran: %s' % ' '.join(plot_cmd))
		except Exception as e:
			logObject.error('Had an issue running R based plotting: %s' % ' '.join(plot_cmd))
			sys.stderr.write('Had an issue running R based plotting: %s\n' % ' '.join(plot_cmd))
			logObject.error(e)
			sys.exit(1)

	except Exception as e:
		logObject.error('Issues with plotting overviews of homologous gene-cluster segments identified.')
		sys.stderr.write('Issues with plotting overviews of homologous gene-cluster segments identified.\n')
		sys.stderr.write(str(e) + '\n')
		raise RuntimeError(traceback.format_exc())
		sys.exit(1)