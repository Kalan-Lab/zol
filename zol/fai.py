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
from scipy.stats import pearsonr

def subsetGenBankForQueryLocus(full_gw_genbank, locus_genbank, locus_proteins, reference_contig, reference_start, reference_end, logObject):
	try:
		util.createBGCGenbank(full_gw_genbank, locus_genbank, reference_contig, reference_start, reference_end)

		locus_proteins_handle = open(locus_proteins, 'w')
		with open(locus_genbank) as olg:
			for rec in SeqIO.parse(olg, 'fasta'):
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

		# maybe make more robust later - assumption that representative is first in clusters file
		cluster_rep = None
		with open(cdhit_cluster_file) as occf:
			for line in occf:
				line = line.strip()
				if line.startswith('>'): continue
				ls = line.split()
				lt = ls[2][1:-3]
				if line.endswith(' *'):
					cluster_rep = lt
				protein_to_hg[lt] = cluster_rep
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
					proteins, nucleotides = util.parseGenbankForCDSProteinsAndDNA(gbk, logObject)
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
		prot_algn_dir = proc_dir + 'HG_Protein_Alignments/'
		phmm_dir = proc_dir + 'HG_Profile_HMMs/'
		cons_dir = proc_dir + 'HG_Consensus_Sequences/'
		consensus_prot_seqs_faa = outdir + 'HG_Consensus_Seqs.faa'
		util.setupReadyDirectory([proc_dir, prot_algn_dir, phmm_dir, cons_dir, hg_prot_dir])
		zol.partitionSequencesByHomologGroups(ortho_matrix_file, prot_dir, nucl_dir, hg_prot_dir, logObject)
		zol.createProteinAlignments(hg_prot_dir, prot_algn_dir, logObject, use_super5=use_super5, cpus=cpus)
		zol.createProfileHMMsAndConsensusSeqs(prot_algn_dir, phmm_dir, cons_dir, logObject, cpus=1)
		consensus_prot_seqs_handle = open(consensus_prot_seqs_faa, 'w')
		for f in os.listdir(cons_dir):
			with open(cons_dir + f) as ocf:
				for rec in SeqIO.parse(ocf, 'fasta'):
					consensus_prot_seqs_handle.write('>' + f.split('.cons.faa')[0] + '\n' + str(rec.seq) + '\n')
		consensus_prot_seqs_handle.close()
		return([ortho_matrix_file, consensus_prot_seqs_faa])

def processGenomeWideGenbanks(target_annot_listing_file, logObject, cpus=1):
	gene_locations = {}
	scaffold_genes = {}
	boundary_genes = {}
	gene_id_to_order = {}
	gene_order_to_id = {}
	try:
		with multiprocessing.Manager() as manager:
			sample_gbk_info = manager.dict()
			genbanks = []
			for sample in target_annot_listing_file:
				sample_genbank = target_annot_listing_file[sample]['genbank']
				genbanks.append([sample, sample_genbank, sample_gbk_info])

			with manager.Pool(cpus) as pool:
				pool.map(util.parseGenbankAndFindBoundaryGenes, genbanks)

			for sample in sample_gbk_info:
				gene_location, scaff_genes, bound_genes, gito, goti = sample_gbk_info[sample]
				gene_locations[sample] = gene_location
				scaffold_genes[sample] = scaff_genes
				boundary_genes[sample] = bound_genes
				gene_id_to_order[sample] = gito
				gene_order_to_id[sample] = goti
	except Exception as e:
		logObject.error('Issues with parsing CDS location information from genes of target genomes.')
		sys.stderr.write('Issues with parsing CDS location information from genes of target genomes.\n')
		sys.exit(1)
		sys.stderr.write(str(e) + '\n')

	target_genome_gene_info = {'gene_locations': gene_locations, 'scaffold_genes': scaffold_genes,
							   'boundary_genes': boundary_genes, 'gene_id_to_order': gene_id_to_order,
							   'gene_order_to_id': gene_order_to_id}
	return(target_genome_gene_info)

def runDiamondBlastp(target_annot_listing_file, query_fasta, work_dir, logObject, evalue_cutoff=1e-10, cpus=1):
	"""
	Function to run Diamond blastp and process results
	"""
	diamond_results = []
	try:
		search_res_dir = work_dir + 'Alignment_Results/'
		util.setupReadyDirectory([search_res_dir])

		alignment_cmds = []
		for sample in target_annot_listing_file:
			sample_proteome = target_annot_listing_file[sample]['predicted_proteome']
			sample_proteome_db = sample_proteome.split('.faa')[0] + '.dmnd'
			result_file = search_res_dir + sample + '.txt'
			diamond_cmd = ['diamond', 'makedb', '--in', sample_proteome, '-d', sample_proteome_db, ';', 'diamond',
						   'blastp', '--threads', '1', '--very-sensitive', '--query', query_fasta, '--db',
						   sample_proteome_db, '--outfmt', '6', '--out', result_file, '--evalue', evalue_cutoff,
						   logObject]
			alignment_cmds.append(diamond_cmd)

		p = multiprocessing.Pool(cpus)
		p.map(util.multiProcess, alignment_cmds)
		p.close()

		for sample in target_annot_listing_file:
			result_file = search_res_dir + sample + '.txt'
			assert (os.path.isfile(result_file))

			best_hit_per_lt = defaultdict(lambda: [set([]), 100000.0])
			with open(result_file) as orf:
				for line in orf:
					line = line.strip()
					ls = line.split()
					hg = ls[0]
					lt = ls[1]
					eval = decimal.Decimal(ls[10])
					if eval < best_hit_per_lt[lt][1]:
						best_hit_per_lt[lt][1] = eval
						best_hit_per_lt[lt][0] = set([hg])
					elif eval == best_hit_per_lt[lt][1]:
						best_hit_per_lt[lt][0].add(hg)

			for lt in best_hit_per_lt:
				for hg in best_hit_per_lt[lt][0]:
					diamond_results[lt].append([hg, best_hit_per_lt[lt][1], sample])
	except Exception as e:
		logObject.error('Issues with running DIAMOND blastp or processing of results.')
		sys.stderr.write('Issues with running DIAMOND blastp or processing of results.\n')
		sys.exit(1)
		sys.stderr.write(str(e) + '\n')

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
		diamond_cmd = ['diamond', 'makedb', '--in', query_dmnd_db, '-d', query_dmnd_db, ';', 'diamond',
					   'blastp', '--threads', '1', '--very-sensitive', '--query', key_protein_queries_fasta, '--db',
					   query_dmnd_db, '--outfmt', '6', '--out', align_result_file, '--evalue', '1e-10']

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
					sys.stderr.write(str(2) + '\n')
					sys.exit(1)

	except Exception as e:
		logObject.error('Issues with mapping key proteins to general proteins or consensus sequence of homolog groups.')
		sys.stderr.write('Issues with mapping key proteins to general proteins or consensus sequence of homolog groups.\n')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)
	return key_hgs

def identifyGCFInstances(query_information, target_information, diamond_results, work_dir, logObject, min_hits=5, min_key_hits=3,
						 gc_to_gc_transition_prob=0.9, bg_to_bg_transition_prob=0.9, gc_emission_prob_with_hit=0.95,
						 gc_emission_prob_without_hit=0.2, syntenic_correlation_threshold=0.8, kq_evalue_threshold=1e-20,
						 cpus=1, block_size=3000):
	"""
	Function to search for instances of Gene Cluster in samples using HMM based approach based on homolog groups as
	characters, "part of Gene Cluster" and "not part of Gene Cluster" as states - all trained on initial BGCs
	constituting GCF as identified by lsaBGC-Cluster.py. This function utilizes the convenient Python library
	Pomegranate.
	:param query_information: Dictionary with 3 keys: protein_to_hg (dict), query_fasta (file path), comp_gene_info (dict)
	:param target_information: Dictionary with 2 keys: target_annotation_information (dict), target_genome_gene_info (dict)
	:param diamond_results: Dictionary of DIAMOND alignment results meeting e-value threshold
	"""

	try:
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
		bg_state = State(bg_distribution, name='Background')

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
			for lt in diamond_results:
				for i, hits in enumerate(sorted(diamond_results[lt], key=itemgetter(1))):
					if i == 0 and hits[2] in block_samp_set:
						sample_lt_to_hg[hits[2]][lt] = hits[0]
						sample_hgs[hits[2]].add(hits[0])
						sample_lt_to_evalue[hits[2]][lt] = decimal.Decimal(hits[1])

			identify_gcf_segments_input = []
			for sample in sample_hgs:
				if len(sample_hgs[sample]) < 3: continue
				hgs_ordered_dict = {}
				lts_ordered_dict = {}
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
					hgs_ordered_dict[scaffold] = hgs_ordered
					lts_ordered_dict[scaffold] = lts_ordered

				identify_gcf_segments_input.append([gc_info_dir, gc_genbanks_dir, sample, target_annotation_info[sample],
													sample_lt_to_evalue[sample], model, lts_ordered_dict, hgs_ordered_dict,
													query_gene_info, dict(gene_locations[sample]),
													dict(gene_id_to_order[sample]), dict(gene_order_to_id[sample]),
													boundary_genes[sample], lt_to_hg, min_hits, min_key_hits, key_hgs,
													kq_evalue_threshold, syntenic_correlation_threshold])
			with multiprocessing.Manager() as manager:
				with manager.Pool(cpus) as pool:
					pool.map(identify_gcf_instances, identify_gcf_segments_input)

		os.system('find %s -type f -name "*.bgcs.txt" -exec cat {} + >> %s' % (gc_info_dir, gc_list_file))

		os.system('find %s -type f -name "*.hg_evalues.txt" -exec cat {} + >> %s' % (gc_info_dir, gc_hmm_evalues_file))

	except Exception as e:
		logObject.error('Issues with managing running of HMM to find gene-cluster homolog segments.')
		sys.stderr.write('Issues with managing running of HMM to find gene-cluster homolog segments.\n')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def identify_gcf_instances(input_args):
	gc_info_dir, gc_genbanks_dir, sample, target_annotation_info, sample_lt_to_evalue, model, lts_ordered_dict, hgs_ordered_dict, query_gene_info, gene_locations, gene_id_to_order, gene_order_to_id, boundary_genes, lt_to_hg, min_hits, min_key_hits, key_hgs, kq_evalue_threshold, syntenic_correlation_threshold = input_args

	sample_gc_predictions = []
	for scaffold in hgs_ordered_dict:
		hgs_ordered = hgs_ordered_dict[scaffold]
		lts_ordered = lts_ordered_dict[scaffold]
		hg_seq = numpy.array(list(hgs_ordered))
		hmm_predictions = model.predict(hg_seq)

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
					if len(key_hgs.intersection(set(gcf_state_hgs).difference('background'))) > 0: features_key_hg = True
					if len(boundary_genes.intersection(set(gcf_state_lts).difference('background'))) > 0: boundary_lt_featured = True
					sample_gc_predictions.append([gcf_state_lts, gcf_state_hgs, len(gcf_state_lts),
												   len(set(gcf_state_hgs).difference("background")),
												   len(set(gcf_state_hgs).difference("background").intersection(key_hgs)),
												   scaffold, boundary_lt_featured, features_key_hg])
				gcf_state_lts = []
				gcf_state_hgs = []

	if len(sample_gc_predictions) == 0: return

	sorted_sample_gc_predictions = [x for x in sorted(sample_gc_predictions, key=itemgetter(3), reverse=True)]

	cumulative_edge_hgs = set([])
	visited_scaffolds_with_edge_gcf_segment = set([])
	sample_gcf_predictions_filtered = []
	sample_edge_gcf_predictions_filtered = []

	for gc_segment in sorted_sample_gc_predictions:
		if (gc_segment[3] >= min_hits and gc_segment[4] >= min_key_hits) or (gc_segment[-1]) or (gc_segment[3] >= 3 and gc_segment[-2] and not gc_segment[5] in visited_scaffolds_with_edge_gcf_segment):
			# code to determine whether syntenically, the considered segment aligns with what is expected.
			# (skipped if input mode was 3)

			input_mode_3 = False
			if query_gene_info == None:
				input_mode_3 = True

			best_corr = 'irrelevant'
			if not input_mode_3:
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
					gene_midpoint = (gene_locations[g]['start'] + gene_locations[g]['end']) / 2.0
					segment_hg_order.append(gene_midpoint)
					segment_hg_direction.append(gene_locations[g]['direction'])

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
						assert (len(segment_hg_order) == len(query_gene_info[gc]))
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
						pass

			if best_corr == None or best_corr < syntenic_correlation_threshold: continue

			if (gc_segment[3] >= min_hits and gc_segment[4] >= min_key_hits) or (gc_segment[-1]):
				sample_gcf_predictions_filtered.append(gc_segment)
				if gc_segment[-2]:
					cumulative_edge_hgs = cumulative_edge_hgs.union(set(gc_segment[1]))
					visited_scaffolds_with_edge_gcf_segment.add(gc_segment[5])
			elif gc_segment[3] >= 3 and gc_segment[-2] and not gc_segment[5] in visited_scaffolds_with_edge_gcf_segment:
				sample_edge_gcf_predictions_filtered.append(gc_segment)
				visited_scaffolds_with_edge_gcf_segment.add(gc_segment[5])
				cumulative_edge_hgs = cumulative_edge_hgs.union(set(gc_segment[1]))

	if len(sample_edge_gcf_predictions_filtered) >= 1:
		if len(cumulative_edge_hgs) >= min_hits and len(cumulative_edge_hgs.intersection(key_hgs)) >= min_key_hits:
			sample_gcf_predictions_filtered += sample_edge_gcf_predictions_filtered

	gc_sample_listing_handle = open(gc_info_dir + sample + '.bgcs.txt', 'w')
	gc_hg_evalue_handle = open(gc_info_dir + sample + '.hg_evalues.txt', 'w')

	sample_gc_id = 1
	for gc_segment in enumerate(sample_gcf_predictions_filtered):
		gc_genbank_file = gc_genbanks_dir + sample + '_fai-gene-cluster_-' + str(sample_gc_id) + '.gbk'
		sample_gc_id += 1

		gc_segment_scaff = gc_segment[5]
		min_gc_pos = min([gene_locations[g]['start'] for g in gc_segment[0]])
		max_gc_pos = max([gene_locations[g]['end'] for g in gc_segment[0]])

		util.createBGCGenbank(target_annotation_info['genbank'], gc_genbank_file, gc_segment_scaff,
							  min_gc_pos, max_gc_pos)
		gc_sample_listing_handle.write('\t'.join([sample, gc_genbank_file]) + '\n')

		for i, lt in enumerate(gc_segment[0]):
			hg = gc_segment[1][i]
			evalue = decimal.Decimal(100000.0)
			if lt in sample_lt_to_evalue: evalue = sample_lt_to_evalue[lt]
			gc_hg_evalue_handle.write('\t'.join([gc_genbank_file, sample, lt, hg, str(evalue), str(hg in key_hgs)]) + '\n')

	gc_hg_evalue_handle.close()
	gc_sample_listing_handle.close()

def concatenateBGCsIntoMultiRecordGenBanks(gc_segments_dir, fin_outdir, logObject):
	try:
		for gcs in os.listdir(gc_segments_dir):
			sample = gcs.split('_fai-gene-cluster')[0]
			sample_gbk_handle = open(fin_outdir + sample + '.gbk', 'a+')
			with open(gc_segments_dir + gcs) as oggbk:
				for line in oggbk:
					sample_gbk_handle.write(line)
			sample_gbk_handle.close()
	except Exception as e:
		logObject.error('Issues concatenating multiple gene-cluster segments per genome/sample into a multi-record GenBank.')
		sys.stderr.write('Issues concatenating multiple gene-cluster segments per genome/sample into a multi-record GenBank.\n')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)
