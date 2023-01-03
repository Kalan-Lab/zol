import copy
import os
import sys
import logging
import traceback
import statistics
import random
import subprocess
import pysam
import gzip
import multiprocessing
from scipy.stats import f_oneway, fisher_exact, pearsonr, median_abs_deviation, entropy
from ete3 import Tree
import numpy as np
from operator import itemgetter
import itertools
from collections import defaultdict
from zol.classes.Pan import Pan
from zol import util
from pandas import DataFrame
from pomegranate import *
import math
import warnings
import decimal
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from zol.classes.BGC import BGC

warnings.filterwarnings("ignore")
# updated 07/22/2022 to have key term be transpos instead of transp - because of transporter now appearing in definitions
# due to switch from Prokka annotation to custom KO/PGAP annotations
mges = set(['transpos', 'integrase'])
purine_alleles = set(['A', 'G'])

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-3])
RSCRIPT_FOR_COLORBREW = lsaBGC_main_directory + '/zol/Rscripts/brewColors.R'
RSCRIPT_FOR_BGSEE = lsaBGC_main_directory + '/zol/Rscripts/clusterHeatmap.R'
RSCRIPT_FOR_CLUSTER_ASSESSMENT_PLOTTING = lsaBGC_main_directory + '/zol/Rscripts/generatePopGenePlots.R'
RSCRIPT_FOR_TAJIMA = lsaBGC_main_directory + '/zol/Rscripts/calculateTajimasD.R'
RSCRIPT_FOR_GENERATE = lsaBGC_main_directory + '/zol/Rscripts/GeneRatePhylogeny.R'
RSCRIPT_FOR_PCA = lsaBGC_main_directory + '/zol/Rscripts/ClusterVisualOfSamples.R'
SEED = 1234

class GCF(Pan):
	def __init__(self, bgc_genbanks_listing, gcf_id='GCF_X', logObject=None, lineage_name='Unnamed lineage'):
		super().__init__(bgc_genbanks_listing, lineage_name=lineage_name, logObject=logObject)
		self.gcf_id = gcf_id

		#######
		## Variables not set during initialization
		#######

		# General variables
		self.hg_to_color = None
		self.hg_order_scores = defaultdict(lambda: ['NA', 'NA'])
		self.specific_core_homologs =set([])
		self.scc_homologs = set([])
		self.core_homologs = set([])
		self.protocluster_core_homologs = set([])

		# Sequence and alignment directories
		self.nucl_seq_dir = None
		self.prot_seq_dir = None
		self.prot_alg_dir = None
		self.codo_alg_dir = None
		self.nucl_filt_seq_dir = None
		self.prot_filt_seq_dir = None
		self.prot_filt_alg_dir = None
		self.codo_filt_alg_dir = None

		# Concatenated HMMER3 HMM profiles database of homolog groups in GCF
		self.concatenated_profile_HMM = None

		# Dictionary of individual genes to haplotype/allelic representative gene
		self.instance_to_haplotype = {}

		# Set of samples with sequencing reads to avoid for reporting alleles and novel SNVs,
		# these samples do not exhibit enough support for harboring a full BGC for the GC
		self.avoid_samples = set([])

	def identifyKeyHomologGroups(self, all_primary=False):
		try:
			initial_samples_with_at_least_one_gcf_hg = set([])
			for hg in self.hg_genes:
				for gene in self.hg_genes[hg]:
					gene_info = self.comp_gene_info[gene]
					if not gene_info['is_expansion_bgc'] or all_primary:
						bgc_id = gene_info['bgc_name']
						sample_id = self.bgc_sample[bgc_id]
						initial_samples_with_at_least_one_gcf_hg.add(sample_id)

			for hg in self.hg_genes:
				sample_counts = defaultdict(int)
				sample_with_hg_as_protocluster_core = 0
				for gene in self.hg_genes[hg]:
					gene_info = self.comp_gene_info[gene]
					if not gene_info['is_expansion_bgc'] or all_primary:
						bgc_id = gene_info['bgc_name']
						sample_id = self.bgc_sample[bgc_id]
						sample_counts[sample_id] += 1
						if gene_info['core_overlap']:
							sample_with_hg_as_protocluster_core += 1


				samples_with_single_copy = set([s[0] for s in sample_counts.items() if s[1] == 1])
				samples_with_any_copy = set([s[0] for s in sample_counts.items() if s[1] > 0])

				# check that hg is single-copy-core or just core
				if len(samples_with_single_copy.symmetric_difference(initial_samples_with_at_least_one_gcf_hg)) == 0:
					self.scc_homologs.add(hg)
				if len(samples_with_any_copy.symmetric_difference(initial_samples_with_at_least_one_gcf_hg)) == 0:
					self.core_homologs.add(hg)

				if len(samples_with_any_copy) > 0 and float(sample_with_hg_as_protocluster_core)/len(samples_with_any_copy) >= 0.5:
					self.protocluster_core_homologs.add(hg)

			self.logObject.info("Conserved Core Set of Homolog Groups:\t" + '; '.join(self.core_homologs))
			self.logObject.info("Conserved Single-Copy Set of Homolog Groups:\t" + '; '.join(self.scc_homologs))
			self.logObject.info("Proto-Core / Rule Based Homolog Groups:\t" + '; '.join(self.protocluster_core_homologs))
		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues with identifying key homolog groups.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def modifyPhylogenyForSamplesWithMultipleBGCs(self, input_phylogeny, result_phylogeny, prune_set=None):
		"""
		Function which takes in an input phylogeny and produces a replicate resulting phylogeny with samples/leafs which
		have multiple BGC instances for a GCF expanded.

		:param input_phylogeny: input newick phylogeny file
		:result result_phylogeny: resulting newick phylogeny file
		"""
		try:
			number_of_added_leaves = 0
			t = Tree(input_phylogeny)
			if prune_set != None:
				t.prune(prune_set)

			for node in t.traverse('postorder'):
				if node.name in self.sample_bgcs and len(self.sample_bgcs[node.name]) > 1:
					og_node_name = node.name
					node.name = node.name + '_INNERNODE'
					for bgc_id in self.sample_bgcs[og_node_name]:
						# if bgc_id == node.name: continue
						node.add_child(name=bgc_id)
						child_node = t.search_nodes(name=bgc_id)[0]
						child_node.dist = 0
						if bgc_id != og_node_name: number_of_added_leaves += 1

			t.write(format=0, outfile=result_phylogeny)
			if self.logObject:
				self.logObject.info(
					"New phylogeny with an additional %d leafs to reflect samples with multiple BGCs can be found at: %s." % (
						number_of_added_leaves, result_phylogeny))
		except Exception as e:
			if self.logObject:
				self.logObject.error(
					"Had difficulties properly editing phylogeny to duplicate leafs for samples with multiple BGCs for GCF.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def assignColorsToHGs(self, gene_to_hg, bgc_genes, outdir):
		"""
		Simple function to associate each homolog group with a color for consistent coloring.

		:param gene_to_hg: gene to HG relationship.
		:param bgc_genes:  set of genes per HG.
		:return: dictionary mapping each HG to a hex color value.
		"""

		hg_bgc_counts = defaultdict(int)
		for b in bgc_genes:
			for g in bgc_genes[b]:
				if g in gene_to_hg:
					hg_bgc_counts[gene_to_hg[g]] += 1

		hgs = set([])
		for c in hg_bgc_counts:
			if hg_bgc_counts[c] > 1:
				hgs.add(c)

		len_hgs = len(hgs)
		color_listing_file = outdir + 'colors_for_hgs.txt'

		rscript_brew_color = ["Rscript", RSCRIPT_FOR_COLORBREW, str(len_hgs), color_listing_file]
		if self.logObject:
			self.logObject.info('Running R-based plotting with the following command: %s' % ' '.join(rscript_brew_color))
		try:
			subprocess.call(' '.join(rscript_brew_color), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert(os.path.isfile(color_listing_file) and os.path.getsize(color_listing_file) > 0)
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(rscript_brew_color))

		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(rscript_brew_color))
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

		# read in list of colors
		colors = []
		with open(color_listing_file) as ocf:
			colors = [x.strip() for x in ocf.readlines()]
		random.Random(SEED).shuffle(colors)

		hg_to_color = {}
		for i, c in enumerate(set(hgs)):
			hg_to_color[c] = colors[i]
		self.hg_to_color = hg_to_color

	def createItolBGCSeeTrack(self, result_track_file):
		"""
		Function to create a track file for visualizing BGC gene architecture across a phylogeny in the interactive tree
		of life (iTol)

		:param result_track_file: The path to the resulting iTol track file for BGC gene visualization.
		"""
		try:
			track_handle = open(result_track_file, 'w')

			if self.logObject:
				self.logObject.info("Writing iTol track file to: %s" % result_track_file)
				self.logObject.info("Track will have label: %s" % self.gcf_id)

			# write header for iTol track file
			track_handle.write('DATASET_DOMAINS\n')
			track_handle.write('SEPARATOR TAB\n')
			track_handle.write('DATASET_LABEL\t%s\n' % self.gcf_id)
			track_handle.write('COLOR\t#000000\n')
			track_handle.write('BORDER_WIDTH\t1\n')
			track_handle.write('BORDER_COLOR\t#000000\n')
			track_handle.write('SHOW_DOMAIN_LABELS\t0\n')
			track_handle.write('DATA\n')

			# write the rest of the iTol track file for illustrating genes across BGC instances
			ref_hg_directions = {}
			bgc_gene_counts = defaultdict(int)
			for bgc in self.bgc_genes:
				bgc_gene_counts[bgc] = len(self.bgc_genes[bgc])

			for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
				bgc = item[0]
				curr_bgc_genes = self.bgc_genes[bgc]
				last_gene_end = max([self.comp_gene_info[lt]['end'] for lt in curr_bgc_genes])
				printlist = [bgc, str(last_gene_end)]
				hg_directions = {}
				hg_lengths = defaultdict(list)
				for lt in curr_bgc_genes:
					ginfo = self.comp_gene_info[lt]
					hg = 'singleton'
					if lt in self.gene_to_hg:
						hg = self.gene_to_hg[lt]
					shape = 'None'
					if ginfo['direction'] == '+':
						shape = 'TR'
					elif ginfo['direction'] == '-':
						shape = 'TL'
					gstart = ginfo['start']
					gend = ginfo['end']
					hg_color = "#dbdbdb"
					if hg in self.hg_to_color:
						hg_color = self.hg_to_color[hg]
					gene_string = '|'.join([str(x) for x in [shape, gstart, gend, hg_color, hg]])
					printlist.append(gene_string)
					if hg != 'singleton':
						hg_directions[hg] = ginfo['direction']
						hg_lengths[hg].append(gend - gstart)
				if i == 0:
					ref_hg_directions = hg_directions
					track_handle.write('\t'.join(printlist) + '\n')
				else:
					flip_support = 0
					keep_support = 0
					for c in ref_hg_directions:
						if not c in hg_directions: continue
						hg_weight = statistics.mean(hg_lengths[c])
						if hg_directions[c] == ref_hg_directions[c]:
							keep_support += hg_weight
						else:
							flip_support += hg_weight

					# flip the genbank visual if necessary, first BGC processed is used as reference guide
					if flip_support > keep_support:
						flip_printlist = printlist[:2]
						for gene_string in printlist[2:]:
							gene_info = gene_string.split('|')
							new_shape = None
							if gene_info[0] == 'TR':
								new_shape = 'TL'
							elif gene_info[0] == 'TL':
								new_shape = 'TR'
							new_gstart = int(last_gene_end) - int(gene_info[2])
							new_gend = int(last_gene_end) - int(gene_info[1])
							new_gene_info = '|'.join([new_shape, str(new_gstart), str(new_gend)] + gene_info[-2:])
							flip_printlist.append(new_gene_info)
						track_handle.write('\t'.join(flip_printlist) + '\n')
					else:
						track_handle.write('\t'.join(printlist) + '\n')
			track_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error("Had difficulties creating iTol track for visualization of BGC gene architecture.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def visualizeGCFViaR(self, gggenes_track_file, heatmap_track_file, phylogeny_file, result_pdf_file):
		"""
		Function to create tracks for visualization of gene architecture of BGCs belonging to GCF and run Rscript clusterHeatmap.R
		to produce automatic PDFs of plots. In addition, clusterHeatmap.R also produces a heatmap to more easily identify homolog
		groups which are conserved across isolates found to feature GCF.

		:param gggenes_track_file: Path to file with gggenes track information (will be created/written to by function, if it doesn't exist!)
		:param heatmap_track_file: Path to file for heatmap visual component (will be created/written to by function, if it doesn't exist!)
		:param phylogeny_file: Phylogeny to use for visualization.
		:param result_pdf_file: Path to PDF file where plots from clusterHeatmap.R will be written to.
		"""
		try:
			if os.path.isfile(gggenes_track_file) or os.path.isfile(heatmap_track_file):
				os.system('rm -f %s %s' % (gggenes_track_file, heatmap_track_file))
			gggenes_track_handle = open(gggenes_track_file, 'w')
			heatmap_track_handle = open(heatmap_track_file, 'w')
			if self.logObject:
				self.logObject.info("Writing gggenes input file to: %s" % gggenes_track_file)
				self.logObject.info("Writing heatmap input file to: %s" % heatmap_track_file)
			# write header for track files
			gggenes_track_handle.write('label\tgene\tstart\tend\tforward\tog\tog_color\n')
			heatmap_track_handle.write('label\tog\tog_presence\tog_count\n')

			ref_hg_directions = {}

			bgc_gene_counts = defaultdict(int)
			for bgc in self.bgc_genes:
				bgc_gene_counts[bgc] = len(self.bgc_genes[bgc])

			tree_obj = Tree(phylogeny_file)
			bgc_weights = defaultdict(int)
			all_bgcs_in_tree = set([])
			for leaf in tree_obj:
				all_bgcs_in_tree.add(str(leaf).strip('\n').lstrip('-'))
				bgc_weights[str(leaf).strip('\n').lstrip('-')] += 1

			bgc_hg_presence = defaultdict(lambda: defaultdict(lambda: 'Absent'))
			hg_counts = defaultdict(int)
			for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
				bgc = item[0]
				if not bgc in all_bgcs_in_tree: continue
				curr_bgc_genes = self.bgc_genes[bgc]
				last_gene_end = max([self.comp_gene_info[lt]['end'] for lt in curr_bgc_genes])
				printlist = []
				hg_directions = {}
				hg_lengths = defaultdict(list)
				for lt in curr_bgc_genes:
					ginfo = self.comp_gene_info[lt]
					hg = 'singleton'
					if lt in self.gene_to_hg:
						hg = self.gene_to_hg[lt]

					gstart = ginfo['start']
					gend = ginfo['end']
					forward = "FALSE"
					if ginfo['direction'] == '+': forward = "TRUE"

					hg_color = '"#dbdbdb"'
					if hg in self.hg_to_color:
						hg_color = '"' + self.hg_to_color[hg] + '"'

					gene_string = '\t'.join([str(x) for x in [bgc, lt, gstart, gend, forward, hg, hg_color]])
					printlist.append(gene_string)
					if hg != 'singleton':
						bgc_hg_presence[bgc][hg] = hg
						hg_counts[hg] += bgc_weights[bgc]
						hg_directions[hg] = ginfo['direction']
						hg_lengths[hg].append(gend - gstart)
				if i == 0:
					ref_hg_directions = hg_directions
					gggenes_track_handle.write('\n'.join(printlist) + '\n')
				else:
					flip_support = 0
					keep_support = 0
					for c in ref_hg_directions:
						if not c in hg_directions: continue
						hg_weight = statistics.mean(hg_lengths[c])
						if hg_directions[c] == ref_hg_directions[c]:
							keep_support += hg_weight
						else:
							flip_support += hg_weight

					# flip the genbank visual if necessary, first BGC processed is used as reference guide
					if flip_support > keep_support:
						flip_printlist = []
						for gene_string in printlist:
							gene_info = gene_string.split('\t')
							new_forward = 'TRUE'
							if gene_info[4] == 'TRUE': new_forward = 'FALSE'
							new_gstart = int(last_gene_end) - int(gene_info[3])
							new_gend = int(last_gene_end) - int(gene_info[2])
							new_gene_string = '\t'.join([str(x) for x in
														 [gene_info[0], gene_info[1], new_gstart, new_gend, new_forward,
															gene_info[-2], gene_info[-1]]])
							flip_printlist.append(new_gene_string)
						gggenes_track_handle.write('\n'.join(flip_printlist) + '\n')
					else:
						gggenes_track_handle.write('\n'.join(printlist) + '\n')

			dummy_hg = None
			for bgc in bgc_hg_presence:
				for hg in hg_counts:
					dummy_hg = hg
					heatmap_track_handle.write('\t'.join([bgc, hg, bgc_hg_presence[bgc][hg], str(hg_counts[hg])]) + '\n')

			for i, bgc in enumerate(all_bgcs_in_tree):
				if not bgc in bgc_gene_counts.keys():
					gggenes_track_handle.write('\t'.join([bgc] + ['NA']*4 + ['Absent', '"#FFFFFF"']) + '\n')
					heatmap_track_handle.write('\t'.join([bgc, dummy_hg, 'Absent', '1']) + '\n')
				elif i == 0:
					gggenes_track_handle.write('\t'.join([bgc] + ['NA']*4 + ['Absent', '"#FFFFFF"']) + '\n')

			gggenes_track_handle.close()
			heatmap_track_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error(
					"Had difficulties creating tracks for visualization of BGC gene architecture along phylogeny using R libraries.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

		rscript_plot_cmd = ["Rscript", RSCRIPT_FOR_BGSEE, phylogeny_file, gggenes_track_file, heatmap_track_file,
							result_pdf_file]
		if self.logObject:
			self.logObject.info('Running R-based plotting with the following command: %s' % ' '.join(rscript_plot_cmd))
		try:
			subprocess.call(' '.join(rscript_plot_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			self.logObject.info('Successfully ran: %s' % ' '.join(rscript_plot_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
				self.logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue running: %s' % ' '.join(rscript_plot_cmd))

		if self.logObject:
			self.logObject.info('Plotting completed (I think successfully)!')

	def constructCodonAlignments(self, outdir, cpus=1, only_scc=False, list_alignments=False, filter_outliers=False, use_magus=True):
		"""
		Function to automate construction of codon alignments. This function first extracts protein and nucleotide sequnces
		from BGC Genbanks, then creates protein alignments for each homolog group using MAFFT, and finally converts those
		into codon alignments using PAL2NAL.

		:param outdir: Path to output/workspace directory. Intermediate files (like extracted nucleotide and protein
						 sequences, protein and codon alignments, will be writen to respective subdirectories underneath this
						 one).
		:param cpus: Number of cpus/threads to use when fake-parallelizing jobs using multiprocessing.
		:param only_scc: Whether to construct codon alignments only for homolog groups which are found to be core and in
						 single copy for samples with the GCF. Note, if working with draft genomes and the BGC is fragmented
						 this should be able to still identify SCC homolog groups across the BGC instances belonging to the
						 GCF.
		"""

		nucl_seq_dir = os.path.abspath(outdir) + '/Nucleotide_Sequences/'
		prot_seq_dir = os.path.abspath(outdir) + '/Protein_Sequences/'
		prot_alg_dir = os.path.abspath(outdir) + '/Protein_Alignments/'
		codo_alg_dir = os.path.abspath(outdir) + '/Codon_Alignments/'
		if filter_outliers:
			nucl_seq_dir = os.path.abspath(outdir) + '/Nucleotide_Sequences_MAD_Refined/'
			prot_seq_dir = os.path.abspath(outdir) + '/Protein_Sequences_MAD_Refined/'
			prot_alg_dir = os.path.abspath(outdir) + '/Protein_Alignments_MAD_Refined/'
			codo_alg_dir = os.path.abspath(outdir) + '/Codon_Alignments_MAD_Refined/'

		if not os.path.isdir(nucl_seq_dir): os.system('mkdir %s' % nucl_seq_dir)
		if not os.path.isdir(prot_seq_dir): os.system('mkdir %s' % prot_seq_dir)
		if not os.path.isdir(prot_alg_dir): os.system('mkdir %s' % prot_alg_dir)
		if not os.path.isdir(codo_alg_dir): os.system('mkdir %s' % codo_alg_dir)

		pool_size = 1
		if cpus > 10:
			pool_size = math.floor(cpus / 10)
			cpus = 10

		all_samples = set(self.bgc_sample.values())
		try:
			inputs = []
			for hg in self.hg_genes:
				# if len(self.hg_genes[hg]) < 2: continue
				sample_counts = defaultdict(int)
				gene_sequences = {}
				for gene in self.hg_genes[hg]:
					gene_info = self.comp_gene_info[gene]
					bgc_id = gene_info['bgc_name']
					sample_id = self.bgc_sample[bgc_id]
					nucl_seq = gene_info['nucl_seq']
					prot_seq = gene_info['prot_seq']
					sample_counts[sample_id] += 1
					gid = sample_id + '|' + gene
					if only_scc:
						gid = sample_id
					gene_sequences[gid] = tuple([nucl_seq, prot_seq])
				samples_with_single_copy = set([s[0] for s in sample_counts.items() if s[1] == 1])
				# check that hg is single-copy-core
				if only_scc and len(samples_with_single_copy.symmetric_difference(all_samples)) > 0:
					continue
				elif only_scc and self.logObject:
					self.logObject.info('Homolog group %s detected as SCC across samples (not individual BGCs).' % hg)
				# check that hg is present in the original instances of GCF
				#if len([x for x in gene_sequences.keys() if len(x.split('|')[1].split('_')[0]) == 3]) == 0: continue
				if filter_outliers:
					gene_sequences = util.determineOutliersByGeneLength(gene_sequences, self.logObject)
				inputs.append([hg, gene_sequences, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir, cpus, use_magus, self.logObject])

			p = multiprocessing.Pool(pool_size)
			p.map(create_codon_msas, inputs)

			if not filter_outliers:
				self.nucl_seq_dir = nucl_seq_dir
				self.prot_seq_dir = prot_seq_dir
				self.prot_alg_dir = prot_alg_dir
				self.codo_alg_dir = codo_alg_dir
			else:
				self.nucl_filt_seq_dir = nucl_seq_dir
				self.prot_filt_seq_dir = prot_seq_dir
				self.prot_filt_alg_dir = prot_alg_dir
				self.codo_filt_alg_dir = codo_alg_dir

			if list_alignments:
				codon_alg_listings_file = outdir + 'Codon_Alignments_Listings.txt'
				if filter_outliers:
					codon_alg_listings_file = outdir + 'Codon_Alignments_Listings.MAD_Refined.txt'
				codon_alg_listings_handle = open(codon_alg_listings_file, 'w')
				for f in os.listdir(codo_alg_dir):
					codon_alg_listings_handle.write(f.split('.msa.fna')[0] + '\t' + codo_alg_dir + f + '\n')
				codon_alg_listings_handle.close()

		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues with create protein/codon alignments of SCC homologs for BGC.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def constructGCFPhylogeny(self, output_alignment, output_phylogeny, only_scc=False, ambiguious_position_cutoff=0.000001):
		"""
		Function to create phylogeny based on codon alignments of SCC homolog groups for GCF.

		:param output_alignment: Path to output file for concatenated SCC homolog group alignment.
		:param output_phylogeny: Path to output file for approximate maximum-likelihood phylogeny produced by FastTree2 from
									 concatenated SCC homolog group alignment.
		"""
		try:
			if only_scc:
				bgc_sccs = defaultdict(lambda: "")
				fasta_data = []
				fasta_data_tr = []

				for f in os.listdir(self.codo_alg_dir):
					hg_align_msa = self.codo_alg_dir + f
					# concatenate gene alignments
					with open(hg_align_msa) as opm:
						for rec in SeqIO.parse(opm, 'fasta'):
							bgc_sccs['>' + rec.id] += str(rec.seq).upper()

				for b in bgc_sccs:
					fasta_data.append([b] + list(bgc_sccs[b]))

				for i, ls in enumerate(zip(*fasta_data)):
					if i == 0:
						fasta_data_tr.append(ls)
					else:
						n_count = len([x for x in ls if x == '-'])
						if (float(n_count) / len(ls)) < 0.1:
							fasta_data_tr.append(list(ls))

				scc_handle = open(output_alignment, 'w')

				for rec in zip(*fasta_data_tr):
					scc_handle.write(rec[0] + '\n' + ''.join(rec[1:]) + '\n')
				scc_handle.close()
			else:
				bgc_sccs = defaultdict(lambda: "")
				fasta_data = []
				fasta_data_tr = []

				for f in os.listdir(self.codo_alg_dir):
					hg_align_msa = self.codo_alg_dir + f
					#print(f)
					# perform consensus calling
					sample_seqs = defaultdict(list)
					with open(hg_align_msa) as opm:
						for rec in SeqIO.parse(opm, 'fasta'):
							sample = rec.id.split('|')[0]
							sample_seqs[sample].append(list(str(rec.seq).upper()))

					for samp in sample_seqs:
						samp_seqs = sample_seqs[samp]
						consensus_seq = []
						for alleles in zip(*samp_seqs):
							valid_alleles = set([a for a in list(alleles) if a in set(['A', 'C', 'G', 'T'])])
							if len(valid_alleles) == 1:
								consensus_seq.append(list(valid_alleles)[0])
							else:
								consensus_seq.append('-')
						bgc_sccs['>' + samp] += "".join(consensus_seq)

				for b in bgc_sccs:
					fasta_data.append([b] + list(bgc_sccs[b]))

				for i, ls in enumerate(zip(*fasta_data)):
					if i == 0:
						fasta_data_tr.append(ls)
					else:
						n_count = len([x for x in ls if x == '-'])
						if (float(n_count) / len(ls)) < ambiguious_position_cutoff:
							fasta_data_tr.append(list(ls))

				scc_handle = open(output_alignment, 'w')

				for rec in zip(*fasta_data_tr):
					scc_handle.write(rec[0] + '\n' + ''.join(rec[1:]) + '\n')
				scc_handle.close()

		except Exception as e:
			if self.logObject:
				self.logObject.error('Had issues with creating concatenated alignment of the SCC homolog groups.')
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

		# use FastTree2 to construct phylogeny
		fasttree_cmd = ['fasttree', '-nt', output_alignment, '>', output_phylogeny]
		if self.logObject:
			self.logObject.info('Running FastTree2 with the following command: %s' % ' '.join(fasttree_cmd))
		try:
			subprocess.call(' '.join(fasttree_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(fasttree_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(fasttree_cmd))
				self.logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue running: %s' % ' '.join(fasttree_cmd))

	def refineBGCGenbanks(self, new_gcf_listing_file, outdir, first_boundary_homolog, second_boundary_homolog):
		"""
		Function to refine BGC Genbanks based on boundaries defined by two single copy core homolog groups. Genbanks
		are filtered to retain only features in between the positions of the two boundary homolog groups. Coordinates
		relevant to the BGC framework are updated (but not all location coordinates!!!).

		:param new_gcf_listing_file: Path to where new GCF listing file will be written.
		:param outdir: Path to workspace directory.
		:param first_boundary_homolog: Identifier of the first boundary homolog group
		:param second_boundary_homolog: Identifier of the second boundary homolog group
		"""
		try:
			refined_gbks_dir = outdir + 'Refined_Genbanks/'
			if not os.path.isdir(refined_gbks_dir): os.system('mkdir %s' % refined_gbks_dir)

			nglf_handle = open(new_gcf_listing_file, 'w')

			first_boundary_homolog_genes = self.hg_genes[first_boundary_homolog]
			second_boundary_homolog_genes = self.hg_genes[second_boundary_homolog]

			for bgc in self.pan_bgcs:
				bgc_genes = set(self.pan_bgcs[bgc].gene_information.keys())
				bgc_fbh_genes = bgc_genes.intersection(first_boundary_homolog_genes)
				bgc_sbh_genes = bgc_genes.intersection(second_boundary_homolog_genes)
				if len(bgc_fbh_genes) == 1 and len(bgc_sbh_genes) == 1:
					refined_gbk = refined_gbks_dir + self.bgc_gbk[bgc].split('/')[-1]
					self.pan_bgcs[bgc].refineGenbank(refined_gbk, list(bgc_fbh_genes)[0], list(bgc_sbh_genes)[0])
					nglf_handle.write('%s\t%s\n' % (self.bgc_sample[bgc], refined_gbk))
				elif self.logObject:
					self.logObject.warning(
						"Dropping the BGC genbank %s from consideration / refinement process because it does not have the boundary homolog groups in single-copy copy." %
						self.bgc_gbk[bgc])
			nglf_handle.close()

		except:
			if self.logObject:
				self.logObject.error('Had an issue refining BGC genbanks associated with GCF')
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def determineHgOrderIndex(self):
		"""
		Function to determine an "ordering" score for homolog groups in GCF. The order score is relative to each GCF,
		even a homolog group has a large or small order score indicates it is on the edges (beginning will be chosen
		arbitrarily).

		To determine this, a Markov chain esque approach is used, whereby gene order is determined by where in the chain
		a homolog group is best positioned.
		"""
		try:
			bgc_gene_counts = defaultdict(int)
			core_bgcs = set([])
			for bgc in self.bgc_genes:
				bgc_gene_counts[bgc] = len(self.bgc_genes[bgc])
				if len(list(self.bgc_genes[bgc])[0].split('_')[0]) == 3:
					core_bgcs.add(bgc)
			ref_bgc = None
			for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
				if item[0] in core_bgcs:
					ref_bgc = item[0]
					break
			if ref_bgc == None:
				for i, item in enumerate(sorted(bgc_gene_counts.items(), key=itemgetter(1), reverse=True)):
					ref_bgc = item[0]
					break

			bgcs_ref_first = [ref_bgc] + sorted(list(set(self.bgc_genes.keys()).difference(set([ref_bgc]))))
			ref_hg_directions = {}
			hg_pair_scpus = defaultdict(int)
			hg_preceding_scpus = defaultdict(lambda: defaultdict(int))
			hg_following_scpus = defaultdict(lambda: defaultdict(int))
			all_hgs = set(['start', 'end'])
			direction_forward_support = defaultdict(int)
			direction_reverse_support = defaultdict(int)
			for i, bgc in enumerate(bgcs_ref_first):
				curr_bgc_genes = self.bgc_genes[bgc]
				hg_directions = {}
				hg_lengths = defaultdict(list)
				hg_starts = {}
				for g in sorted(curr_bgc_genes):
					ginfo = self.comp_gene_info[g]
					gstart = ginfo['start']
					gend = ginfo['end']
					if g in self.gene_to_hg:
						hg = self.gene_to_hg[g]
						hg_directions[hg] = ginfo['direction']
						hg_lengths[hg].append(abs(gend - gstart))
						hg_starts[hg] = ginfo['start']

				reverse_flag = False
				if i == 0:
					ref_hg_directions = hg_directions
				else:
					flip_support = 0
					keep_support = 0
					for c in ref_hg_directions:
						if not c in hg_directions: continue
						hg_weight = statistics.mean(hg_lengths[c])
						if hg_directions[c] == ref_hg_directions[c]:
							keep_support += hg_weight
						else:
							flip_support += hg_weight

					# reverse ordering
					if flip_support > keep_support:
						reverse_flag = True

				hgs = []
				for c in sorted(hg_starts.items(), key=itemgetter(1), reverse=reverse_flag):
					hgs.append(c[0])
					if reverse_flag == False:
						if hg_directions[c[0]] == '+':
							direction_forward_support[c[0]] += 1
						elif hg_directions[c[0]] == '-':
							direction_reverse_support[c[0]] += 1
					else:
						if hg_directions[c[0]] == '+':
							direction_reverse_support[c[0]] += 1
						elif hg_directions[c[0]] == '-':
							direction_forward_support[c[0]] += 1

				for j, hg in enumerate(hgs):
					all_hgs.add(hg)
					if j == 0:
						hg_previ = "start"
						hg_preceding_scpus[hg][hg_previ] += 1
						hg_following_scpus[hg_previ][hg] += 1
						hg_pair_scpus[tuple([hg_previ, hg])] += 1
					try:
						hg_after = hgs[j+1]
						# make sure you don't get lost with broken/fragmented genes in BGCs that might be
						# in the process being lost.
						if hg != hg_after:
							hg_preceding_scpus[hg_after][hg] += 1
							hg_following_scpus[hg][hg_after] += 1
							hg_pair_scpus[tuple([hg, hg_after])] += 1
					except:
						hg_after = 'end'
						hg_preceding_scpus[hg_after][hg] += 1
						hg_following_scpus[hg][hg_after] += 1
						hg_pair_scpus[tuple([hg, hg_after])] += 1

			anchor_edge = None
			for hps in sorted(hg_pair_scpus.items(), key=itemgetter(1), reverse=True):
				if hps[0][0] in self.protocluster_core_homologs or hps[0][1] in self.protocluster_core_homologs:
					anchor_edge = hps[0]
					break
			try:
				assert(anchor_edge != None)
			except:
				raise RuntimeError("Unexpected error, no anchor edge found, could be because no protocore homolog group exists, which shouldn't be the case!")
				sys.exit(1)

			# use to keep track of which HGs have been accounted for already at different steps of assigning order
			accounted_hgs = set([anchor_edge[0], anchor_edge[1]])

			# primary expansion left
			curr_hg = anchor_edge[0]
			left_expansion = [curr_hg]
			while not curr_hg == 'start':
				new_hg = None
				for i, hg in enumerate(sorted(hg_preceding_scpus[curr_hg].items(), key=itemgetter(1), reverse=True)):
					if not hg[0] in accounted_hgs:
						new_hg = hg[0]
						left_expansion = [new_hg] + left_expansion
						accounted_hgs.add(new_hg)
						break
				if new_hg != None:
					curr_hg = new_hg
				else:
					# shouldn't ever be the case, but breaking just in case
					break

			# primary expansion right
			curr_hg = anchor_edge[1]
			right_expansion = [curr_hg]
			while not curr_hg == 'end':
				new_hg = None
				for i, hg in enumerate(sorted(hg_following_scpus[curr_hg].items(), key=itemgetter(1), reverse=True)):
					if not hg[0] in accounted_hgs:
						new_hg = hg[0]
						right_expansion.append(new_hg)
						accounted_hgs.add(new_hg)
						break
				if new_hg != None:
					curr_hg = new_hg
				else:
					# shouldn't ever be the case, but breaking just in case
					break

			primary_path_ordered = left_expansion + right_expansion
			ordered_hgs_list = primary_path_ordered

			# figure out where non-accounted for HGs belong best in the primary path.
			not_accounted_hgs = all_hgs.difference(accounted_hgs)
			while len(not_accounted_hgs) > 0:
				progress_made = False
				for hg in sorted(not_accounted_hgs):
					best_score = 0
					relative_pos = None
					neighboriest_hg = None
					for phg in sorted(hg_preceding_scpus[hg].items(), key=itemgetter(1), reverse=True):
						if best_score < phg[1] and phg[0] in accounted_hgs:
							best_score = phg[1]
							relative_pos = 'after'
							neighboriest_hg = phg[0]
							break
					for fhg in sorted(hg_following_scpus[hg].items(), key=itemgetter(1), reverse=True):
						if best_score < fhg[1] and fhg[0] in accounted_hgs:
							best_score = fhg[1]
							relative_pos = 'before'
							neighboriest_hg = fhg[0]
							break
					if best_score > 0:
						neighboriest_hg_index = ordered_hgs_list.index(neighboriest_hg)
						#print(hg + '\t' + str(best_score) + '\t'+ relative_pos + '\t' + str(neighboriest_hg) + '\t' + str(neighboriest_hg_index))

						if relative_pos == 'before':
							ordered_hgs_list.insert(neighboriest_hg_index, hg)
						elif relative_pos == 'after':
							ordered_hgs_list.insert(neighboriest_hg_index+1, hg)
						accounted_hgs.add(hg)
						not_accounted_hgs = all_hgs.difference(accounted_hgs)
						progress_made = True
						break

				if not progress_made:
					break
			# these shouldn't really exist but just append them to the end if they do
			unaccountable_hgs = all_hgs.difference(accounted_hgs)
			ordered_hgs_list += list(sorted(unaccountable_hgs))

			i = 1
			for hg in ordered_hgs_list:
				if not hg in set(['start', 'end']):
					consensus_direction = '0'
					if direction_forward_support[hg] >= direction_reverse_support[hg]: consensus_direction = '1'
					self.hg_order_scores[hg] = [i, consensus_direction]
					i+=1

		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues in attempting to calculate order score for each homolog group.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def runPopulationGeneticsAnalysis(self, outdir, cpus=1, population=None, filter_outliers=False, population_analysis_on=False, gw_pairwise_similarities=None, use_translation=True, species_phylogeny=None, sample_size=1000):
		"""
		Wrapper function which serves to parallelize population genetics analysis.

		:param outdir: The path to the workspace / output directory.
		:param cpus: The number of cpus (will be used for parallelizing)
		"""

		popgen_dir = outdir + 'Codon_PopGen_Analyses'
		plots_dir = outdir + 'Codon_MSA_Plots'
		final_output_file = outdir + 'Homolog_Group_Information'

		if filter_outliers:
			final_output_file = outdir + 'Homolog_Group_Information_MAD_Refined'
			popgen_dir = outdir + 'Codon_PopGen_Analyses_MAD_Refined'
			plots_dir = outdir + 'Codon_MSA_Plots_MAD_Refined'

		if population:
			final_output_file = final_output_file + '_Pop-' + str(population) + '.txt'
			popgen_dir += '_Pop-' + str(population) + '/'
			plots_dir += '_Pop-' + str(population) + '/'
		else:
			final_output_file = final_output_file + '.txt'
			popgen_dir += '/'
			plots_dir += '/'

		if not os.path.isdir(popgen_dir): os.system('mkdir %s' % popgen_dir)
		if not os.path.isdir(plots_dir): os.system('mkdir %s' % plots_dir)

		final_output_handle = open(final_output_file, 'w')
		header = ['gcf_id', 'gcf_annotation', 'homolog_group', 'annotation', 'hg_order_index', 'hg_consensus_direction',
				  'hg_proportion_multi-copy_genome-wide', 'single-copy_in_GCF_context', 'median_gene_length', 'is_core_to_bgc',
				  'num_of_hg_instances', 'samples_with_hg', 'proportion_of_samples_with_hg',
				  'ambiguous_sites_proporition', 'Tajimas_D', 'proportion_variable_sites',
				  'proportion_nondominant_major_allele', 'median_beta_rd', 'max_beta_rd', 'median_dn_ds', 'mad_dn_ds',
				  'all_domains']
		if population:
			header = ['population'] + header
		elif population_analysis_on:
			header = header[:-1]
			header += ['populations_with_hg', 'proportion_of_total_populations_with_hg',
					   'most_significant_Fisher_exact_pvalues_presence_absence', 'median_Tajimas_D_per_population',
					   'mad_Tajimas_D_per_population', 'most_negative_population_Tajimas_D',
					   'most_positive_population_Tajimas_D', 'population_entropy', 'median_fst_like_estimate',
					   'population_proportion_of_members_with_hg', 'all_domains']

		final_output_handle.write('\t'.join(header) + '\n')

		inputs = []
		input_codon_dir = self.codo_alg_dir
		if filter_outliers:
			input_codon_dir = self.codo_filt_alg_dir

		products = defaultdict(float)
		for bgc in self.bgc_product:
			for prod in self.bgc_product[bgc]:
				products[prod] += 1.0 / len(self.bgc_product[bgc])
		gcf_product_summary = '; '.join([x[0] + ':' + str(x[1]) for x in products.items()])

		for f in os.listdir(input_codon_dir):
			hg = f.split('.msa.fna')[0]
			codon_alignment_fasta = input_codon_dir + f
			if gw_pairwise_similarities != None:
				gw_pairwise_similarities = dict(gw_pairwise_similarities)
			sample_population_local = self.sample_population
			if sample_population_local != None:
				sample_population_local = dict(sample_population_local)
			inputs.append([self.gcf_id, gcf_product_summary, hg, codon_alignment_fasta, popgen_dir, plots_dir, self.comp_gene_info,
						   self.hg_genes, self.bgc_sample, self.hg_prop_multi_copy, dict(self.hg_order_scores),
						   gw_pairwise_similarities, use_translation, sample_population_local, population,
						   species_phylogeny, sample_size, self.logObject])

		p = multiprocessing.Pool(cpus)
		p.map(popgen_analysis_of_hg, inputs)

		final_output_handle = open(final_output_file, 'a+')
		data = []
		data_for_sorting = []
		pvals = []
		for f in os.listdir(popgen_dir):
			if not f.endswith('_stats.txt'): continue
			with open(popgen_dir + f) as opf:
				for line in opf:
					line = line.rstrip('\n')
					ls = line.split('\t')
					if population_analysis_on:
						pvals.append(float(ls[-9]))
						data.append(ls)
					else:
						con_order = 1E10
						if util.is_numeric(ls[4]):
							con_order = int(ls[4])
						data_for_sorting.append([con_order, line])

		adj_pvals = util.p_adjust_bh(pvals)

		for i, ls in enumerate(data):
			newline = '\t'.join(ls[:-9]) + '\t' + str(adj_pvals[i]) + '\t' + '\t'.join(ls[-8:])
			con_order = 1E10
			if util.is_numeric(ls[4]):
				con_order = int(ls[4])
			data_for_sorting.append([con_order, newline])

		for ls in sorted(data_for_sorting, key=itemgetter(0)):
			final_output_handle.write(ls[1] + '\n')
		final_output_handle.close()

	def identifyGCFInstances(self, outdir, sample_prokka_data, orthofinder_matrix_file, min_size=5, min_core_size=3,
							 gcf_to_gcf_transition_prob=0.9, background_to_background_transition_prob=0.9,
							 syntenic_correlation_threshold=0.8, surround_gene_max=5, no_orthogroup_matrix=False,
							 bgc_prediction_method='ANTISMASH', loose_flag=False, cpus=1, block_size=3000):
		"""
		Function to search for instances of GCF in sample using HMM based approach based on homolog groups as characters,
		"part of GCF" and "not part of GCF" as states - all trained on initial BGCs constituting GCF as identified by
		zol-Cluster.py. This function utilizes the convenient Python library Pomegranate.

		:param outdir: The path to the workspace / output directory.
		:param sample_prokka_data: Dictionary with keys being sample identifiers and values being paths to genbank and
									 proteome files from Prokka based annotation
		:param cpus: The number of cpus (will be used for parallelizing)
		"""

		bgc_genbanks_dir = os.path.abspath(outdir + 'BGC_Genbanks') + '/'
		bgc_info_dir = os.path.abspath(outdir + 'BGC_Sample_Info') + '/'
		if not os.path.isdir(bgc_genbanks_dir): os.system('mkdir %s' % bgc_genbanks_dir)
		if not os.path.isdir(bgc_info_dir): os.system('mkdir %s' % bgc_info_dir)

		# Estimate HMM parameters
		gcf_hg_probabilities = defaultdict(lambda: 0.0)
		other_hg_probabilities = defaultdict(lambda: 0.0)
		number_gcf_hgs = 0
		number_other_hgs = 0
		specific_hgs = set([])
		with open(orthofinder_matrix_file) as ofmf:
			for i, line in enumerate(ofmf):
				if i == 0: continue
				line = line.strip('\n')
				ls = line.split('\t')
				hg = ls[0]
				if hg in self.hg_genes.keys():
					number_gcf_hgs += 1
					total_genes = 0
					gcf_genes = 0
					for samp_genes in ls[1:]:
						for gene in samp_genes.split(', '):
							if not gene.strip(): continue
							total_genes += 1
							if gene in self.pan_genes:
								gcf_genes += 1

					if float(gcf_genes) == float(total_genes) and self.hg_max_self_evalue[hg][1] == True:
						specific_hgs.add(hg)

					other_hg_probabilities[hg] = 0.01
					gcf_hg_probabilities[hg] = 0.99

					if self.hg_max_self_evalue[hg][1] == False:
						other_hg_probabilities[hg] = min(1.0 - (gcf_genes / float(total_genes)), 0.2)
						gcf_hg_probabilities[hg] = 1.0 - other_hg_probabilities[hg]
				else:
					number_other_hgs += 1

		gcf_hg_probabilities['other'] = 0.2
		other_hg_probabilities['other'] = 0.8

		gcf_distribution = DiscreteDistribution(dict(gcf_hg_probabilities))
		other_distribution = DiscreteDistribution(dict(other_hg_probabilities))

		gcf_state = State(gcf_distribution, name='GCF')
		other_state = State(other_distribution, name='Non GCF')

		model = HiddenMarkovModel()
		model.add_states(gcf_state, other_state)

		# estimate transition probabilities
		gcf_to_gcf = gcf_to_gcf_transition_prob  # float(number_gcf_hgs - 1) / float(number_gcf_hgs)
		gcf_to_other = 1.0 - gcf_to_gcf_transition_prob  # 1.0 - gcf_to_gcf
		other_to_other = background_to_background_transition_prob  # float(number_other_hgs - 1) / float(number_other_hgs)
		other_to_gcf = 1.0 - background_to_background_transition_prob  # 1.0 - other_to_other

		start_to_gcf = 0.5	# float(number_gcf_hgs)/float(number_gcf_hgs + number_other_hgs)
		start_to_other = 0.5  # 1.0 - start_to_gcf
		gcf_to_end = 0.5  # float(number_gcf_hgs)/float(number_gcf_hgs + number_other_hgs)
		other_to_end = 0.5	# 1.0 - gcf_to_end

		model.add_transition(model.start, gcf_state, start_to_gcf)
		model.add_transition(model.start, other_state, start_to_other)
		model.add_transition(gcf_state, model.end, gcf_to_end)
		model.add_transition(other_state, model.end, other_to_end)
		model.add_transition(gcf_state, gcf_state, gcf_to_gcf)
		model.add_transition(gcf_state, other_state, gcf_to_other)
		model.add_transition(other_state, gcf_state, other_to_gcf)
		model.add_transition(other_state, other_state, other_to_other)

		model.bake()

		bgc_hmm_evalues_file = outdir + 'GCF_NewInstances_HMMEvalues.txt'
		expanded_gcf_list_file = outdir + 'GCF_Expanded.txt'

		# open handle to file where expanded GCF listings will be written
		expanded_gcf_list_handle = open(expanded_gcf_list_file, 'w')
		with open(self.bgc_genbanks_listing) as obglf:
			for line in obglf:
				expanded_gcf_list_handle.write(line)
		expanded_gcf_list_handle.close()

		total_samples = sorted(sample_prokka_data.keys())
		for block_samp_list in util.chunks(total_samples, block_size):
			block_samp_set = set(block_samp_list)

			# start process of finding new BGCs
			sample_lt_to_hg = defaultdict(dict)
			sample_hgs = defaultdict(set)
			sample_lt_to_evalue = defaultdict(dict)
			for lt in self.hmmscan_results:
				for i, hits in enumerate(sorted(self.hmmscan_results[lt], key=itemgetter(1))):
					if i == 0 and hits[2] in block_samp_set:
						sample_lt_to_hg[hits[2]][lt] = hits[0]
						sample_hgs[hits[2]].add(hits[0])
						sample_lt_to_evalue[hits[2]][lt] = decimal.Decimal(hits[1])

			simplified_comp_gene_info = defaultdict(dict)
			for g in self.comp_gene_info:
				simplified_comp_gene_info[g]['start'] = self.comp_gene_info[g]['start']
				simplified_comp_gene_info[g]['end'] = self.comp_gene_info[g]['end']
				simplified_comp_gene_info[g]['direction'] = self.comp_gene_info[g]['direction']

			identify_gcf_segments_input = []
			for sample in sample_hgs:
				if len(sample_hgs[sample]) < 3: continue
				hgs_ordered_dict = {}
				lts_ordered_dict = {}
				for scaffold in self.scaffold_genes[sample]:
					lts_with_start = []
					for lt in self.scaffold_genes[sample][scaffold]:
						lts_with_start.append([lt, self.gene_location[sample][lt]['start']])

					hgs_ordered = []
					lts_ordered = []
					for lt_start in sorted(lts_with_start, key=itemgetter(1)):
						lt, start = lt_start
						lts_ordered.append(lt)
						if lt in sample_lt_to_hg[sample]:
							hgs_ordered.append(sample_lt_to_hg[sample][lt])
						else:
							hgs_ordered.append('other')
					hgs_ordered_dict[scaffold] = hgs_ordered
					lts_ordered_dict[scaffold] = lts_ordered

				identify_gcf_segments_input.append([bgc_info_dir, bgc_genbanks_dir, sample, sample_prokka_data[sample],
													sample_lt_to_evalue[sample], dict(self.hmmscan_results_lenient[sample]),
													model, lts_ordered_dict, hgs_ordered_dict, dict(simplified_comp_gene_info),
													dict(self.gene_location[sample]), dict(self.gene_id_to_order[sample]),
													dict(self.gene_order_to_id[sample]), self.protocluster_core_homologs,
													self.core_homologs, self.boundary_genes[sample], specific_hgs,
													dict(self.bgc_genes), dict(self.gene_to_hg), min_size, min_core_size,
													surround_gene_max, syntenic_correlation_threshold, loose_flag])
			with multiprocessing.Manager() as manager:
				with manager.Pool(cpus) as pool:
					pool.map(identify_gcf_instances, identify_gcf_segments_input)

		os.system('find %s -type f -name "*.bgcs.txt" -exec cat {} + >> %s' % (bgc_info_dir, expanded_gcf_list_file))

		os.system('find %s -type f -name "*.hg_evalues.txt" -exec cat {} + >> %s' % (bgc_info_dir, bgc_hmm_evalues_file))
		if not os.path.isfile(bgc_hmm_evalues_file): os.system('touch %s' % bgc_hmm_evalues_file)

		if not no_orthogroup_matrix:
			bgc_lt_to_hg = defaultdict(lambda: defaultdict(dict))
			with open(bgc_hmm_evalues_file) as oghef:
				for line in oghef:
					line = line.strip()
					bgc_gbk_path, sample, lt, hg, eval, hg_is_functionally_core = line.split('\t')
					bgc_lt_to_hg[bgc_gbk_path][lt] = hg

			sample_hg_proteins = defaultdict(lambda: defaultdict(set))
			all_samples = set([])
			with open(expanded_gcf_list_file) as obhef:
				for line in obhef:
					line = line.strip()
					sample, bgc_gbk_path = line.split('\t')
					is_expansion_flag = False
					if '_Expansion_BGC_' in bgc_gbk_path:
						is_expansion_flag = True
					BGC_Object = BGC(bgc_gbk_path, bgc_gbk_path, is_expansion_flag, prediction_method=bgc_prediction_method)
					BGC_Object.parseGenbanks(comprehensive_parsing=False)
					curr_bgc_lts = set(BGC_Object.gene_information.keys())
					sample = util.cleanUpSampleName(sample)
					for lt in curr_bgc_lts:
						if lt in bgc_lt_to_hg[bgc_gbk_path]:
							hg = bgc_lt_to_hg[bgc_gbk_path][lt]
							sample_hg_proteins[sample][hg].add(lt)
					all_samples.add(sample)

			original_samples = []
			all_hgs = set([])
			with open(orthofinder_matrix_file) as omf:
				for i, line in enumerate(omf):
					line = line.strip('\n')
					ls = line.split('\t')
					if i == 0:
						original_samples = [util.cleanUpSampleName(x) for x in ls[1:]]
						all_samples = all_samples.union(set(original_samples))
					else:
						hg = ls[0]
						all_hgs.add(hg)
						for j, prot in enumerate(ls[1:]):
							sample_hg_proteins[original_samples[j]][hg] = sample_hg_proteins[original_samples[j]][hg].union(set(prot.split(', ')))

			expanded_orthofinder_matrix_file = outdir + 'Orthogroups.expanded.tsv'
			expanded_orthofinder_matrix_handle = open(expanded_orthofinder_matrix_file, 'w')

			header = [''] + [s for s in sorted(all_samples)]
			expanded_orthofinder_matrix_handle.write('\t'.join(header) + '\n')
			for hg in sorted(all_hgs):
				printlist = [hg]
				for s in sorted(all_samples):
					printlist.append(', '.join(sample_hg_proteins[s][hg]))
				expanded_orthofinder_matrix_handle.write('\t'.join(printlist) + '\n')
			expanded_orthofinder_matrix_handle.close()

	def extractGenesAndCluster(self, genes_representative_fasta, genes_fasta, codon_alignments_file, bowtie2_db_prefix):
		"""
		Function to cluster gene sequences for each homolog group using codon alignments and then select representative
		sequences which will be used to build a Bowtie2 reference database.

		:param genes_representative_fasta: Path to FASTA file which will harbor gene + flanks
		:param codon_alignments_file: File listing paths to codon alignments (2nd column) for each homolog group (1st
									  column).
		:param bowtei2_db_prefix: Path to prefix of Bowtie2 refernece database/index to be used for aligning downstream
															in the zol-DiscoVary workflow
		"""

		try:
			grf_handle = open(genes_representative_fasta, 'w')
			gf_handle = open(genes_fasta, 'w')
			with open(codon_alignments_file) as ocaf:
				for line in ocaf:
					line = line.strip()
					hg, cod_alignment_fasta = line.split('\t')
					alleles_clustered, pair_matching = util.determineAllelesFromCodonAlignment(cod_alignment_fasta)
					all_genes_in_codon_alignments = set([])
					for ac in alleles_clustered:
						for gid in alleles_clustered[ac]:
							all_genes_in_codon_alignments.add(gid.split('|')[-1])
					# NOTE, the following is a warning because only genes which fall within one median absolute deviation
					# from the median length of genes for a homolog group are considered eligible as serving for
					# representative alleles. This is the default, but will consider making this a user option in the
					# future! - it was also over-reporting warnings which has been corrected as of Aug 16, 2022.
					try:
						assert(len(self.hg_genes[hg].symmetric_difference(all_genes_in_codon_alignments)) == 0)
					except:
						self.logObject.warning("Not all genes featured in codon alignment for homolog group %s, these will be excluded." % hg)

					for allele_cluster in alleles_clustered:
						best_rep_score = defaultdict(float)
						for ag1 in alleles_clustered[allele_cluster]:
							gf_handle.write('>' + hg + '|' + allele_cluster + '|' + ag1 + '\n' + str(self.comp_gene_info[ag1.split('|')[-1]]['nucl_seq']) + '\n')
							for ag2 in alleles_clustered[allele_cluster]:
								best_rep_score[ag1] += pair_matching[ag1][ag2]
						representative_gene = sorted([x for x in best_rep_score if best_rep_score[x] == max(best_rep_score.values())])[0]
						for ag in alleles_clustered[allele_cluster]:
							self.instance_to_haplotype[hg + '|' + allele_cluster + '|' + ag] = hg + '|' + allele_cluster + '|' + representative_gene
						grf_handle.write('>' + hg + '|' + allele_cluster + '|' + representative_gene + '\n' + str(self.comp_gene_info[representative_gene.split('|')[-1]]['nucl_seq']) + '\n')
			grf_handle.close()
			gf_handle.close()

			bowtie2_build = ['bowtie2-build', genes_representative_fasta, bowtie2_db_prefix]
			if self.logObject:
				self.logObject.info('Running the following command: %s' % ' '.join(bowtie2_build))
			try:
				subprocess.call(' '.join(bowtie2_build), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
								executable='/bin/bash')
				if self.logObject:
					self.logObject.info('Successfully ran: %s' % ' '.join(bowtie2_build))
			except:
				if self.logObject:
					self.logObject.error('Had an issue running: %s' % ' '.join(bowtie2_build))
				raise RuntimeError('Had an issue running: %s' % ' '.join(bowtie2_build))
			if self.logObject:
				self.logObject.info('Build Bowtie2 database/index for %s' % genes_representative_fasta)

		except Exception as e:
			if self.logObject:
				self.logObject.error("Unable to create FASTA of genes and/or subsequent bowtie2 database.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def extractGeneWithFlanksAndCluster(self, genes_with_flanks_fasta, cd_hit_clusters_fasta_file, cd_hit_nr_fasta_file, bowtie2_db_prefix):
		"""
		Function to extract gene sequences and surrounding flanking sequences into a FASTA file, which will then be
		clustered using CD-HIT at both a coarse and very granular (just remove 100% redundancy) and to construct
		a Bowtie2 reference database of the granular clustering. From the coarse clustering, representative genes
		will be selected to depict different alleles of the gene and stored as a dictionary.

		:param genes_with_flanks_fasta: Path to FASTA file which will harbor gene + flanks
		:param cd_hit_clusters_fasta_file: Path to FASTA file which will be used to output coarse level clustering by
											 CD-HIT
		:param cd_hit_nr_fasta_file: Path to FASTA file which will be used to output granular level clustering by CD-HIT
		:param bowtei2_db_prefix: Path to prefix of Bowtie2 refernece database/index to be used for aligning downstream
															in the zol-DiscoVary workflow
		"""
		try:
			gwff_handle = open(genes_with_flanks_fasta, 'w')
			for bgc in self.bgc_genes:
				for gene in self.bgc_genes[bgc]:
					if gene in self.gene_to_hg:
						gwff_handle.write('>' + gene + '|' + bgc + '|' + self.gene_to_hg[gene] + '\n' + self.comp_gene_info[gene]['nucl_seq_with_flanks'] + '\n')
			gwff_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error("Unable to extract flanking sequences of gene into FASTA file.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

		cd_hit_nr = ['cd-hit-est', '-i', genes_with_flanks_fasta, '-o', cd_hit_nr_fasta_file, '-G', '1', '-g',
					 '1', '-d', '0', '-n', '10', '-M', '2000', '-c', '1.0', '-aL', '1.0', '-aS', '1.0', '-T', '1']
		if self.logObject:
			self.logObject.info('Running the following command: %s' % ' '.join(cd_hit_nr))
		try:
			subprocess.call(' '.join(cd_hit_nr), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(cd_hit_nr))
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(cd_hit_nr))
			raise RuntimeError('Had an issue running: %s' % ' '.join(cd_hit_nr))
		if self.logObject:
			self.logObject.info('Ran CD-HIT for collapsing redundancy.')

		cd_hit_cluster = ['cd-hit-est', '-i', genes_with_flanks_fasta, '-o', cd_hit_clusters_fasta_file, '-G', '1',
							'-g',
							'1', '-d', '0', '-n', '10', '-M', '2000', '-c', '0.98', '-aL', '0.95', '-aS', '0.95', '-T',
							'1']
		if self.logObject:
			self.logObject.info('Running the following command: %s' % ' '.join(cd_hit_cluster))
		try:
			subprocess.call(' '.join(cd_hit_cluster), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(cd_hit_cluster))
		except:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(cd_hit_cluster))
			raise RuntimeError('Had an issue running: %s' % ' '.join(cd_hit_cluster))
		if self.logObject:
			self.logObject.info('Ran CD-HIT for clustering genes, with their flanks, into haplotype groups.')

		bowtie2_build = ['bowtie2-build', cd_hit_nr_fasta_file, bowtie2_db_prefix]
		if self.logObject:
			self.logObject.info('Running the following command: %s' % ' '.join(bowtie2_build))
		try:
			subprocess.call(' '.join(bowtie2_build), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(bowtie2_build))
		except:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(bowtie2_build))
			raise RuntimeError('Had an issue running: %s' % ' '.join(bowtie2_build))
		if self.logObject:
			self.logObject.info('Build Bowtie2 database/index for %s' % cd_hit_nr_fasta_file)

		try:
			cd_hit_clusters_cltr_file = cd_hit_clusters_fasta_file + '.clstr'
			assert (os.path.isfile(cd_hit_clusters_cltr_file))

			cluster = []
			with open(cd_hit_clusters_cltr_file) as off:
				for line in off:
					line = line.strip()
					ls = line.split()
					if line.startswith('>'):
						if len(cluster) > 0:
							for g in cluster:
								self.instance_to_haplotype[g] = rep
						cluster = []
						rep = None
					else:
						gene_id = ls[2][1:-3]
						cluster.append(gene_id)
						if line.endswith('*'): rep = gene_id
			if len(cluster) > 0:
				if len(cluster) > 0:
					for g in cluster:
						self.instance_to_haplotype[g] = rep

		except Exception as e:
			if self.logObject:
				self.logObject.error("Unable to parse CD-HIT clustering of gene sequences (with flanks) to obtain representative sequence per cluster.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def runSNVMining(self, paired_end_sequencing_file, bowtie2_ref_fasta, codon_alignment_file, bowtie2_alignment_dir, results_dir, debug_mode=False, cpus=1):
		"""
		Wrapper function for mining for novel SNVs across genes of GCF.

		:param paired_end_sequencing_file: tab delimited file with three columns: (1) sample name (2) path to forward
											 reads and (3) path to reverse reads
		:param bowtie2_ref_fasta: FASTA file corresponding
		:param bowtie2_alignment_dir: Path to directory where Bowtie2 alignments were written. This directory should
										include BAM files ending in *.filtered.sorted.bam which are sorted and indexed.
		:param results_dir: Path to directory where results of SNV mining will be written.
		:param cpus: The number of processes to be run in parallel.
		"""

		try:

			gene_pos_to_msa_pos = defaultdict(lambda: defaultdict(dict))
			gene_pos_to_allele = defaultdict(lambda: defaultdict(dict))
			codon_alignment_lengths = defaultdict(int)
			with open(codon_alignment_file) as ocaf:
				for line in ocaf:
					line = line.strip()
					hg, cod_alignment = line.split('\t')
					seq_count = 0
					with open(cod_alignment) as oca:
						for j, rec in enumerate(SeqIO.parse(oca, 'fasta')):
							sample_id, gene_id = rec.id.split('|')
							real_pos = 1
							for msa_pos, bp in enumerate(str(rec.seq)):
								if j == 0:
									codon_alignment_lengths[hg] += 1
								if bp != '-':
									gene_pos_to_allele[hg][gene_id][real_pos] = bp.upper()
									gene_pos_to_msa_pos[hg][gene_id][real_pos] = msa_pos+1
									real_pos += 1

							seq_count += 1

			process_args = []
			with open(paired_end_sequencing_file) as opesf:
				for line in opesf:
					line = line.strip()
					sample = line.split('\t')[0]
					process_args.append([sample, bowtie2_alignment_dir + sample + '.sorted.bam',
										 bowtie2_ref_fasta, self.instance_to_haplotype, results_dir, self.hg_genes,
										 self.comp_gene_info, dict(gene_pos_to_msa_pos), dict(gene_pos_to_allele),
										 dict(codon_alignment_lengths), debug_mode, self.logObject])

			p = multiprocessing.Pool(cpus)
			p.map(snv_miner_single, process_args)
			p.close()

		except Exception as e:
			if self.logObject:
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def calculatePairwiseDifferences(self, paired_end_sequencing_file, snv_mining_outdir, outdir):
		try:
			sample_profiles = defaultdict(lambda: defaultdict(list))
			sample_depths = defaultdict(lambda: defaultdict(float))
			with open(paired_end_sequencing_file) as opesf:
				for line in opesf:
					line = line.strip()
					pe_sample = line.strip().split('\t')[0]
					result_file = snv_mining_outdir + pe_sample + '.filt.txt'
					if not os.path.isfile(result_file): continue

					with open(result_file) as orf:
						for i, l in enumerate(orf):
							if i == 0: continue
							l = l.strip().split(',')
							hg, pos = l[:2]
							base_counts = [int(x) for x in l[2:]]
							tot_count = sum(base_counts)
							if tot_count > 0:
								base_freqs = [float(b)/float(tot_count) for b in base_counts]
							else:
								base_freqs = [0.0]*4
							sample_profiles[pe_sample][hg + '_|_' + pos] = base_freqs
							sample_depths[pe_sample][hg + '_|_' + pos] = tot_count

			pairwise_distance_file = outdir + 'sample_pairwise_differences.txt'
			sample_information_file = outdir + 'sample_information.txt'
			pairwise_distance_handle = open(pairwise_distance_file, 'w')
			sample_information_handle = open(sample_information_file, 'w')

			sample_information_handle.write('\t'.join(['sample_id', 'sample_depth']) + '\n')
			pairwise_distances_storage = defaultdict(lambda: defaultdict(lambda: 0.0))

			for si1, s1 in enumerate(sample_profiles):
				s1_hps = set(sample_profiles[s1].keys())
				s1_depth = sum(sample_depths[s1].values())/float(len(sample_depths[s1].keys()))
				if s1_depth < 10.0: continue
				sample_information_handle.write('\t'.join([s1, str(s1_depth)]) + '\n')
				for si2, s2 in enumerate(sample_profiles):
					if si1 == si2: continue
					s2_depth = sum(sample_depths[s2].values()) / float(len(sample_depths[s2].keys()))
					if s2_depth < 10.0: continue
					s2_hps = set(sample_profiles[s2].keys())
					total_intersect_positions = 0
					stat_pos = 0.0
					for hp in s1_hps.intersection(s2_hps):
						for bi, s1b in enumerate(sample_profiles[s1][hp]):
							s2b = sample_profiles[s2][hp][bi]
							stat_pos += abs(s1b - s2b)
						total_intersect_positions += 1
					distance_stat = 1.0
					if total_intersect_positions > 0:
						distance_stat = float(stat_pos)/float(total_intersect_positions)
					pairwise_distances_storage[s1][s2] = distance_stat
			sample_information_handle.close()

			pairwise_distance_handle.write('samples\t' + '\t'.join(sorted(sample_profiles)) + '\n')
			for s1 in sorted(sample_profiles):
				printlist = [s1]
				for s2 in sorted(sample_profiles):
					printlist.append(str(pairwise_distances_storage[s1][s2]))
				pairwise_distance_handle.write('\t'.join(printlist) + '\n')
			pairwise_distance_handle.close()

			# use Rscript to plot phylogeny and showcase how new sequences identified ("Query") relate to known ones
			# ("Database")
			cluster_pdf_file = outdir + 'MDS_Visualization.pdf'
			plot_cmd = ['Rscript', RSCRIPT_FOR_PCA, pairwise_distance_file, sample_information_file, cluster_pdf_file]
			if self.logObject:
				self.logObject.info('Running Rscript with the following command: %s' % ' '.join(plot_cmd))
			try:
				subprocess.call(' '.join(plot_cmd), shell=True, stdout=subprocess.DEVNULL,
								stderr=subprocess.DEVNULL,
								executable='/bin/bash')
				if self.logObject:
					self.logObject.info('Successfully ran: %s' % ' '.join(plot_cmd))
			except Exception as e:
				if self.logObject:
					self.logObject.error('Had an issue running: %s' % ' '.join(plot_cmd))
					self.logObject.error(traceback.format_exc())
				raise RuntimeError('Had an issue running: %s' % ' '.join(plot_cmd))

		except Exception as e:
			if self.logObject:
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def generateGenePhylogenies(self, codon_alignments_file, phased_alleles_outdir, comp_hg_phylo_outdir, hg_nonunique_positions, ambiguity_filter=0.00001, sequence_filter=0.25, min_number_of_sites=10):
		try:
			codon_alignment_paths = {}
			with open(codon_alignments_file) as ocaf:
				for line in ocaf:
					line = line.strip()
					ls = line.split('\t')
					codon_alignment_paths[ls[0]] = ls[1]

			for f in os.listdir(phased_alleles_outdir):
				if not f.endswith('.fasta'): continue
				hg = f.split('.fasta')[0]

				seqs = []
				ids = []
				types = []
				with open(codon_alignment_paths[hg]) as of:
					for rec in SeqIO.parse(of, 'fasta'):
						ids.append(rec.id)
						seqs.append(list(str(rec.seq).upper()))
						types.append('Database')

				cod_alg_len = len(seqs[0])
				ambiguous_positions_in_og_alignment = set([])
				for i, pos_bases in enumerate(zip(*seqs)):
					pos = i+1
					if pos <= 50 or pos >= (cod_alg_len-50):
						ambiguous_positions_in_og_alignment.add(pos)
						continue
					pos_bases = list(pos_bases)
					tot_seq_count = len(pos_bases)
					gap_seq_count = len([a for a in pos_bases if a == '-'])
					amb_prop = float(gap_seq_count)/float(tot_seq_count)
					if amb_prop >= 0.1:
						for p in range(pos - 50, pos + 51):
							ambiguous_positions_in_og_alignment.add(p)
				ambiguous_positions_in_og_alignment = ambiguous_positions_in_og_alignment.union(hg_nonunique_positions[hg])

				with open(phased_alleles_outdir + f) as of:
					for rec in SeqIO.parse(of, 'fasta'):
						seqlist = list(str(rec.seq).upper())
						gap_count = 0
						tot_count = 0
						for i, bp in enumerate(seqlist):
							pos = i+1
							if not pos in ambiguous_positions_in_og_alignment:
								tot_count += 1
								if bp == '-':
									gap_count += 1
						amb_prop = float(gap_count)/float(tot_count)
						if amb_prop < sequence_filter:
							ids.append(rec.id)
							seqs.append(seqlist)
							types.append('Query')

				ambiguous_positions_to_filter = ambiguous_positions_in_og_alignment
				for i, pos_bases in enumerate(zip(*seqs)):
					pos = i+1
					pos_bases = list(pos_bases)
					tot_seq_count = len(pos_bases)
					gap_seq_count = len([a for a in pos_bases if a == '-'])
					amb_prop = float(gap_seq_count)/float(tot_seq_count)
					if amb_prop >= ambiguity_filter:
						ambiguous_positions_to_filter.add(pos)

				gene_alignment_with_refs_filtered_file = comp_hg_phylo_outdir + hg + '.fasta'
				gene_phylogeny_with_refs_filtered_file = comp_hg_phylo_outdir + hg + '.tre'
				gene_phylogeny_track_file = comp_hg_phylo_outdir + hg + '.txt'
				gene_phylogeny_pdf_file = comp_hg_phylo_outdir + hg + '.pdf'

				gene_alignment_with_refs_filtered_handle = open(gene_alignment_with_refs_filtered_file, 'w')
				gene_phylogeny_track_handle = open(gene_phylogeny_track_file, 'w')
				gene_phylogeny_track_handle.write('name\ttype\n')

				too_few_sites_flag = False
				for i, seq in enumerate(seqs):
					id = ids[i]
					seq_filt = ''.join([a for j, a in enumerate(seq) if not (j+1) in ambiguous_positions_to_filter])
					if len(seq_filt) < min_number_of_sites:
						too_few_sites_flag = True
						break
					gene_alignment_with_refs_filtered_handle.write('>' + id + '\n' + str(seq_filt) + '\n')
					gene_phylogeny_track_handle.write(id + '\t' + types[i] + '\n')

				gene_alignment_with_refs_filtered_handle.close()
				gene_phylogeny_track_handle.close()
				if too_few_sites_flag:
					os.system('rm -f %s %s' % (gene_alignment_with_refs_filtered_file, gene_phylogeny_track_file))
					continue

				# use FastTree2 to construct gene-specific phylogeny
				fasttree_cmd = ['fasttree', '-nt', gene_alignment_with_refs_filtered_file, '>',
								gene_phylogeny_with_refs_filtered_file]
				if self.logObject:
					self.logObject.info('Running FastTree2 with the following command: %s' % ' '.join(fasttree_cmd))
				try:
					subprocess.call(' '.join(fasttree_cmd), shell=True, stdout=subprocess.DEVNULL,
									stderr=subprocess.DEVNULL,
									executable='/bin/bash')
					if self.logObject:
						self.logObject.info('Successfully ran: %s' % ' '.join(fasttree_cmd))
				except Exception as e:
					if self.logObject:
						self.logObject.error('Had an issue running: %s' % ' '.join(fasttree_cmd))
						self.logObject.error(traceback.format_exc())
					raise RuntimeError('Had an issue running: %s' % ' '.join(fasttree_cmd))

				# use Rscript to plot phylogeny and showcase how new sequences identified ("Query") relate to known ones
				# ("Database")
				plot_cmd = ['Rscript', RSCRIPT_FOR_GENERATE, gene_phylogeny_with_refs_filtered_file, gene_phylogeny_track_file, gene_phylogeny_pdf_file]
				if self.logObject:
					self.logObject.info('Running Rscript with the following command: %s' % ' '.join(plot_cmd))
				try:
					subprocess.call(' '.join(plot_cmd), shell=True, stdout=subprocess.DEVNULL,
									stderr=subprocess.DEVNULL,
									executable='/bin/bash')
					if self.logObject:
						self.logObject.info('Successfully ran: %s' % ' '.join(plot_cmd))
				except Exception as e:
					if self.logObject:
						self.logObject.error('Had an issue running: %s' % ' '.join(plot_cmd))
						self.logObject.error(traceback.format_exc())
					raise RuntimeError('Had an issue running: %s' % ' '.join(plot_cmd))

		except Exception as e:
			if self.logObject:
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def phaseAndSummarize(self, paired_end_sequencing_file, codon_alignment_file, snv_mining_outdir, phased_alleles_outdir, outdir, hg_nonunique_positions, min_hetero_prop=0.05, min_allele_depth = 5, allow_phasing=True, metagenomic=True, cpus=1):
		try:
			specific_homolog_groups = set([])
			for hg in self.hg_differentiation_stats:
				if self.hg_differentiation_stats[hg]['able_to_differentiate']:
					specific_homolog_groups.add(hg)

			gene_ignore_positions = defaultdict(set)
			gene_core_positions = defaultdict(set)
			gene_pos_to_msa_pos = defaultdict(lambda: defaultdict(dict))
			msa_pos_to_gene_allele = defaultdict(lambda: defaultdict(dict))
			gene_pos_to_allele = defaultdict(lambda: defaultdict(dict))
			msa_pos_alleles = defaultdict(lambda: defaultdict(set))
			msa_pos_ambiguous_freqs = defaultdict(lambda: defaultdict(float))
			total_core_genomes = set([])
			hg_core_genome_count = {}

			with open(codon_alignment_file) as ocaf:
				for line in ocaf:
					line = line.strip()
					hg, cod_alignment = line.split('\t')
					seq_count = 0
					msa_pos_ambiguous_counts = defaultdict(int)
					msa_pos_non_ambiguous_counts = defaultdict(int)
					seqlen_information = {}
					msa_positions = set([])
					core_genomes_with_hg = set([])
					with open(cod_alignment) as oca:
						for rec in SeqIO.parse(oca, 'fasta'):
							sequence_without_gaps = str(rec.seq).upper().replace('-', '')
							sample_id, gene_id = rec.id.split('|')
							if len(gene_id.split('_')[0]) == 3:
								core_genomes_with_hg.add(gene_id.split('_')[0])
								total_core_genomes.add(gene_id.split('_')[0])
							seqlen = len(sequence_without_gaps)
							seqlen_lower = 50
							seqlen_upper = seqlen - seqlen_lower
							seqlen_information[gene_id] = [seqlen, seqlen_lower, seqlen_upper]

							real_pos = 1
							for msa_pos, bp in enumerate(str(rec.seq)):
								msa_positions.add(msa_pos+1)
								msa_pos_to_gene_allele[hg][gene_id][msa_pos + 1] = bp.upper()
								if bp != '-':
									msa_pos_non_ambiguous_counts[msa_pos + 1] += 1
									gene_pos_to_allele[hg][gene_id][real_pos] = bp.upper()
									gene_pos_to_msa_pos[hg][gene_id][real_pos] = msa_pos + 1
									real_pos += 1
									msa_pos_alleles[hg][msa_pos + 1].add(bp.upper())
								else:
									msa_pos_ambiguous_counts[msa_pos + 1] += 1
							seq_count += 1

					for pos, seqs_with_al in msa_pos_non_ambiguous_counts.items():
						if seqs_with_al/float(seq_count) >= 0.9:
							gene_core_positions[hg].add(pos)

					cod_algn_len = max(msa_positions)
					for p in (list(range(1, 51)) + list(range(cod_algn_len-50, cod_algn_len+1))):
						gene_ignore_positions[hg].add(p)

					for pos in msa_pos_ambiguous_counts:
						msa_pos_ambiguous_freqs[hg][pos] = msa_pos_ambiguous_counts[pos] / float(seq_count)
						if (msa_pos_ambiguous_counts[pos]/float(seq_count)) >= 0.1:
							for p in range(pos - 50, pos + 51):
								gene_ignore_positions[hg].add(p)

					gene_ignore_positions[hg] = gene_ignore_positions[hg].union(hg_nonunique_positions[hg])
					hg_core_genome_count[hg] = len(core_genomes_with_hg)

			rare_hgs_in_core_genomes = set([])
			for hg in hg_core_genome_count:
				if float(hg_core_genome_count[hg])/len(total_core_genomes) < 0.5:
					rare_hgs_in_core_genomes.add(hg)

			parallel_inputs = []
			with open(paired_end_sequencing_file) as ossf:
				for line in ossf:
					pe_sample = line.strip().split('\t')[0]
					pe_sample_reads = line.strip().split('\t')[1:]
					parallel_inputs.append([pe_sample, pe_sample_reads, snv_mining_outdir, phased_alleles_outdir,
											dict(gene_ignore_positions), dict(gene_core_positions), dict(gene_pos_to_msa_pos),
											dict(msa_pos_to_gene_allele), dict(gene_pos_to_allele), dict(msa_pos_alleles),
											dict(msa_pos_ambiguous_freqs), min_hetero_prop, min_allele_depth,
											allow_phasing, metagenomic, specific_homolog_groups, set(self.core_homologs),
											dict(self.hg_genes), dict(self.comp_gene_info),
											dict(self.hg_prop_multi_copy), set(self.protocluster_core_homologs),
											rare_hgs_in_core_genomes, self.gcf_id, self.logObject])

			p = multiprocessing.Pool(cpus)
			p.map(phase_and_id_snvs, parallel_inputs)
			p.close()

			novelty_report_file = outdir + 'Potentially_Novel_SNVs.txt'
			homolog_presence_report_file = outdir + 'Sample_Homolog_Group_Coverage.txt'
			no_handle = open(novelty_report_file, 'w')
			hpr_handle = open(homolog_presence_report_file, 'w')
			no_handle.write('\t'.join(['gcf_id', 'sample', 'homolog_group', 'position_along_msa', 'site_total_coverage_standardized',
									   'site_allele_coverage_standardized', 'alternate_allele_is_major_allele', 'alternate_allele',
									   'codon_position', 'alternate_codon', 'alternate_aa', 'alternate_aa_is_novel', 'dn_or_ds', 'ts_or_tv',
									   'reference_allele', 'reference_sample', 'reference_gene', 'reference_position',
									   'ref_codon', 'ref_aa', 'snv_support', 'snv_support_read_ids']) + '\n')

			hpr_handle.write('\t'.join(['sample', 'homolog_group', 'outlier_in_coverage',
										'homolog_group_is_differentiable_from_paralogs',
										'homolog_group_proportion_initial_samples_with_multi_copy',
										'product_has_MGE_related_term', 'homolog_group_trimmed_median_depth',
										'homolog_group_early_truncation_stop_codon_position', 'difficult_to_resolve_positions']) + '\n')

			with open(paired_end_sequencing_file) as ossf:
				for line in ossf:
					pe_sample = line.strip().split('\t')[0]
					pes_novelty_report_file = snv_mining_outdir + pe_sample + '.novel_snvs_report.txt'
					pes_group_coverage_file = snv_mining_outdir + pe_sample + '.homolog_group_coverage.txt'

					with open(pes_novelty_report_file) as opnrf:
						for i, line in enumerate(opnrf):
							if i == 0: continue
							no_handle.write(line)

					with open(pes_group_coverage_file) as opgcf:
						for i, line in enumerate(opgcf):
							if i == 0: continue
							hpr_handle.write(line)

			no_handle.close()
			hpr_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues with phasing and identifying potentially novel SNVs.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())


def phase_and_id_snvs(input_args):
	pe_sample, pe_sample_reads, snv_mining_outdir, phased_alleles_outdir, gene_ignore_positions, gene_core_positions, gene_pos_to_msa_pos, msa_pos_to_gene_allele, gene_pos_to_allele, msa_pos_alleles, msa_pos_ambiguous_freqs, min_hetero_prop, min_allele_depth, allow_phasing, metagenomic, specific_homolog_groups, core_homologs, hg_genes, comp_gene_info, hg_prop_multi_copy, protocluster_core_homologs, rare_hgs_in_core_genomes, gcf_id, logObject = input_args
	try:
		result_file = snv_mining_outdir + pe_sample + '.txt'
		snv_file = snv_mining_outdir + pe_sample + '.snvs'
		homolog_group_positions = defaultdict(list)
		homolog_group_depths = defaultdict(list)
		previous_positions = []

		novelty_report_file = snv_mining_outdir + pe_sample + '.novel_snvs_report.txt'
		homolog_presence_report_file = snv_mining_outdir + pe_sample + '.homolog_group_coverage.txt'
		no_handle = open(novelty_report_file, 'w')
		hpr_handle = open(homolog_presence_report_file, 'w')
		no_handle.write('\t'.join(['gcf_id', 'sample', 'homolog_group', 'position_along_msa', 'site_total_coverage_standardized',
									'site_allele_coverage_standardized', 'alternate_allele_is_major_allele', 'alternate_allele',
								   'codon_position', 'alternate_codon', 'alternate_aa', 'alternate_aa_is_novel', 'dn_or_ds', 'ts_or_tv',
								   'reference_allele', 'reference_sample', 'reference_gene', 'reference_position',
								   'ref_codon', 'ref_aa', 'snv_support', 'snv_support_read_ids']) + '\n')
		hpr_handle.write('\t'.join(['sample', 'homolog_group', 'outlier_in_coverage',
									'homolog_group_is_differentiable_from_paralogs',
									'homolog_group_proportion_initial_samples_with_multi_copy',
									'product_has_MGE_related_term', 'homolog_group_trimmed_median_depth',
									'homolog_group_early_truncation_stop_codon_position',
									'difficult_to_resolve_positions']) + '\n')

		hg_median_depths = defaultdict(float)
		hg_first_position_of_stop_codon = defaultdict(lambda: None)

		if not os.path.isfile(result_file): return
		hg_hetero_sites = defaultdict(set)
		with open(result_file) as orf:
			for i, line in enumerate(orf):
				if i == 0: continue
				line = line.strip()
				homolog_group = line.split(',')[0]
				position, a, c, g, t = [int(x) for x in line.split(',')[1:]]
				total_depth = sum([a,c,g,t])

				adequate_coverage_alleles = []
				if a >= min_allele_depth:
					adequate_coverage_alleles.append('A')
				if c >= min_allele_depth:
					adequate_coverage_alleles.append('C')
				if g >= min_allele_depth:
					adequate_coverage_alleles.append('G')
				if t >= min_allele_depth:
					adequate_coverage_alleles.append('T')

				if len(adequate_coverage_alleles) > 1:
					hg_hetero_sites[homolog_group].add(position)

				if position % 3 == 0:
					if len(previous_positions) == 2:
						for al1 in previous_positions[0]:
							for al2 in previous_positions[1]:
								for al3 in adequate_coverage_alleles:
									codon = al1 + al2 + al3
									if codon in set(['TAG', 'TGA', 'TAA']):
										if not homolog_group in hg_first_position_of_stop_codon:
											hg_first_position_of_stop_codon[homolog_group] = position-2
					previous_positions = []
				else:
					previous_positions.append(adequate_coverage_alleles)

				homolog_group_depths[homolog_group].append(total_depth)
				homolog_group_positions[homolog_group].append(position)

		present_homolog_groups = set([])
		for hg in homolog_group_depths:
			hg_depths = homolog_group_depths[hg]
			hg_positions = homolog_group_positions[hg]

			if hg in rare_hgs_in_core_genomes: continue
			if not hg in gene_core_positions.keys(): continue
			core_positions_covered = 0
			for i, pos in enumerate(hg_positions):
				dp = hg_depths[i]
				if pos in gene_core_positions[hg] and dp >= 1:
					core_positions_covered += 1

			if float(core_positions_covered)/len(gene_core_positions[hg]) >= 0.9:
				present_homolog_groups.add(hg)
				filtered_hg_depths = []
				for pos in range(1, len(hg_depths)+1):
					if pos in gene_ignore_positions[hg]: continue
					filtered_hg_depths.append(hg_depths[pos-1])
				if len(filtered_hg_depths) > 10:
					hg_filtered_depth_median = statistics.median(filtered_hg_depths)
					hg_median_depths[hg] = hg_filtered_depth_median

		if len(hg_median_depths) < 5:
			no_handle.close()
			hpr_handle.close()
			return
		median_of_medians = statistics.median(list(hg_median_depths.values()))
		mad_of_medians = median_abs_deviation(list(hg_median_depths.values()), scale="normal")

		outlier_homolog_groups = set([])
		for hg in present_homolog_groups:
			hg_median = hg_median_depths[hg]
			if hg_median > (median_of_medians + (2*mad_of_medians)) or hg_median < (median_of_medians - (2*mad_of_medians)):
				outlier_homolog_groups.add(hg)

		depths_at_all_refined_present_hgs = []
		hetero_sites = 0
		total_sites = 0
		mge_hgs = set([])
		refined_present_homolog_groups = set([])
		report_lines = []
		for hg in hg_genes:
			product_has_mge_term = False
			for gene in hg_genes[hg]:
				for word in mges:
					if word in comp_gene_info[gene]['product'].lower():
						product_has_mge_term = True
					for domain_dict in comp_gene_info[gene]['gene_domains']:
						if word in domain_dict['description'].lower():
							product_has_mge_term = True

			report_lines.append('\t'.join([str(x) for x in [pe_sample, hg, (hg in outlier_homolog_groups),
														 (hg in specific_homolog_groups),
														 hg_prop_multi_copy[hg], product_has_mge_term,
														 hg_median_depths[hg],
														 hg_first_position_of_stop_codon[hg],
														 ','.join([str(x) for x in sorted(gene_ignore_positions[hg])])]]))

			if product_has_mge_term: mge_hgs.add(hg)
			if not product_has_mge_term and not hg in outlier_homolog_groups and hg_median_depths[hg] > 0.0:
				refined_present_homolog_groups.add(hg)
				if hg_prop_multi_copy[hg] < 0.05:
					hg_depths = homolog_group_depths[hg]

					for pos in range(1, len(hg_depths)+1):
						if pos in gene_ignore_positions[hg]: continue
						depths_at_all_refined_present_hgs.append(hg_depths[pos-1])
						total_sites += 1
						if pos in hg_hetero_sites[hg]:
							hetero_sites += 1

		#hpr_handle.write('\n'.join(report_lines) + '\n')

		if (len(refined_present_homolog_groups) < 5) or (len(refined_present_homolog_groups.intersection(protocluster_core_homologs)) == 0) or (len(refined_present_homolog_groups.intersection(core_homologs))/float(len(core_homologs.difference(mge_hgs))) < 0.7 and len(refined_present_homolog_groups.intersection(specific_homolog_groups)) == 0):
			no_handle.close()
			hpr_handle.close()
			return

		hpr_handle.write('\n'.join(report_lines) + '\n')

		filt_result_file = snv_mining_outdir + pe_sample + '.filt.txt'
		filt_result_handle = open(filt_result_file, 'w')
		pos_allele_support = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
		with open(result_file) as orf:
			for linenum, line in enumerate(orf):
				line = line.strip()
				if linenum == 0:
					filt_result_handle.write(line + '\n')
				else:
					hg = line.split(',')[0]
					if not hg in refined_present_homolog_groups or hg_prop_multi_copy[hg] >= 0.05: continue
					pos, a, c, g, t = [int(x) for x in line.split(',')[1:]]
					pos_allele_support[hg][pos]['A'] = a
					pos_allele_support[hg][pos]['C'] = c
					pos_allele_support[hg][pos]['G'] = g
					pos_allele_support[hg][pos]['T'] = t
					pos_allele_support[hg][pos]['TOTAL'] = sum([a,c,g,t])
					filt_result_handle.write(line + '\n')
		filt_result_handle.close()

		trimmed_depth_median = statistics.median(depths_at_all_refined_present_hgs)
		trimmed_depth_mad = median_abs_deviation(depths_at_all_refined_present_hgs, scale="normal")

		# Check if phasing needed
		homolog_variable_positions = defaultdict(set)
		haplotype_allele_at_position = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: None)))
		number_of_haplotypes = 0

		if float(hetero_sites)/total_sites >= min_hetero_prop and allow_phasing and metagenomic:
			# perform phasing using desman
			cwd = os.getcwd()

			desman_general_dir = snv_mining_outdir + pe_sample + '_Desman_Dir/'
			desman_variants_dir = desman_general_dir + 'Variants/'
			desman_inferstrains_dir = desman_general_dir + 'InferStrains/'

			if not os.path.isdir(desman_general_dir): os.system('mkdir %s' % desman_general_dir)
			if not os.path.isdir(desman_variants_dir): os.system('mkdir %s' % desman_variants_dir)
			if not os.path.isdir(desman_inferstrains_dir): os.system('mkdir %s' % desman_inferstrains_dir)

			desman_variant_filter_cmd = ['cd', desman_variants_dir, ';', 'Variant_Filter.py', filt_result_file,
										 ';', 'cd', cwd]

			if logObject:
				logObject.info(
					'Running Desman variant filtering with the following command: %s' % ' '.join(desman_variant_filter_cmd))
			desman_issue = False
			try:
				subprocess.call(' '.join(desman_variant_filter_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
				logObject.info('Successfully ran: %s' % ' '.join(desman_variant_filter_cmd))
			except Exception as e:
				if logObject:
					logObject.error('Had an issue running: %s' % ' '.join(desman_variant_filter_cmd))
					logObject.error(traceback.format_exc())
				desman_issue = True
				#raise RuntimeWarning('Had an issue running: %s' % ' '.join(desman_variant_filter_cmd))

			if not desman_issue:
				try:
					for g in [2, 3, 4, 5, 6, 7, 8]:
						for r in [0, 1, 2, 3, 4]:
							desmand_inferstrains_cmd = ['cd', desman_inferstrains_dir, ';', 'desman', '-e',
														desman_variants_dir + 'outputtran_df.csv', '-o',
														'ClusterEC_' + str(g) + '_' + str(r), '-r', '1000', '-i',
														'100', '-g', str(g), '-s', str(r),
														desman_variants_dir + 'outputsel_var.csv', '>',
														desman_inferstrains_dir + 'ClusterEC_' + str(g) + '_' + str(r) + '.out']

							if logObject:
								logObject.info(
									'Running Desman for strain inference with the following command: %s' % ' '.join(desmand_inferstrains_cmd))
							try:
								subprocess.call(' '.join(desmand_inferstrains_cmd), shell=True, stdout=sys.stderr, stderr=sys.stderr, executable='/bin/bash')
								logObject.info('Successfully ran: %s' % ' '.join(desmand_inferstrains_cmd))
							except Exception as e:
								if logObject:
									logObject.error(
										'Had an issue running: %s' % ' '.join(desmand_inferstrains_cmd))
									logObject.error(traceback.format_exc())
								raise RuntimeError('Had an issue running: %s' % ' '.join(desmand_inferstrains_cmd))

					desman_resolvehap_cmd = ['cd', desman_inferstrains_dir, ';', 'resolvenhap.py', 'ClusterEC', '>',
											 desman_general_dir + 'Best_Parameter_Combo.txt']
					if logObject:
						logObject.info(
							'Assessing Desman runs for strain inference with the following command: %s' % ' '.join(
								desman_resolvehap_cmd))
					try:
						subprocess.call(' '.join(desman_resolvehap_cmd), shell=True, stdout=sys.stderr, stderr=sys.stderr, executable='/bin/bash')
						logObject.info('Successfully ran: %s' % ' '.join(desman_resolvehap_cmd))
					except Exception as e:
						if logObject:
							logObject.error('Had an issue running: %s' % ' '.join(desman_resolvehap_cmd))
							logObject.error(traceback.format_exc())
						raise RuntimeError('Had an issue running: %s' % ' '.join(desman_resolvehap_cmd))

					haps, conf_haps, seed, avg_error = [None]*4
					with open(desman_general_dir + 'Best_Parameter_Combo.txt') as obpc:
						for line in obpc:
							line = line.strip()
							haps, conf_haps, seed, avg_error, _ = line.split(',')

					haplotype_spec_file = desman_inferstrains_dir + 'ClusterEC_' + str(haps) + '_' + str(seed) + '/Filtered_Tau_star.csv'
					with open(haplotype_spec_file) as obhpf:
						for linenum, line in enumerate(obhpf):
							if linenum == 0: continue
							line = line.strip()
							ls = line.split(',')
							hg, position = ls[:2]
							position = int(position)
							homolog_variable_positions[hg].add(position)

							ambiguity_region_flag = False
							if position in gene_ignore_positions[hg]: ambiguity_region_flag = True
							haplotype_calls = ls[2:]

							total_depth = pos_allele_support[hg][position]['TOTAL']
							depth_above_expectation = False
							depth_below_expectation = False
							if total_depth > (trimmed_depth_median + (2 * trimmed_depth_mad)):
								depth_above_expectation = True
							if total_depth < (trimmed_depth_median - (3 * trimmed_depth_mad)):
								depth_below_expectation = True

							for i, allele_call in enumerate(haplotype_calls):
								haplotype_num = math.floor(i/4)
								if haplotype_num > number_of_haplotypes: number_of_haplotypes = haplotype_num
								if allele_call == '1':
									if (i+1) % 4 == 0:
										allele_base = 'T'
									elif (i+1) % 3 == 0:
										allele_base = 'G'
									elif (i+1) % 2 == 0:
										allele_base = 'C'
									else:
										allele_base = 'A'

									allele_depth = pos_allele_support[hg][position][allele_base]

									if allele_depth >= min_allele_depth and not depth_below_expectation and \
											not depth_above_expectation and not ambiguity_region_flag and \
											haplotype_allele_at_position[hg][position] == None:
										haplotype_allele_at_position[hg][haplotype_num][position] = allele_base
									else:
										haplotype_allele_at_position[hg][haplotype_num][position] = '-'
				except:
					if logObject:
						logObject.error('Had an issue with parsing DESMAN results.')
						logObject.error(traceback.format_exc())


		haplotype_sequences = defaultdict(lambda: defaultdict(lambda: ""))
		with open(filt_result_file) as ofrf:
			for linenum, line in enumerate(ofrf):
				if linenum == 0: continue
				line = line.strip()
				ls = line.split(',')
				hg, position = ls[:2]
				position = int(position)

				ambiguity_region_flag = False
				if position in gene_ignore_positions[hg]: ambiguity_region_flag = True
				haplotype_calls = [int(x) for x in ls[2:]]

				total_depth = pos_allele_support[hg][position]['TOTAL']
				depth_above_expectation = False
				depth_below_expectation = False
				if total_depth > (trimmed_depth_median + (2 * trimmed_depth_mad)):
					depth_above_expectation = True
				if total_depth < (trimmed_depth_median - (3 * trimmed_depth_mad)):
					depth_below_expectation = True

				max_allele_depth = max(haplotype_calls)

				tie_exists = len([x for x in haplotype_calls if x == max_allele_depth]) > 1

				for i, allele_depth in enumerate(haplotype_calls):
					if allele_depth == max_allele_depth:
						if i == 0: allele_base = 'A'
						elif i == 1: allele_base = 'C'
						elif i == 2: allele_base = 'G'
						elif i == 3: allele_base = 'T'

						allele_call = '-'
						if allele_depth >= min_allele_depth and not depth_below_expectation and not depth_above_expectation and \
								not ambiguity_region_flag and not tie_exists:
							allele_call = allele_base

						for hi in range(0, number_of_haplotypes+1):
							if position in homolog_variable_positions[hg]:
								haplotype_sequences[hg][hi] += haplotype_allele_at_position[hg][hi][position]
							else:
								haplotype_sequences[hg][hi] += allele_call
						break

		for hg in haplotype_sequences:
			bgc_fasta_file = phased_alleles_outdir + hg + '.fasta'
			bgc_fasta_handle = open(bgc_fasta_file, 'a+')
			for hi in haplotype_sequences[hg]:
				seq = haplotype_sequences[hg][hi]
				codons = [str(seq)[i:i + 3] for i in range(0, len(str(seq)), 3)]
				first_stop_codon = None
				for cod_i, cod in enumerate(codons):
					if cod in set(['TAG', 'TGA', 'TAA']):
						first_stop_codon = 3*(cod_i+1)
						break
				if first_stop_codon is not None:
					seq = seq[:first_stop_codon] + ''.join(['-']*len(seq[first_stop_codon:]))
				bgc_fasta_handle.write('>' + pe_sample + '_|_' + str(hi+1) + '\n' + seq + '\n')
			bgc_fasta_handle.close()


		all_snv_supporting_reads = set([])
		with open(snv_file) as of:
			for i, line in enumerate(of):
				line = line.strip()
				ls = line.split('\t')
				snv_id, snv_support_count, snv_support_reads = ls

				gsh, ref_pos, ref_al, alt_al = snv_id.split('_|_')
				hg, allele_cluster, sample, gene = gsh.split('|')

				ref_pos = int(ref_pos)

				# check whether to regard this SNV as potentially legit
				if not hg in refined_present_homolog_groups or hg_prop_multi_copy[hg] >= 0.05: continue
				if not ref_pos in gene_pos_to_msa_pos[hg][gene]: continue
				msa_pos = gene_pos_to_msa_pos[hg][gene][ref_pos]
				msa_pos_als = msa_pos_alleles[hg][msa_pos]
				if hg_first_position_of_stop_codon[hg] != None and msa_pos >= hg_first_position_of_stop_codon[hg]: continue
				if int(snv_support_count) >= min_allele_depth and int(homolog_group_depths[hg][msa_pos-1]) <= (trimmed_depth_median+(2*trimmed_depth_mad)) and int(homolog_group_depths[hg][msa_pos-1]) >= trimmed_depth_median-(3*trimmed_depth_mad):
					assert (ref_al in msa_pos_alleles[hg][msa_pos] and ref_al == gene_pos_to_allele[hg][gene][ref_pos])
					if not alt_al in msa_pos_als and not msa_pos in gene_ignore_positions[hg]:
						#if not alt_al in msa_pos_als and not msa_pos in gene_edgy_positions[hg] and msa_pos_ambiguous_freqs[hg][msa_pos] <= 0.1:
						codon_position = None
						ref_codon = None
						ref_aa = None
						alt_codon = None
						alt_aa = None
						alt_aa_novel = True
						dn_or_ds = None
						ts_or_tv = "transition"
						if (ref_al in purine_alleles) != (alt_al in purine_alleles):
							ts_or_tv = "transversion"
						if ref_pos % 3 == 1:
							codon_position = 1
							ref_codon = ref_al + gene_pos_to_allele[hg][gene][ref_pos + 1] + \
										gene_pos_to_allele[hg][gene][ref_pos + 2]
							alt_codon = alt_al + gene_pos_to_allele[hg][gene][ref_pos + 1] + \
										gene_pos_to_allele[hg][gene][ref_pos + 2]
							ref_aa = str(Seq(ref_codon).translate())
							alt_aa = str(Seq(alt_codon).translate())

							for g in gene_pos_to_allele[hg]:
								gene_codon = msa_pos_to_gene_allele[hg][g][msa_pos] + \
										   msa_pos_to_gene_allele[hg][g][msa_pos + 1] + \
										   msa_pos_to_gene_allele[hg][g][msa_pos + 2]
								if g == gene:
									assert(gene_codon == ref_codon)
								gene_aa = str(Seq(gene_codon).translate())
								if alt_aa == gene_aa: alt_aa_novel = False
						elif ref_pos % 3 == 2:
							codon_position = 2
							ref_codon = gene_pos_to_allele[hg][gene][ref_pos - 1] + ref_al + \
										gene_pos_to_allele[hg][gene][ref_pos + 1]
							alt_codon = gene_pos_to_allele[hg][gene][ref_pos - 1] + alt_al + \
										gene_pos_to_allele[hg][gene][ref_pos + 1]
							ref_aa = str(Seq(ref_codon).translate())
							alt_aa = str(Seq(alt_codon).translate())

							for g in gene_pos_to_allele[hg]:
								gene_codon = msa_pos_to_gene_allele[hg][g][msa_pos - 1] + \
										   msa_pos_to_gene_allele[hg][g][msa_pos] + \
										   msa_pos_to_gene_allele[hg][g][msa_pos + 1]
								if g == gene:
									assert(gene_codon == ref_codon)
								gene_aa = str(Seq(gene_codon).translate())
								if alt_aa == gene_aa: alt_aa_novel = False
						elif ref_pos % 3 == 0:
							codon_position = 3
							ref_codon = gene_pos_to_allele[hg][gene][ref_pos - 2] + \
										gene_pos_to_allele[hg][gene][ref_pos - 1] + ref_al
							alt_codon = gene_pos_to_allele[hg][gene][ref_pos - 2] + \
										gene_pos_to_allele[hg][gene][ref_pos - 1] + alt_al
							ref_aa = str(Seq(ref_codon).translate())
							alt_aa = str(Seq(alt_codon).translate())

							for g in gene_pos_to_allele[hg]:
								gene_codon = msa_pos_to_gene_allele[hg][g][msa_pos - 2] + \
										   msa_pos_to_gene_allele[hg][g][msa_pos - 1] + \
										   msa_pos_to_gene_allele[hg][g][msa_pos]
								if g == gene:
									assert(gene_codon == ref_codon)
								gene_aa = str(Seq(gene_codon).translate())
								if alt_aa == gene_aa: alt_aa_novel = False

						if ref_aa != alt_aa:
							dn_or_ds = "non-synonymous"
						else:
							dn_or_ds = "synonymous"

						site_total_coverage = pos_allele_support[hg][msa_pos]['TOTAL']
						site_allele_coverage = pos_allele_support[hg][msa_pos][alt_al]

						site_total_coverage_standardized = (site_total_coverage - trimmed_depth_median)/float(trimmed_depth_mad)
						site_allele_coverage_standardized = (site_allele_coverage - trimmed_depth_median) / float(trimmed_depth_mad)

						site_major_allele_count = max([pos_allele_support[hg][msa_pos]['A'],
													   pos_allele_support[hg][msa_pos]['C'],
													   pos_allele_support[hg][msa_pos]['G'],
													   pos_allele_support[hg][msa_pos]['T']])
						snv_is_major_allele = False
						if site_major_allele_count == site_allele_coverage:
							snv_is_major_allele = True

						no_handle.write('\t'.join([str(x) for x in [gcf_id, pe_sample, hg, msa_pos,
															site_total_coverage_standardized,
															site_allele_coverage_standardized, snv_is_major_allele,
															alt_al, codon_position, alt_codon, alt_aa, alt_aa_novel, dn_or_ds,
															ts_or_tv, ref_al, sample, gene, ref_pos,
															ref_codon, ref_aa, snv_support_count, snv_support_reads]]) + '\n')
						all_snv_supporting_reads = all_snv_supporting_reads.union(set(snv_support_reads.split(', ')))

		snv_support_fastq_file = snv_mining_outdir + pe_sample + '.snv_support.fastq'
		snv_support_fastq_handle = open(snv_support_fastq_file, 'w')
		visited = set([])
		for read_file in pe_sample_reads:
			if read_file.endswith('.gz'):
				fastq_handle = gzip.open(read_file, 'rt')
			else:
				fastq_handle = open(read_file)

			lines = []
			for line in fastq_handle:
				lines.append(line.rstrip())
				if len(lines) == 4:
					if lines[0][1:] in all_snv_supporting_reads:
						visited.add(lines[0][1:])
						snv_support_fastq_handle.write('\n'.join(lines) + '\n')
					lines = []

			fastq_handle.close()

		"""
		for r in all_snv_supporting_reads:
			if not r in visited:
				print(pe_sample + '\t' + r)
		"""
		snv_support_fastq_handle.close()
		os.system('gzip %s' % snv_support_fastq_file)
		no_handle.close()
		hpr_handle.close()
	except Exception as e:
		if logObject:
			logObject.error('Difficulties with phasing and identifying potentially novel SNVs.')
			logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

def snv_miner_single(input_args):
	"""
	Function to parse BAM alignment files (assumes reads aligned independently - not paired - CURRENT DEFAULT).
	"""
	sample, bam_alignment, ref_fasta, hg_gene_to_rep, res_dir, bgc_hg_genes, comp_gene_info, gene_pos_to_msa_pos, gene_pos_to_allele, codon_alignment_lengths, debug_mode, logObject = input_args
	try:
		hg_rep_genes = defaultdict(set)
		for g, r in hg_gene_to_rep.items():
			hg_rep_genes[r].add(g)

		if not os.path.isfile(bam_alignment): return
		snvs_file = res_dir + sample + '.snvs'
		result_file = res_dir + sample + '.txt'
		details_file = res_dir + sample + '.full.txt'
		snv_outf = open(snvs_file, 'w')
		res_outf = open(result_file, 'w')
		det_outf = None
		if debug_mode:
			det_outf = open(details_file, 'w')

		res_outf.write('Contig,Position,Sample-A,Sample-C,Sample-G,Sample-T\n')

		bam_handle = pysam.AlignmentFile(bam_alignment, 'rb')

		topaligns_file = res_dir + sample + '_topaligns.bam'
		topaligns_file_sorted = res_dir + sample + '_topaligns.sorted.bam'
		topaligns_handle = pysam.AlignmentFile(topaligns_file, "wb", template=bam_handle)

		for hg, hg_genes in bgc_hg_genes.items():
			read_ascpus_per_allele = defaultdict(list)
			hg_genes_covered = 0
			gene_sequence = {}
			total_reads = set([])
			with open(ref_fasta) as opff:
				for rec in SeqIO.parse(opff, 'fasta'):
					if rec.id.split('|')[0] != hg: continue
					_, allele_cluster, _, g = rec.id.split('|')
					ginfo = comp_gene_info[g]
					gene_sequence[g] = str(rec.seq)
					gstart = ginfo['start']
					gend = ginfo['end']

					gene_length = gend - gstart + 1

					gene_covered_1 = 0

					try:
						for pileupcolumn in bam_handle.pileup(contig=rec.id, stepper="nofilter"):
							pos_depth = 0
							for pileupread in pileupcolumn.pileups:
								if pileupread.is_del or pileupread.is_refskip: continue
								read = pileupread.alignment
								if read.query_qualities[pileupread.query_position] < 30: continue
								pos_depth += 1
							if pos_depth >= 1:
								gene_covered_1 += 1
					except:
						pass

					gene_coverage_1 = gene_covered_1 / float(gene_length)
					if gene_coverage_1 < 0.90: continue
					hg_genes_covered += 1
					#print('\t'.join([sample, hg, rec.id, str(gene_coverage_1), str(gene_coverage_3)]))

					for read_alignment in bam_handle.fetch(rec.id):
						read_name = read_alignment.query_name
						total_reads.add(read_name)
						read_ascore = read_alignment.tags[0][1]

						read_ref_positions = set(read_alignment.get_reference_positions())

						first_real_alignment_pos = None
						last_real_alignment_pos = None
						indel_positions = set([])
						matches = set([])
						for b in read_alignment.get_aligned_pairs(with_seq=True):
							if not (b[0] == None or b[1] == None):
								if first_real_alignment_pos == None:
									first_real_alignment_pos = b[1]
								last_real_alignment_pos = b[1]
								if b[2].isupper():
									matches.add(b[1])
							else:
								indel_positions.add(b[1])

						main_alignment_positions = set(range(first_real_alignment_pos, last_real_alignment_pos + 1))
						sum_indel_len = len(main_alignment_positions.intersection(indel_positions))
						matching_percentage = float(len(matches))/float(len(main_alignment_positions))

						read_ascpus_per_allele[read_name].append([g, read_ascore, matching_percentage, len(main_alignment_positions), sum_indel_len, read_alignment])

			accounted_reads = set([])
			hg_align_pos_alleles = defaultdict(lambda: defaultdict(set))
			supported_snvs = defaultdict(lambda: defaultdict(set))
			for read in read_ascpus_per_allele:
				top_score = -1000000
				score_sorted_alignments = sorted(read_ascpus_per_allele[read], key=itemgetter(1), reverse=True)
				for i, align in enumerate(score_sorted_alignments):
					if i == 0: top_score = align[1]
					if align[1] == top_score and ((align[2] >= 0.99 and align[3] >= 60) or (align[2] >= 0.95 and align[3] >= 100)) and align[4] <= 5:
						read_alignment = align[-1]
						topaligns_handle.write(read_alignment)

						min_read_ref_pos = min(read_alignment.get_reference_positions())
						read_referseq = read_alignment.get_reference_sequence().upper()
						read_queryseq = read_alignment.query_sequence
						read_queryqua = read_alignment.query_qualities

						for b in read_alignment.get_aligned_pairs(with_seq=True):
							if b[0] == None or b[1] == None: continue
							ref_pos = b[1]+1
							alt_al = read_queryseq[b[0]].upper()
							ref_al = read_referseq[b[1] - min_read_ref_pos].upper()
							assert (ref_al == str(gene_sequence[align[0]]).upper()[b[1]])
							if b[2] == 'n' or ref_al == 'N' or alt_al == 'N': continue
							que_qual = read_queryqua[b[0]]
							if (que_qual >= 30) and ((ref_pos+3) < len(gene_sequence[align[0]])):
								cod_pos = gene_pos_to_msa_pos[hg][align[0]][ref_pos]
								hg_align_pos_alleles[cod_pos][alt_al].add(read)
								accounted_reads.add(read)
								if debug_mode:
									det_outf.write('\t'.join([str(x) for x in [sample, hg, align[0], ref_pos, cod_pos, ref_al, alt_al, read, align[1], align[2], align[3], align[4]]]) + '\n')
								if b[2].islower():
									assert (alt_al != ref_al)
									# because reads with at most 5 indel positions allowed above, a base might appear in the consensus/phased haplotypes
									# that differs from known reference alleles, but not be called as an SNV
									if align[4] == 0:
										snv_id = str(read_alignment.reference_name) + '_|_' + str(ref_pos) + '_|_' + ref_al + '_|_' + alt_al
										supported_snvs[snv_id][read].add(align[1])

			for pos in range(1, codon_alignment_lengths[hg]+1):
				printlist = [hg, str(pos)]
				for al in ['A', 'C', 'G', 'T']:
					printlist.append(str(len(hg_align_pos_alleles[pos][al])))
				res_outf.write(','.join(printlist) + '\n')

			for snv in supported_snvs:
				support_info = []
				for read in supported_snvs[snv]:
					support_info.append(read) # + '_|_' + str(max(supported_snvs[snv][read])))
				snv_outf.write('\t'.join([snv, str(len(support_info)), ', '.join(support_info)]) + '\n')

			#print('\t'.join([hg, str(len(total_reads)), str(len(accounted_reads))]))

		snv_outf.close()
		res_outf.close()
		if debug_mode:
			det_outf.close()
		topaligns_handle.close()
		bam_handle.close()

		os.system("samtools sort -@ %d %s -o %s" % (1, topaligns_file, topaligns_file_sorted))
		os.system("samtools index %s" % topaligns_file_sorted)

	except Exception as e:
		if logObject:
			logObject.error("Issues with mining for SNVs/parsing alignment file.")
			logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())


def snv_miner_paired(input_args):
	"""
	Function to mine for novel SNVs and identify alleles of homolog groups for GCF present in paired-end sequencing
	dataset.
	"""
	sample, bam_alignment, ref_fasta, hg_gene_to_rep, res_dir, bgc_hg_genes, comp_gene_info, gene_pos_to_msa_pos, gene_pos_to_allele, codon_alignment_lengths, logObject = input_args
	try:
		hg_rep_genes = defaultdict(set)
		for g, r in hg_gene_to_rep.items():
			hg_rep_genes[r].add(g)

		if not os.path.isfile(bam_alignment): return
		snvs_file = res_dir + sample + '.snvs'
		result_file = res_dir + sample + '.txt'
		details_file = res_dir + sample + '.full.txt'
		snv_outf = open(snvs_file, 'w')
		res_outf = open(result_file, 'w')
		det_outf = open(details_file, 'w')

		res_outf.write('Contig,Position,Sample-A,Sample-C,Sample-G,Sample-T\n')

		bam_handle = pysam.AlignmentFile(bam_alignment, 'rb')

		topaligns_file = res_dir + sample + '_topaligns.bam'
		topaligns_file_sorted = res_dir + sample + '_topaligns.sorted.bam'
		topaligns_handle = pysam.AlignmentFile(topaligns_file, "wb", template=bam_handle)

		for hg, hg_genes in bgc_hg_genes.items():
			read_ascpus_per_allele = defaultdict(list)
			hg_genes_covered = 0
			gene_sequence = {}
			total_reads = set([])
			with open(ref_fasta) as opff:
				for rec in SeqIO.parse(opff, 'fasta'):
					if rec.id.split('|')[0] != hg: continue
					_, allele_cluster, _, g = rec.id.split('|')
					ginfo = comp_gene_info[g]
					gene_sequence[g] = str(rec.seq)
					gstart = ginfo['start']
					gend = ginfo['end']

					gene_length = gend - gstart + 1

					gene_covered_1 = 0

					try:
						for pileupcolumn in bam_handle.pileup(contig=rec.id, stepper="nofilter"):
							pos_depth = 0
							for pileupread in pileupcolumn.pileups:
								if pileupread.is_del or pileupread.is_refskip: continue
								read = pileupread.alignment
								if read.query_qualities[pileupread.query_position] < 30: continue
								pos_depth += 1
							if pos_depth >= 1:
								gene_covered_1 += 1
					except:
						pass

					gene_coverage_1 = gene_covered_1 / float(gene_length)
					if gene_coverage_1 < 0.90: continue
					hg_genes_covered += 1
					#print('\t'.join([sample, hg, rec.id, str(gene_coverage_1), str(gene_coverage_3)]))

					for read1_alignment, read2_alignment in util.read_pair_generator(bam_handle, rec.id, gene_length):
						if read1_alignment and read2_alignment:
							read_name = read1_alignment.query_name
							total_reads.add(read_name)
							combined_ascore = read1_alignment.tags[0][1] + read2_alignment.tags[0][1]
							read1_ref_positions = set(read1_alignment.get_reference_positions())
							read2_ref_positions = set(read2_alignment.get_reference_positions())

							align_union = len(read1_ref_positions.union(read2_ref_positions))
							align_intersect = len(read1_ref_positions.intersection(read2_ref_positions))
							min_align_length = min(len(read1_ref_positions), len(read2_ref_positions))
							align_overlap_prop = float(align_intersect) / float(align_union) # / min_align_length

							matches_1 = set([])
							first_real_alignment_pos = None
							last_real_alignment_pos = None
							indel_positions = set([])
							for b in read1_alignment.get_aligned_pairs(with_seq=True):
								if not (b[0] == None or b[1] == None):
									if first_real_alignment_pos == None:
										first_real_alignment_pos = b[1]
									last_real_alignment_pos = b[1]
									if b[2].isupper():
										matches_1.add(b[1])
								else:
									indel_positions.add(b[1])

							main_alignment_positions_1 = set(range(first_real_alignment_pos, last_real_alignment_pos + 1))
							read1_has_indel = len(main_alignment_positions_1.intersection(indel_positions)) > 0

							matches_2 = set([])
							first_real_alignment_pos = None
							last_real_alignment_pos = None
							indel_positions = set([])
							for b in read2_alignment.get_aligned_pairs(with_seq=True):
								if not (b[0] == None or b[1] == None):
									if first_real_alignment_pos == None:
										first_real_alignment_pos = b[1]
									last_real_alignment_pos = b[1]
									if b[2].isupper():
										matches_2.add(b[1])
								else:
									indel_positions.add(b[1])


							main_alignment_positions_2 = set(range(first_real_alignment_pos, last_real_alignment_pos + 1))
							read2_has_indel = len(main_alignment_positions_2.intersection(indel_positions)) > 0

							matches = matches_1.union(matches_2)
							main_alignment_positions = main_alignment_positions_1.union(main_alignment_positions_2)

							matching_percentage = float(len(matches))/float(len(main_alignment_positions))
							matching_percentage_read1 = float(len(matches_1))/float(len(main_alignment_positions_1))
							matching_percentage_read2 = float(len(matches_2))/float(len(main_alignment_positions_2))

							#if matching_percentage < 0.95: continue
							#if read1_has_indel or read2_has_indel or matching_percentage < 0.95 or align_overlap_prop > 0.75: continue

							read_ascpus_per_allele[read_name].append(
								[g, combined_ascore, matching_percentage, # max([matching_percentage, matching_percentage_read1, matching_percentage_read2]) ,
								 len(main_alignment_positions),
								 (read1_has_indel or read2_has_indel),
								 abs(matching_percentage_read1 - matching_percentage_read2), abs(len(main_alignment_positions_1) - len(main_alignment_positions_2))/float(len(main_alignment_positions)),
								 [read1_alignment, read2_alignment]])
							"""
							if max([matching_percentage, matching_percentage_read1, matching_percentage_read2]) == matching_percentage:
								read_ascpus_per_allele[read_name].append([g, combined_ascore, matching_percentage, len(main_alignment_positions), (read1_has_indel or read2_has_indel), abs(matching_percentage_read1-matching_percentage_read2), [read1_alignment, read2_alignment]])
							elif max([matching_percentage, matching_percentage_read1, matching_percentage_read2]) == matching_percentage_read1:
								read_ascpus_per_allele[read_name].append([g, combined_ascore, matching_percentage_read1, len(main_alignment_positions), (read1_has_indel or read2_has_indel), 0.0, [read1_alignment]])
							elif max([matching_percentage, matching_percentage_read1, matching_percentage_read2]) == matching_percentage_read2:
								read_ascpus_per_allele[read_name].append([g, combined_ascore, matching_percentage_read2, len(main_alignment_positions), (read1_has_indel or read2_has_indel), 0.0, [read2_alignment]])
							"""
						elif read1_alignment or read2_alignment:
							read_alignment = None
							if read1_alignment != None:
								read_alignment = read1_alignment
							else:
								read_alignment = read2_alignment

							read_name = read_alignment.query_name
							total_reads.add(read_name)

							matches = set([])
							first_real_alignment_pos = None
							last_real_alignment_pos = None
							indel_positions = set([])
							for b in read.get_aligned_pairs(with_seq=True):
								if b[0] != None and b[1] != None and b[2] != None:
									if first_real_alignment_pos == None:
										first_real_alignment_pos = b[1]
									last_real_alignment_pos = b[1]
									if b[2].isupper():
										matches.add(b[1])
								else:
									indel_positions.add(b[1])
							main_alignment_positions = set(range(first_real_alignment_pos, last_real_alignment_pos + 1))
							has_indel = len(main_alignment_positions.intersection(indel_positions)) > 0

							read_length = len(set(read_alignment.get_reference_positions()))
							matching_percentage = float(len(matches))/float(len(main_alignment_positions))

							# 0 added just to signify that it is a single mate contributing to the paired end combined ascore
							combined_ascore = 0 + read_alignment.tags[0][1]

							read_ascpus_per_allele[read_name].append([g, combined_ascore, matching_percentage, len(main_alignment_positions), has_indel, 0.0, 0.0, [read_alignment]])

			accounted_reads = set([])
			hg_align_pos_alleles = defaultdict(lambda: defaultdict(set))
			supported_snvs = defaultdict(lambda: defaultdict(set))
			for read in read_ascpus_per_allele:
				top_score = -1000000
				score_sorted_alignments = sorted(read_ascpus_per_allele[read], key=itemgetter(1), reverse=True)
				for i, align in enumerate(score_sorted_alignments):
					if i == 0: top_score = align[1]
					if align[1] == top_score and ((align[2] >= 0.99 and align[3] >= 100) or (align[2] >= 0.95 and align[3] >= 180)) and align[4] == False and align[5] < 0.02: # and align[6] < 0.25:
						for read_alignment in align[-1]:
							topaligns_handle.write(read_alignment)

							min_read_ref_pos = min(read_alignment.get_reference_positions())
							read_referseq = read_alignment.get_reference_sequence().upper()
							read_queryseq = read_alignment.query_sequence
							read_queryqua = read_alignment.query_qualities

							for b in read_alignment.get_aligned_pairs(with_seq=True):
								if b[0] == None or b[1] == None: continue
								ref_pos = b[1]+1
								alt_al = read_queryseq[b[0]].upper()
								ref_al = read_referseq[b[1] - min_read_ref_pos].upper()
								assert (ref_al == str(gene_sequence[align[0]]).upper()[b[1]])
								if b[2] == 'n' or ref_al == 'N' or alt_al == 'N': continue
								que_qual = read_queryqua[b[0]]
								if (que_qual >= 30) and ((ref_pos+3) < len(gene_sequence[align[0]])):
									cod_pos = gene_pos_to_msa_pos[hg][align[0]][ref_pos]
									hg_align_pos_alleles[cod_pos][alt_al].add(read)
									accounted_reads.add(read)
									det_outf.write('\t'.join([str(x) for x in [sample, hg, align[0], ref_pos, cod_pos, ref_al, alt_al, read, align[1], align[2], align[3], align[4]]]) + '\n')
									if b[2].islower():
										assert (alt_al != ref_al)
										snv_id = str(read_alignment.reference_name) + '_|_' + str(ref_pos) + '_|_' + ref_al + '_|_' + alt_al
										supported_snvs[snv_id][read].add(align[1])

			for pos in range(1, codon_alignment_lengths[hg]+1):
				printlist = [hg, str(pos)]
				for al in ['A', 'C', 'G', 'T']:
					printlist.append(str(len(hg_align_pos_alleles[pos][al])))
				res_outf.write(','.join(printlist) + '\n')

			for snv in supported_snvs:
				support_info = []
				for read in supported_snvs[snv]:
					support_info.append(read + '_|_' + str(max(supported_snvs[snv][read])))
				snv_outf.write('\t'.join([snv, str(len(support_info))] + support_info) + '\n')

			#print('\t'.join([hg, str(len(total_reads)), str(len(accounted_reads))]))

		snv_outf.close()
		res_outf.close()
		det_outf.close()
		topaligns_handle.close()
		bam_handle.close()

		os.system("samtools sort -@ %d %s -o %s" % (1, topaligns_file, topaligns_file_sorted))
		os.system("samtools index %s" % topaligns_file_sorted)

	except Exception as e:
		if logObject:
			logObject.error("Issues with mining for SNVs/parsing alignment file.")
			logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

def popgen_analysis_of_hg(inputs):
	"""
	Helper function which is to be called from the runPopulationGeneticsAnalysis() function to parallelize population
	genetics analysis of each homolog group.

	:param inputs: list of inputs passed in by GCF.runPopulationGeneticsAnalysis().
	"""

	gcf_id, gcf_annot, hg, codon_alignment_fasta, popgen_dir, plots_dir, comp_gene_info, hg_genes, bgc_sample, hg_prop_multi_copy, hg_order_scpus, gw_pairwise_similarities, comparem_used, sample_population, population, species_phylogeny, sample_size, logObject = inputs

	domain_plot_file = plots_dir + hg + '_domain.txt'
	position_plot_file = plots_dir + hg + '_position.txt'
	plot_pdf_file = plots_dir + hg + '.pdf'

	domain_plot_handle = open(domain_plot_file, 'w')
	position_plot_handle = open(position_plot_file, 'w')

	seqs = []
	samples_ordered = []
	genes_ordered = []
	bgc_codons = defaultdict(list)
	num_codons = None
	samples = set([])
	gene_lengths = []
	gene_locs = defaultdict(dict)
	core_counts = defaultdict(int)
	products = set([])
	sample_leaf_names = defaultdict(list)
	updated_codon_alignment_fasta = popgen_dir + codon_alignment_fasta.split('/')[-1]
	updated_codon_alignment_handle = open(updated_codon_alignment_fasta, 'w')
	with open(codon_alignment_fasta) as ocaf:
		for rec in SeqIO.parse(ocaf, 'fasta'):
			sample_id, gene_id = rec.id.split('|')
			if sample_population != None and population != None and population != sample_population[sample_id]: continue
			if len(gene_id.split('_')[0]) == 3:
				if comp_gene_info[gene_id]['core_overlap']:
					core_counts['core'] += 1
				else:
					core_counts['auxiliary'] += 1
			updated_codon_alignment_handle.write('>' + rec.description + '\n' + str(rec.seq) + '\n')
			products.add(comp_gene_info[gene_id]['product'])
			real_pos = 1
			seq = str(rec.seq).upper().replace('N', '-')
			#seqlen = len(seq)
			#gapless_seqlen = len(b for b in seqlen if b != '-')
			sample_leaf_names[sample_id].append(rec.id)
			seqs.append(list(seq))
			codons = [str(rec.seq)[i:i + 3] for i in range(0, len(str(rec.seq)), 3)]
			num_codons = len(codons)
			bgc_codons[rec.id] = codons
			samples.add(sample_id)
			samples_ordered.append(sample_id)
			genes_ordered.append(gene_id)
			for msa_pos, bp in enumerate(str(rec.seq)):
				if bp != '-':
					gene_locs[gene_id][real_pos] = msa_pos + 1
					real_pos += 1
			gene_lengths.append(len(str(rec.seq).replace('-', '')))
	updated_codon_alignment_handle.close()
	codon_alignment_fasta = updated_codon_alignment_fasta

	if len(seqs) == 0: return

	single_copy_in_gcf = True
	for sample in sample_leaf_names:
		if len(sample_leaf_names[sample]) > 1:
			single_copy_in_gcf = False

	median_beta_rd = "NA"
	if comparem_used != None:
		if gw_pairwise_similarities:
			beta_rd_stats = []
			hg_pairwise_similarities = util.determineSeqSimCodonAlignment(codon_alignment_fasta, use_translation=comparem_used)
			for i, s1 in enumerate(sorted(samples)):
				for j, s2 in enumerate(sorted(samples)):
					if i >= j: continue
					if s1 in hg_pairwise_similarities and s2 in hg_pairwise_similarities:
						hg_seq_sim = hg_pairwise_similarities[s1][s2]
						gw_seq_sim = gw_pairwise_similarities[s1][s2]
						if gw_seq_sim != 0.0:
							beta_rd = hg_seq_sim / float(gw_seq_sim)
							beta_rd_stats.append(beta_rd)
			median_beta_rd = "NA"
			max_beta_rd = "NA"
			if len(beta_rd_stats) >= 1:
				median_beta_rd = round(statistics.median(beta_rd_stats), 2)
				max_beta_rd = round(max(beta_rd_stats), 2)

	is_core = False
	if (sum(core_counts.values()) > 0.0):
		if (float(core_counts['core']) / sum(core_counts.values()) >= 0.5):
			is_core = True

	median_gene_length = statistics.median(gene_lengths)

	variable_sites = set([])
	conserved_sites = set([])
	nondominant_sites  =  set([])
	position_plot_handle.write('\t'.join(['pos', 'num_seqs', 'num_alleles', 'num_gaps', 'maj_allele_freq']) + '\n')

	# TODO: consider out-souring filtering to use phykit
	sample_differences_to_consensus = defaultdict(lambda: defaultdict(int))
	ambiguous_sites = 0
	nonambiguous_sites = 0
	ambiguous_sites_pos = set([])
	for i, ls in enumerate(zip(*seqs)):
		al_counts = defaultdict(int)
		for al in ls:
			al_counts[al] += 1
		maj_allele_count = 0
		if not (len(al_counts) == 1 and '-' in al_counts):
			maj_allele_count = max([al_counts[al] for al in al_counts if al != '-'])
		tot_count = sum(al_counts.values())
		num_alleles = len(al_counts.keys())
		num_gaps = al_counts['-']
		if num_gaps > 0:
			num_alleles -= 1
		gap_allele_freq = float(num_gaps) / tot_count
		maj_allele_freq = 0.0
		if float(tot_count-num_gaps) > 0.0:
			maj_allele_freq = float(maj_allele_count) / float(tot_count-num_gaps)
		position_plot_handle.write('\t'.join([str(x) for x in [i + 1, tot_count, num_alleles, num_gaps, maj_allele_freq]]) + '\n')
		if gap_allele_freq < 0.10:
			if maj_allele_freq >= 0.95:
				conserved_sites.add(i)
			else:
				variable_sites.add(i)
			if maj_allele_freq < 0.75:
				nondominant_sites.add(i)
			nonambiguous_sites += 1
		else:
			ambiguous_sites += 1
			ambiguous_sites_pos.add(i+1)
		maj_allele_count = max(al_counts.values())
		maj_alleles = set([a[0] for a in al_counts.items() if maj_allele_count == a[1]])
		maj_allele = sorted(list(maj_alleles))[0]
		for j, al in enumerate(ls):
			sid = samples_ordered[j]
			gid = genes_ordered[j]
			if al != maj_allele:
				sample_differences_to_consensus[sid][gid] += 1
			else:
				sample_differences_to_consensus[sid][gid] += 0
	position_plot_handle.close()

	ambiguous_prop = float(ambiguous_sites)/float(ambiguous_sites + nonambiguous_sites)

	# avoid looking at domains if
	domain_positions_msa = defaultdict(set)
	domain_min_position_msa = defaultdict(lambda: 1e8)
	all_domains = set([])
	domain_plot_handle.write('\t'.join(['domain', 'domain_index', 'min_pos', 'max_pos']) + '\n')

	issue_with_domain_coords = False
	at_least_one_gene_is_multi_part = False
	at_least_one_domain_is_multi_part = False
	for gene in gene_locs:
		gene_is_multi_part = comp_gene_info[gene]['is_multi_part']
		if gene_is_multi_part or at_least_one_domain_is_multi_part:
			at_least_one_gene_is_multi_part = True
			break
		gene_start = comp_gene_info[gene]['start']
		gene_end = comp_gene_info[gene]['end']
		gene_direction = comp_gene_info[gene]['direction']
		for domain_info in comp_gene_info[gene]['gene_domains']:
			dom_is_multi_part = domain_info['is_multi_part']
			if dom_is_multi_part:
				at_least_one_domain_is_multi_part = True
				break
			if gene_direction == '-':
				dom_start = domain_info['start']
				dom_end = domain_info['end']
				dom_length = dom_end - dom_start
				dist_to_dom_start = max(dom_start - gene_start, 0)
				flip_dom_end = gene_end - dist_to_dom_start
				flip_dom_start = flip_dom_end - dom_length
				domain_info['start'] = flip_dom_start
				domain_info['end'] = flip_dom_end
			domain_start = max(domain_info['start'], gene_start)
			domain_end = min(domain_info['end'], gene_end)
			domain_name = domain_info['aSDomain'] + '_|_' + domain_info['description']

			relative_start = domain_start - gene_start
			try:
				assert (len(gene_locs[gene]) + 3 >= (domain_end - gene_start))
			except:
				issue_with_domain_coords = True
			relative_end = min([len(gene_locs[gene]), domain_end - gene_start])
			domain_range = range(relative_start, relative_end)
			for pos in domain_range:
				msa_pos = gene_locs[gene][pos + 1]
				if domain_info['type'] == 'PFAM_domain':
					domain_positions_msa[domain_name].add(msa_pos)
					if domain_min_position_msa[domain_name] > msa_pos:
						domain_min_position_msa[domain_name] = msa_pos
			all_domains.add(domain_info['type'] + '_|_' + domain_info['aSDomain'] + '_|_' + domain_info['description'])

	if not at_least_one_domain_is_multi_part and not at_least_one_gene_is_multi_part:
		if issue_with_domain_coords:
			raise RuntimeError("Issue with fitting domain coordinates to gene coordinates")
		for i, dom in enumerate(sorted(domain_min_position_msa.items(), key=itemgetter(1))):
			tmp = []
			old_pos = None
			for j, pos in enumerate(sorted(domain_positions_msa[dom[0]])):
				if j == 0:
					old_pos = pos - 1
				if pos - 1 != old_pos:
					if len(tmp) > 0:
						min_pos = min(tmp)
						max_pos = max(tmp)
						domain_plot_handle.write('\t'.join([str(x) for x in [dom[0], i, min_pos, max_pos]]) + '\n')
					tmp = []
				tmp.append(pos)
				old_pos = pos
			if len(tmp) > 0:
				min_pos = min(tmp)
				max_pos = max(tmp)
				domain_plot_handle.write('\t'.join([str(x) for x in [dom[0], i, min_pos, max_pos]]) + '\n')

	domain_plot_handle.close()

	rscript_plot_cmd = ["Rscript", RSCRIPT_FOR_CLUSTER_ASSESSMENT_PLOTTING, domain_plot_file, position_plot_file,
						plot_pdf_file]
	if logObject:
		logObject.info('Running R-based plotting with the following command: %s' % ' '.join(rscript_plot_cmd))
	try:
		subprocess.call(' '.join(rscript_plot_cmd), shell=True, stdout=subprocess.DEVNULL,
						stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		if logObject:
			logObject.info('Successfully ran: %s' % ' '.join(rscript_plot_cmd))
	except Exception as e:
		if logObject:
			logObject.error('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

	#### Perform population genetics analyses, including dN/dS calculation and Tajima's D calculation

	high_ambiguity_sequences = set([])
	with open(codon_alignment_fasta) as ocaf:
		for rec in SeqIO.parse(ocaf, 'fasta'):
			if sample_population != None and population != None and population != sample_population[sample_id]: continue
			total_nonambiguous_positions = 0
			gap_nonambiguous_positions = 0
			for msa_pos, bp in enumerate(str(rec.seq)):
				if not (msa_pos+1) in ambiguous_sites_pos:
					total_nonambiguous_positions += 1
					if bp == '-':
						gap_nonambiguous_positions += 1
			if total_nonambiguous_positions > 0:
				seq_ambiguous_prop = float(gap_nonambiguous_positions)/float(total_nonambiguous_positions)
				if seq_ambiguous_prop >= 0.25:
					high_ambiguity_sequences.add(rec.id)
			else:
				high_ambiguity_sequences.add(rec.id)

	sequences_filtered = defaultdict(lambda: '')
	for cod_index in range(0, num_codons):
		cod_count = defaultdict(int)
		for bgc in bgc_codons:
			if bgc in high_ambiguity_sequences: continue
			cod = bgc_codons[bgc][cod_index].replace('N', '-')
			if not '-' in cod:
				cod_count[cod] += 1

		#singleton_codons = set([])
		#for cod in cod_count:
		#	if cod_count[cod] == 1: singleton_codons.add(cod)

		aa_count = defaultdict(int)
		cod_count = defaultdict(int)
		cods = set([])
		for bgc in bgc_codons:
			if bgc in high_ambiguity_sequences: continue
			cod = bgc_codons[bgc][cod_index].replace('N', '-')
			if '-' in cod: cod = '---'
			aa = None
			if '-' in cod: aa = '-'
			else: cod_obj = Seq(cod); aa = str(cod_obj.translate())
			cods.add(cod)
			cod_count[cod] += 1
			aa_count[aa] += 1
		if len(cod_count) > 0:
			major_codon_count = max(cod_count.values())
			#major_codon_freq = float(major_codon_count) / float(sum(cod_count.values()))
			gap_residue_freq = float(aa_count['-']) / float(sum(aa_count.values()))
			if ((len(cod_count) >= 3) or (len(cod_count) >= 2 and gap_residue_freq == 0.0)) and gap_residue_freq < 0.1:
				for bgc in bgc_codons:
					if bgc in high_ambiguity_sequences: continue
					cod = bgc_codons[bgc][cod_index].replace('N', '-')
					if '-' in cod: cod = '---'
					sequences_filtered[bgc] += cod


	median_dnds = "NA"
	mad_dnds = "NA"
	tajimas_d = "NA"
	if len(sequences_filtered) >= 4 and len(list(sequences_filtered.values())[0]) >= 21:
		if species_phylogeny:
			pass
			"""
			# If species phylogeny was provided we can perform dn/ds analysis using ancestral recosntructed sequence for 
			# the homolog group
	
			# First, we filter the phylogeny for samples with the HG and create duplicate leafs for samples with multiple
			# instances of the HG identified in BGCs belonging to the GCF.
			"""
		else:


			all_median_dnds = []
			for iter in range(0, 20):
				combos = list(itertools.combinations(list(sequences_filtered.values()), 2))
				random.Random(iter).shuffle(combos)

				all_dNdS = []
				for i, pair in enumerate(combos):
					if i >= sample_size: continue
					seqA = pair[0]
					seqB = pair[1]
					assert(len(seqA) == len(seqB))
					csA = CodonSeq(seqA)
					csB = CodonSeq(seqB)
					try:
						dN, dS = cal_dn_ds(csA, csB)
						if dN != -1 and dS != -1 and dS != 0.0:
							all_dNdS.append(float(dN)/float(dS))
					except:
						# ignore errors. Introduced because of an error caused by no synonymous sites being identified
						# between a pair of sequences resulting in Division by Zero error being raised.
						pass
				if len(all_dNdS) > 0:
					all_median_dnds.append(statistics.median(all_dNdS))

			if len(all_median_dnds) >= 10:
				median_dnds = statistics.median(all_median_dnds)
				mad_dnds = median_abs_deviation(all_median_dnds, scale="normal")

		# calculate Tajima's D
		tajimas_d = util.calculateTajimasD(list(sequences_filtered.values()))
		if tajimas_d != 'NA':
			tajimas_d = round(tajimas_d, 2)

	prop_samples_with_hg = len(samples) / float(len(set(bgc_sample.values())))
	prop_conserved = "NA"
	prop_majallele_nondominant = "NA"
	if (len(conserved_sites) + len(variable_sites)) > 0:
		prop_conserved = round(float(len(variable_sites))/float(len(conserved_sites) + len(variable_sites)), 2)
		prop_majallele_nondominant = round(float(len(nondominant_sites))/float(len(conserved_sites) + len(variable_sites)), 2)
	else:
		# All sites in alignment have >= 10% ambiguity
		prop_conserved = "NA"
		prop_majallele_nondominant = "NA"

	hg_ord = 'NA'; hg_dir = 'NA'
	if hg in hg_order_scpus:
		hg_ord, hg_dir = hg_order_scpus[hg]

	hg_info = []
	if population:
		hg_info = [population]
	hg_info += [gcf_id, gcf_annot, hg, '; '.join(products), hg_ord, hg_dir, hg_prop_multi_copy[hg], single_copy_in_gcf,
				median_gene_length, is_core, len(seqs), len(samples), round(prop_samples_with_hg,2),
				round(ambiguous_prop,2), tajimas_d, prop_conserved, prop_majallele_nondominant, median_beta_rd,
				max_beta_rd, median_dnds, mad_dnds]

	hg_consim_handle = open(popgen_dir + hg + '_sim_to_consensus.txt', 'w')
	for s in sample_differences_to_consensus:
		min_diff_to_consensus = 1e100
		for g in sample_differences_to_consensus[s]:
			if min_diff_to_consensus > sample_differences_to_consensus[s][g]:
				min_diff_to_consensus = sample_differences_to_consensus[s][g]
		if min_diff_to_consensus < 1e100:
			hg_consim_handle.write(hg + '\t' + s + '\t' + str(min_diff_to_consensus / float(len(seqs[0]))) + '\n')
	hg_consim_handle.close()

	if sample_population and not population:
		population_counts = defaultdict(int)
		population_samples = defaultdict(set)
		for s, p in sample_population.items():
			population_counts[p] += 1
			population_samples[p].add(s)

		most_positive_tajimas_d = [['NA'], 0.0]
		most_negative_tajimas_d = [['NA'], 0.0]
		all_tajimas_d = []

		for p in population_counts:
			if population_counts[p] < 4: continue
			population_sequences = []
			for seq_id in sequences_filtered:
				if sample_population[seq_id.split('|')[0]] == p:
					population_sequences.append(sequences_filtered[seq_id])
			if len(population_sequences) >= 4 and len(list(sequences_filtered.values())[0]) >= 21:
				p_tajimas_d = util.calculateTajimasD(population_sequences)
				if p_tajimas_d != 'NA':
					p_tajimas_d = round(p_tajimas_d, 3)
					all_tajimas_d.append(p_tajimas_d)
					if p_tajimas_d > most_positive_tajimas_d[1]:
						most_positive_tajimas_d = [[p], p_tajimas_d]
					elif p_tajimas_d == most_positive_tajimas_d[1]:
						most_positive_tajimas_d[0].append(p)

					if p_tajimas_d < most_negative_tajimas_d[1]:
						most_negative_tajimas_d = [[p], p_tajimas_d]
					elif p_tajimas_d == most_negative_tajimas_d[1]:
						most_negative_tajimas_d[0].append(p)

		median_tajimas_d = "NA"
		mad_tajimas_d = "NA"
		if len(all_tajimas_d) > 0:
			median_tajimas_d = statistics.median(all_tajimas_d)
			mad_tajimas_d = median_abs_deviation(all_tajimas_d, scale="normal")

		pops_with_hg = set([])
		pop_count_with_hg = defaultdict(int)
		for s in sample_differences_to_consensus:
			min_diff_to_consensus = 1e100
			for g in sample_differences_to_consensus[s]:
				if min_diff_to_consensus > sample_differences_to_consensus[s][g]:
					min_diff_to_consensus = sample_differences_to_consensus[s][g]
			if min_diff_to_consensus < 1e100:
				pop_count_with_hg[sample_population[s]] += 1
				pops_with_hg.add(sample_population[s])

		fishers_pvals = []
		fst_like_estimates = []
		for pi, pop in enumerate(pop_count_with_hg):
			other_count = sum([x[1] for x in pop_count_with_hg.items() if x[0] != pop])
			other_total = sum([x[1] for x in population_counts.items() if x[0] != pop])
			odds, pval = fisher_exact([[pop_count_with_hg[pop], population_counts[pop]], [other_count, other_total]])
			fishers_pvals.append(pval)
			pds_within = []
			for si1, s1 in enumerate(sorted(population_samples[pop])):
				for si2, s2 in enumerate(sorted(population_samples[pop])):
					if si1 >= si2: continue
					min_diff = 1e100
					for g1 in sample_differences_to_consensus[s1]:
						for g2 in sample_differences_to_consensus[s2]:
							samp_gene_diff = abs(sample_differences_to_consensus[s1][g1] - sample_differences_to_consensus[s2][g2])
							if samp_gene_diff < min_diff:
								min_diff = samp_gene_diff
					if min_diff < 1e100:
						pds_within.append(min_diff)
			pi_within = 0.0
			if len(pds_within) > 0:
				pi_within = statistics.median(pds_within)
			for cpi, comp_pop in enumerate(pop_count_with_hg):
				if pop != comp_pop:
					pds_between = []
					for si1, s1 in enumerate(sorted(population_samples[pop])):
						for si2, s2 in enumerate(sorted(population_samples[comp_pop])):
							min_diff = 1e100
							for g1 in sample_differences_to_consensus[s1]:
								for g2 in sample_differences_to_consensus[s2]:
									samp_gene_diff = abs(sample_differences_to_consensus[s1][g1] - sample_differences_to_consensus[s2][g2])
									if samp_gene_diff < min_diff:
										min_diff = samp_gene_diff
							if min_diff < 1e100:
								pds_between.append(min_diff)
					if len(pds_between) > 0:
						pi_between = statistics.median(pds_between)
						if pi_between > 0:
							fst_like_est = float(pi_between - pi_within)/float(pi_between)
							fst_like_estimates.append(fst_like_est)
						else:
							fst_like_estimates.append(0.0)

		median_fst_like_est = "NA"
		if len(fst_like_estimates) > 0:
			median_fst_like_est = statistics.median(fst_like_estimates)

		fisher_pval = "NA"
		if len(fishers_pvals) > 0:
			fisher_pval = min(fishers_pvals)*len(fishers_pvals)

		pop_rel_freqs = []
		for p in population_counts:
			prf = float(pop_count_with_hg[p])/float(sum(pop_count_with_hg.values()))
			pop_rel_freqs.append(prf)

		population_entropy = entropy(pop_rel_freqs, base=2)

		most_neg_taj_d = most_negative_tajimas_d[1]
		if most_negative_tajimas_d[0][0] == 'NA': most_neg_taj_d = 'NA'

		most_pos_taj_d = most_positive_tajimas_d[1]
		if most_positive_tajimas_d[0][0] == 'NA': most_pos_taj_d = 'NA'

		hg_population_info = [len(pops_with_hg), float(len(pops_with_hg))/float(len(population_counts)), fisher_pval,
							  median_tajimas_d, mad_tajimas_d, str(most_neg_taj_d) + '|' + ','.join(most_negative_tajimas_d[0]),
							  str(most_pos_taj_d) + '|' + ','.join(most_positive_tajimas_d[0]),
							  population_entropy, median_fst_like_est, '|'.join([str(x[0]) + '=' + str(float(x[1])/population_counts[x[0]]) for x in pop_count_with_hg.items()])]
		hg_info += hg_population_info

	hg_info += ['; '.join(all_domains)]
	hg_stats_handle = open(popgen_dir + hg + '_stats.txt', 'w')
	hg_stats_handle.write('\t'.join([str(x) for x in hg_info]) + '\n')
	hg_stats_handle.close()

def create_codon_msas(inputs):
	"""
	Helper function which is to be called from the constructCodonAlignments() function to parallelize construction
	of codon alignments for each homolog group of interest in the GCF.
	:param inputs: list of inputs passed in by GCF.constructCodonAlignments().
	"""
	hg, gene_sequences, nucl_seq_dir, prot_seq_dir, prot_alg_dir, codo_alg_dir, cpus, use_magus, logObject = inputs

	hg_nucl_fasta = nucl_seq_dir + '/' + hg + '.fna'
	hg_prot_fasta = prot_seq_dir + '/' + hg + '.faa'
	hg_prot_msa = prot_alg_dir + '/' + hg + '.msa.faa'
	hg_codo_msa = codo_alg_dir + '/' + hg + '.msa.fna'

	hg_nucl_handle = open(hg_nucl_fasta, 'w')
	hg_prot_handle = open(hg_prot_fasta, 'w')
	# sorted added because order of sequences might be influencing MAFFT alignment
	# TODO: confirm this makes MAFFT results reproducible, no option for seeding apparent
	for s in sorted(gene_sequences):
		hg_nucl_handle.write('>' + s + '\n' + str(gene_sequences[s][0]) + '\n')
		hg_prot_handle.write('>' + s + '\n' + str(gene_sequences[s][1]) + '\n')
	hg_nucl_handle.close()
	hg_prot_handle.close()

	mafft_cmd = ['mafft', '--thread', str(cpus), '--maxiterate', '1000', '--localpair', hg_prot_fasta, '>', hg_prot_msa]
	rm_cmd = None
	if use_magus:
		tmp_outdir = prot_alg_dir + 'tmp_' + hg + '/'
		mafft_cmd = ['magus', '-d', tmp_outdir, '-i', hg_prot_fasta, '-o', hg_prot_msa, '-np', str(cpus), '-r', '1']
		rm_cmd = ['rm', '-r', tmp_outdir]
	pal2nal_cmd = ['pal2nal.pl', hg_prot_msa, hg_nucl_fasta, '-output', 'fasta', '>', hg_codo_msa]

	if logObject:
		logObject.info('Running mafft/magus with the following command: %s' % ' '.join(mafft_cmd))
	try:
		subprocess.call(' '.join(mafft_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		if logObject:
			logObject.info('Successfully ran: %s' % ' '.join(mafft_cmd))
	except Exception as e:
		if logObject:
			logObject.error('Had an issue running: %s' % ' '.join(mafft_cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(mafft_cmd))

	if logObject:
		logObject.info('Running PAL2NAL with the following command: %s' % ' '.join(pal2nal_cmd))
	try:
		subprocess.call(' '.join(pal2nal_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		if logObject:
			logObject.info('Successfully ran: %s' % ' '.join(pal2nal_cmd))
	except Exception as e:
		if logObject:
			logObject.error('Had an issue running: %s' % ' '.join(pal2nal_cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(pal2nal_cmd))

	if rm_cmd != None:
		logObject.info("Removing tmp output directory for MAGUS.")
		os.system(' '.join(rm_cmd))
	if logObject:
		logObject.info('Achieved codon alignment for homolog group %s' % hg)

def identify_gcf_instances(input_args):
	bgc_info_dir, bgc_genbanks_dir, sample, sample_prokka_data, sample_lt_to_evalue, hmmscan_results_lenient, model, lts_ordered_dict, hgs_ordered_dict, comp_gene_info, gene_location, gene_id_to_order, gene_order_to_id, protocluster_core_homologs, core_homologs, boundary_genes, specific_hgs, bgc_genes, gene_to_hg, min_size, min_core_size, surround_gene_max, syntenic_correlation_threshold, loose_flag = input_args

	sample_gcf_predictions = []
	sample_bgc_ids = 1
	for scaffold in hgs_ordered_dict:
		hgs_ordered = hgs_ordered_dict[scaffold]
		lts_ordered = lts_ordered_dict[scaffold]
		hg_seq = numpy.array(list(hgs_ordered))
		hmm_predictions = model.predict(hg_seq)

		gcf_state_lts = []
		gcf_state_hgs = []
		for i, hg_state in enumerate(hmm_predictions):
			lt = lts_ordered[i]
			hg = hgs_ordered[i]
			if hg_state == 0:
				gcf_state_lts.append(lt)
				gcf_state_hgs.append(hg)
			if hg_state == 1 or i == (len(hmm_predictions) - 1):
				if len(set(gcf_state_hgs).difference("other")) >= 3:
					boundary_lt_featured = False
					features_specific_hg = False
					features_protocluster_hg = False
					if len(protocluster_core_homologs.intersection(set(gcf_state_hgs).difference('other'))) > 0: features_protocluster_hg = True
					if len(boundary_genes.intersection(set(gcf_state_lts).difference('other'))) > 0: boundary_lt_featured = True
					if len(specific_hgs.intersection(set(gcf_state_hgs).difference('other'))) > 0: features_specific_hg = True
					sample_gcf_predictions.append([gcf_state_lts, gcf_state_hgs, len(gcf_state_lts),
												   len(set(gcf_state_hgs).difference("other")),
												   len(set(gcf_state_hgs).difference("other").intersection(core_homologs)), scaffold, boundary_lt_featured,
												   features_specific_hg, features_protocluster_hg])
				gcf_state_lts = []
				gcf_state_hgs = []

	if len(sample_gcf_predictions) == 0: return

	sorted_sample_gcf_predictions = [x for x in sorted(sample_gcf_predictions, key=itemgetter(3), reverse=True)]

	cumulative_edge_hgs = set([])
	visited_scaffolds_with_edge_gcf_segment = set([])
	sample_gcf_predictions_filtered = []
	sample_edge_gcf_predictions_filtered = []

	for gcf_segment in sorted_sample_gcf_predictions:
		if (gcf_segment[3] >= min_size and gcf_segment[4] >= min_core_size) or (gcf_segment[-1]) or (gcf_segment[-2]) or (gcf_segment[3] >= 3 and gcf_segment[-3] and not gcf_segment[5] in visited_scaffolds_with_edge_gcf_segment):
			# code to determine whether syntenically, the considered segment aligns with what is expected.
			segment_hg_order = []
			segment_hg_direction = []
			bgc_hg_orders = defaultdict(list)
			bgc_hg_directions = defaultdict(list)

			copy_count_of_hgs_in_segment = defaultdict(int)
			for hg in gcf_segment[1]:
				copy_count_of_hgs_in_segment[hg] += 1

			for gi, g in enumerate(gcf_segment[0]):
				hg = gcf_segment[1][gi]
				if copy_count_of_hgs_in_segment[hg] != 1: continue
				gene_midpoint = (gene_location[g]['start'] + gene_location[g]['end']) / 2.0
				segment_hg_order.append(gene_midpoint)
				segment_hg_direction.append(gene_location[g]['direction'])

				for bgc in bgc_genes:
					bg_matching = []
					for bg in bgc_genes[bgc]:
						if bg in gene_to_hg:
							hg_of_bg = gene_to_hg[bg]
							if hg_of_bg == hg: bg_matching.append(bg)
					if len(bg_matching) == 1:
						bgc_gene_midpoint = (comp_gene_info[bg_matching[0]]['start'] + comp_gene_info[bg_matching[0]]['end']) / 2.0
						bgc_hg_orders[bgc].append(bgc_gene_midpoint)
						bgc_hg_directions[bgc].append(comp_gene_info[bg_matching[0]]['direction'])
					else:
						bgc_hg_orders[bgc].append(None)
						bgc_hg_directions[bgc].append(None)

			best_corr = None
			for bgc in bgc_genes:
				try:
					assert (len(segment_hg_order) == len(bgc_hg_orders[bgc]))

					list1_same_dir = []
					list2_same_dir = []
					list1_comp_dir = []
					list2_comp_dir = []
					for iter, hgval1 in enumerate(segment_hg_order):
						hgdir1 = segment_hg_direction[iter]
						hgval2 = bgc_hg_orders[bgc][iter]
						hgdir2 = bgc_hg_directions[bgc][iter]
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

			if not best_corr or best_corr < syntenic_correlation_threshold: continue

			if (gcf_segment[3] >= min_size and gcf_segment[4] >= min_core_size) or (gcf_segment[-1]) or (gcf_segment[-2]):
				sample_gcf_predictions_filtered.append(gcf_segment)
				if gcf_segment[-3]:
					cumulative_edge_hgs = cumulative_edge_hgs.union(set(gcf_segment[1]))
					visited_scaffolds_with_edge_gcf_segment.add(gcf_segment[5])
			elif gcf_segment[3] >= 3 and gcf_segment[-3] and not gcf_segment[5] in visited_scaffolds_with_edge_gcf_segment:
				sample_edge_gcf_predictions_filtered.append(gcf_segment)
				visited_scaffolds_with_edge_gcf_segment.add(gcf_segment[5])
				cumulative_edge_hgs = cumulative_edge_hgs.union(set(gcf_segment[1]))

	if len(sample_edge_gcf_predictions_filtered) >= 1:
		if len(cumulative_edge_hgs) >= min_size and len(cumulative_edge_hgs.intersection(core_homologs)) >= min_core_size:
			sample_gcf_predictions_filtered += sample_edge_gcf_predictions_filtered

	if not loose_flag:
		protocore_gene_found = False
		for gcf_segment in sample_gcf_predictions_filtered:
			if gcf_segment[-1]: protocore_gene_found = True
		if not protocore_gene_found: return

	bgc_sample_listing_handle = open(bgc_info_dir + sample + '.bgcs.txt', 'w')
	bgc_hg_evalue_handle = open(bgc_info_dir + sample + '.hg_evalues.txt', 'w')

	for gcf_segment in sample_gcf_predictions_filtered:
		bgc_genbank_file = bgc_genbanks_dir + sample + '_Expansion_BGC-' + str(sample_bgc_ids) + '.gbk'
		sample_bgc_ids += 1

		gcf_segment_scaff = gcf_segment[5]
		# check if you can expand and name more hgs
		for i, lt in enumerate(gcf_segment[0]):
			hg = gcf_segment[1][i]
			if hg == 'other' and lt in hmmscan_results_lenient.keys():
				gcf_segment[1][i] = hmmscan_results_lenient[lt][0]

		min_bgc_order = min([gene_id_to_order[gcf_segment_scaff][g] for g in gcf_segment[0]])
		max_bgc_order = max([gene_id_to_order[gcf_segment_scaff][g] for g in gcf_segment[0]])

		for oi in range(min_bgc_order-surround_gene_max, min_bgc_order):
			if oi in gene_order_to_id[gcf_segment_scaff].keys():
				lt = gene_order_to_id[gcf_segment_scaff][oi]
				if lt in hmmscan_results_lenient.keys():
					gcf_segment[0].append(lt)
					gcf_segment[1].append(hmmscan_results_lenient[lt][0])

		for oi in range(max_bgc_order+1, max_bgc_order+surround_gene_max+1):
			if oi in gene_order_to_id[gcf_segment_scaff].keys():
				lt = gene_order_to_id[gcf_segment_scaff][oi]
				if lt in hmmscan_results_lenient.keys():
					gcf_segment[0].append(lt)
					gcf_segment[1].append(hmmscan_results_lenient[lt][0])

		min_bgc_pos = min([gene_location[g]['start'] for g in gcf_segment[0]])
		max_bgc_pos = max([gene_location[g]['end'] for g in gcf_segment[0]])

		util.createBGCGenbank(sample_prokka_data['genbank'], bgc_genbank_file, gcf_segment_scaff,
							  min_bgc_pos, max_bgc_pos)
		bgc_sample_listing_handle.write('\t'.join([sample, bgc_genbank_file]) + '\n')

		for i, lt in enumerate(gcf_segment[0]):
			hg = gcf_segment[1][i]
			evalue = decimal.Decimal(100000.0)
			if lt in sample_lt_to_evalue: evalue = sample_lt_to_evalue[lt]
			elif lt in hmmscan_results_lenient.keys(): evalue = hmmscan_results_lenient[lt][1]
			bgc_hg_evalue_handle.write('\t'.join([bgc_genbank_file, sample, lt, hg, str(evalue), str(hg in protocluster_core_homologs)]) + '\n')

	bgc_hg_evalue_handle.close()
	bgc_sample_listing_handle.close()
