import os
import sys
import logging
import traceback
import subprocess
from collections import defaultdict
from Bio import SeqIO
from zol.classes.BGC import BGC
from zol import util
import statistics
import random
import multiprocessing
import math
from operator import itemgetter
from scipy import stats
import decimal
import _pickle as cPickle

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-3])
RSCRIPT_FOR_CLUSTER_ASSESSMENT_PLOTTING = lsaBGC_main_directory + '/zol/Rscripts/plotParameterImpactsOnGCF.R'

class Pan:
	def __init__(self, bgc_genbanks_listing, logObject=None, lineage_name='Unnamed lineage'):
		self.bgc_genbanks_listing = bgc_genbanks_listing
		self.lineage_name = lineage_name
		self.logObject = logObject

		#######
		## Variables not set during initialization
		#######

		# General variables
		self.pan_bgcs = {}
		self.comp_gene_info = {}
		self.bgc_info = {}
		self.bgc_gbk = {}
		self.bgc_genes = {}
		self.pan_genes = set([])
		self.bgc_sample = {}
		self.sample_bgcs = defaultdict(set)
		self.bgc_product = {}
		self.bgc_core_counts = {}
		self.bgc_population = None
		self.sample_population = None

		# homology related variables
		self.gene_to_hg = None
		self.hg_genes = None
		self.hg_median_copy_count = None
		self.hg_prop_multi_copy = None
		self.bgc_hgs = defaultdict(set)

		### other variables
		# Jaccard similarity distance between bgcs
		self.pairwise_relations = None
		# Absolute pearson correlation coefficient for syntenic similarity between bgcs
		self.pairwise_syntenic_relations = None
		# Concatenated HMMER3 HMM profiles database of homolog groups in GCF
		self.concatenated_profile_HMM = None
		# Consensus sequence in fasta format of concatentated profile HMM file
		self.consensus_sequence_HMM = None
		self.consensus_sequence_HMM_db = None
		self.hg_max_self_evalue = defaultdict(lambda: [100000.0, False])
		self.lowerbound_hg_count = 100000
		self.hg_differentiation_stats = {}

		# HMMScan results - dictionary with keys as gene locus tag ids and values as
		# list of lists: [homolog group, evalue, sample, scaffold]
		self.hmmscan_results = defaultdict(list)
		self.hmmscan_results_lenient = defaultdict(dict)
		self.boundary_genes = defaultdict(set)
		self.scaffold_genes = defaultdict(lambda: defaultdict(set))
		self.gene_location = defaultdict(dict)
		self.gene_id_to_order = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
		self.gene_order_to_id = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

		# variables containing location of files
		self.final_stats_file = None
		self.final_stats_expanded_singletons_file = None
		self.pair_relations_txt_file = None
		self.bgc_to_gcf_map_file = None

	def readInBGCGenbanks(self, comprehensive_parsing=True, prune_set=None, prediction_method='ANTISMASH'):
		"""
		Function to parse file listing location of BGC Genbanks.

		:param comprehensive_parsing (optional): flag specifying whether to perform comprehensive extraction of information from Genbanks. default is True.
		"""
		sample_index = defaultdict(int)
		with open(self.bgc_genbanks_listing) as obsf:
			for i, line in enumerate(obsf):
				line = line.strip()
				try:
					assert (len(line.split('\t')) == 2)
				except Exception as e:
					if self.logObject:
						self.logObject.error("More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (i + 1))
						self.logObject.error(traceback.format_exc())
					raise RuntimeError("More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (i + 1))
				sample, gbk = line.split('\t')
				if prune_set != None and not sample in prune_set: continue
				sample = util.cleanUpSampleName(sample)
				try:
					if prune_set != None and not sample in prune_set: continue
					assert (util.is_genbank(gbk))
					bgc_id = sample
					if sample_index[sample] > 0:
						bgc_id = sample + '_' + str(sample_index[sample] + 1)
					sample_index[sample] += 1

					is_expansion_bgc = False
					if '_Expansion_BGC' in gbk:
						is_expansion_bgc = True

					# Parse genbank using

					BGC_Object = BGC(gbk, bgc_id, is_expansion_bgc=is_expansion_bgc, prediction_method=prediction_method)
					BGC_Object.parseGenbanks(comprehensive_parsing=comprehensive_parsing)
					self.pan_bgcs[bgc_id] = BGC_Object
					self.comp_gene_info.update(BGC_Object.gene_information)
					self.bgc_info[bgc_id] = BGC_Object.cluster_information
					self.bgc_genes[bgc_id] = set(BGC_Object.gene_information)
					self.pan_genes = self.pan_genes.union(BGC_Object.gene_information.keys())
					self.bgc_gbk[bgc_id] = gbk
					self.bgc_sample[bgc_id] = sample
					self.sample_bgcs[sample].add(bgc_id)
					self.bgc_product[bgc_id] = [x['product'] for x in BGC_Object.cluster_information]
					self.bgc_core_counts[bgc_id] = BGC_Object.cluster_information[0]['count_core_gene_groups']

					if self.logObject:
						self.logObject.info("Incorporating genbank %s for sample %s into analysis." % (gbk, sample))
				except Exception as e:
					if self.logObject:
						self.logObject.warning("Unable to validate %s as Genbank. Skipping its incorporation into analysis." % gbk)
						self.logObject.warning(traceback.format_exc())
					raise RuntimeWarning(traceback.format_exc())

	def inputHomologyInformation(self, gene_to_hg, hg_genes, hg_median_copy_count, hg_prop_multi_copy):
		"""
		Simple function to store OrthoFinder homology information (parsed inititally in zol.utils)

		:param gene_to_hg: dictionary mapping gene locus tags to homolog group identifiers
		:param hg_genes: dictionary of sets with keys corresponding to homolog groups and values corresponding to set
						 of gene locus tags belonging to that homolog group.
		:param hg_median_copy_count: dictionary for the median copy count of each homolog group
		:param hg_prop_multi_copy: dictionary for the proportion of samples with homolog group which have multipmle
								   genes assigned to homolog group (paralogs).
		"""
		try:
			self.gene_to_hg = dict(gene_to_hg)
			self.hg_genes = dict(hg_genes)
			self.hg_median_copy_count = dict(hg_median_copy_count)
			self.hg_prop_multi_copy = dict(hg_prop_multi_copy)

			for bgc_id in self.bgc_genes:
				for gene in self.bgc_genes[bgc_id]:
					if gene in self.gene_to_hg:
						self.bgc_hgs[bgc_id].add(self.gene_to_hg[gene])

			if self.logObject:
				self.logObject.info(
					"Successfully inputted homolog information from OrthoFinder into variables of Pan object!")
		except Exception as e:
			if self.logObject:
				self.logObject.error("Had issues inputting homology information into Pan object instance.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def openStatsFile(self, outdir, run_parameter_tests=False):
		"""
		Simple function to initialize final report file with statistics on GCF clustering - main output from zol-Cluster.

		:param outdir: path to workspace directory
		:param run_parameter_tests: whether zol-Cluster is being run in "run_parameter_tests mode".
		"""
		try:
			final_stats_file = outdir + 'GCF_Details.txt'
			final_stats_expanded_singletons_file = outdir + 'GCF_Details_Expanded_Singletons.txt'
			sf_handle = open(final_stats_file, 'w')
			if run_parameter_tests:
				sf_handle.write('\t'.join(
					['MCL inflation parameter', 'Jaccard similarity cutoff', 'GCF id', 'number of BGCs',
					 'number of samples',
					 'samples with multiple BGCs in GCF', 'size of the SCC', 'mean number of OGs',
					 'stdev for number of OGs', 'number of core gene aggregates',
					 'min pairwise Jaccard similarity', 'max pairwise Jaccard similarity', 'annotations']) + '\n')
			else:
				sfes_handle = open(final_stats_expanded_singletons_file, 'w')
				sf_handle.write(
					'\t'.join(['GCF id', 'number of BGCs', 'number of samples', 'samples with multiple BGCs in GCF',
							   'size of the SCC', 'mean number of OGs', 'stdev for number of OGs',
							   'min pairwise Jaccard similarity', 'max pairwise Jaccard similarity', 'number of core gene aggregates', 'annotations']) + '\n')
				sfes_handle.write(
					'\t'.join(['GCF id', 'number of BGCs', 'number of samples', 'samples with multiple BGCs in GCF',
							   'size of the SCC', 'mean number of OGs', 'stdev for number of OGs',
							   'min pairwise Jaccard similarity', 'max pairwise Jaccard similarity', 'number of core gene aggregates', 'annotations']) + '\n')
				sfes_handle.close()

			sf_handle.close()
			self.final_stats_file = final_stats_file
			self.final_stats_expanded_singletons_file = final_stats_expanded_singletons_file
			if self.logObject:
				self.logObject.info(
					"Will be writing final stats report on GCF clustering to file %s" % final_stats_file)
		except Exception as e:
			if self.logObject:
				self.logObject.error("Had issues initializing final stats report file.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def calculateBGCPairwiseRelations(self, outdir, split_by_annotation=False):
		"""
		Function to calculate the Jaccard Similarity between pairs of BGCs based on homolog groups shared.

		:param outdir: path to workspace directory.
		:param split_by_annotation: partition BGCs foremost by general annotation category of BGC by AntiSMASH (similar
		                            to what BiG-Scape performs.
		"""
		try:
			pair_relations_txt_file = outdir + 'bgc_pair_relationships.txt'
			prf_handle = open(pair_relations_txt_file, 'w')

			pairwise_relations = defaultdict(lambda: defaultdict(float))
			pairwise_syntenic_relations = defaultdict(lambda: defaultdict(lambda: 0.0))

			bgc_hg_ranked_orders = defaultdict(dict)
			for bgc in self.bgc_hgs:
				sc_hg_starts = []
				for hg in self.bgc_hgs[bgc]:
					gene_counts = 0
					hg_start = 0
					for g in self.hg_genes[hg]:
						if self.comp_gene_info[g]['bgc_name'] == bgc:
							hg_start = self.comp_gene_info[g]['start']
							gene_counts += 1
					if gene_counts == 1:
						sc_hg_starts.append([hg, hg_start])

				for hg_ord, hg_its in enumerate(sorted(sc_hg_starts, key=itemgetter(1))):
					bgc_hg_ranked_orders[bgc][hg_its[0]] = (hg_ord + 1)

			for i, bgc1 in enumerate(self.bgc_hgs):
				bgc1_hgs = self.bgc_hgs[bgc1]
				bgc1_hgs_ranked_orders = bgc_hg_ranked_orders[bgc1]
				for j, bgc2 in enumerate(self.bgc_hgs):
					if i < j:
						bgc2_hgs = self.bgc_hgs[bgc2]
						bgc2_hgs_ranked_orders = bgc_hg_ranked_orders[bgc2]

						overlap_metric = float(len(bgc1_hgs.intersection(bgc2_hgs))) / float(len(bgc1_hgs.union(bgc2_hgs)))
						overlap_metric_scaled = 100.00 * overlap_metric
						pairwise_relations[bgc1][bgc2] = overlap_metric_scaled
						pairwise_relations[bgc2][bgc1] = overlap_metric_scaled

						sc_hgs_intersect = set(bgc1_hgs_ranked_orders.keys()).intersection(set(bgc2_hgs_ranked_orders.keys()))
						if len(sc_hgs_intersect) >= 3:
							bgc1_orders = []
							bgc2_orders = []
							for hg in sc_hgs_intersect:
								bgc1_orders.append(bgc1_hgs_ranked_orders[hg])
								bgc2_orders.append(bgc2_hgs_ranked_orders[hg])
							r, pval = stats.pearsonr(bgc1_orders, bgc2_orders)
							abs_correlation_coefficient = abs(r)
							pairwise_syntenic_relations[bgc1][bgc2] = abs_correlation_coefficient
							pairwise_syntenic_relations[bgc2][bgc1] = abs_correlation_coefficient

						products_union_not_empty = len(set(self.bgc_product[bgc1]).union(set(self.bgc_product[bgc2]))) > 0
						products_intersection_matches = float(len(set(self.bgc_product[bgc1]).intersection(set(self.bgc_product[bgc2]))))/float(len(set(self.bgc_product[bgc1]).union(set(self.bgc_product[bgc2])))) == 1.0
						if (not split_by_annotation) or (split_by_annotation and products_union_not_empty and products_intersection_matches):
							prf_handle.write('%s\t%s\t%f\n' % (bgc1, bgc2, overlap_metric_scaled))
			prf_handle.close()

			if self.logObject:
				self.logObject.info("Calculated pairwise relations and wrote to: %s" % pair_relations_txt_file)
			self.pairwise_relations = pairwise_relations
			self.pairwise_syntenic_relations = pairwise_syntenic_relations
			self.pair_relations_txt_file = pair_relations_txt_file
			self.bgc_to_gcf_map_file = outdir + 'BGC_to_GCF_Mapping.txt'
		except Exception as e:
			if self.logObject:
				self.logObject.error("Problem with creating relations file between pairs of BGCs.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def runMCLAndReportGCFs(self, mip, jcp, sccp, outdir, run_parameter_tests=False, cpus=1):
		"""
		Function to run MCL and report the GCFs (gene-cluster families) of homologous BGCs identified.

		:param mip: MCL inflation parameter.
		:param jcp: Jaccard similarity threshold for homology between two BGCs to be considered.
		:param outdir: path to workspace directory.
		:param run_parameter_tests: True
		:param cpus: number of cpus/threads to use for MCL.
		"""
		pair_relations_filt_txt_file = outdir + 'bgc_pair_relationships.%f.txt' % jcp
		try:
			prftf_handle = open(pair_relations_filt_txt_file, 'w')
			with open(self.pair_relations_txt_file) as oprtf:
				for line in oprtf:
					line = line.strip('\n')
					b1, b2, jaccard_sim = line.split('\t')
					if float(jaccard_sim) >= jcp and self.pairwise_syntenic_relations[b1][b2] >= sccp:
						prftf_handle.write(line + '\n')
			prftf_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error("Problem with parsing paired sample Jaccard similarity relationship file.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

		pair_relations_mci_file = outdir + 'bgc_pair_relationships.mci'
		pair_relations_tab_file = outdir + 'bgc_pair_relationships.tab'
		relations_mcl_file = outdir + 'mcl.' + str(mip).replace('.', '_') + '.out'
		mcxdump_out_file = outdir + 'final_mcl.' + str(mip).replace('.', '_') + '.out'

		mcxload_cmd = ['mcxload', '-abc', pair_relations_filt_txt_file, '--stream-mirror', '-write-tab',
					   pair_relations_tab_file,
					   '-o', pair_relations_mci_file]
		mcl_cmd = ['mcl', pair_relations_mci_file, '-I', str(mip), '-o', relations_mcl_file, '-te', str(cpus)]
		mcxdump_cmd = ['mcxdump', '-icl', relations_mcl_file, '-tabr', pair_relations_tab_file, '-o',
					   mcxdump_out_file]

		if self.logObject:
			self.logObject.info('Running the following command: %s' % ' '.join(mcxload_cmd))
		try:
			subprocess.call(' '.join(mcxload_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(mcxload_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(mcxload_cmd))
				self.logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue running: %s' % ' '.join(mcxload_cmd))

		if self.logObject:
			self.logObject.info('Converted format of pair relationship file via mxcload.')

		if self.logObject:
			self.logObject.info('Running MCL and MCXDUMP with inflation parameter set to %f' % mip)
			self.logObject.info('Running the following command: %s' % ' '.join(mcl_cmd))
		try:
			subprocess.call(' '.join(mcl_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(relations_mcl_file) and os.path.getsize(relations_mcl_file) > 50)
			self.logObject.info('Successfully ran: %s' % ' '.join(mcl_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(mcl_cmd))
				self.logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue running: %s' % ' '.join(mcl_cmd))

		if self.logObject:
			self.logObject.info('Running the following command: %s' % ' '.join(mcxdump_cmd))
		try:
			subprocess.call(' '.join(mcxdump_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(mcxdump_out_file) and os.path.getsize(mcxdump_out_file) > 50)
			if self.logObject:
				self.logObject.info('Successfully ran: %s' % ' '.join(mcl_cmd))
		except Exception as e:
			if self.logObject:
				self.logObject.error('Had an issue running: %s' % ' '.join(mcxdump_cmd))
				self.logObject.error(traceback.format_exc())
			raise RuntimeError('Had an issue running: %s' % ' '.join(mcxdump_cmd))

		if self.logObject:
			self.logObject.info('Successfully ran MCL and MCXDUMP with inflatimulti_same_sampleon parameter set to %f!' % mip)

		try:
			sf_handle = open(self.final_stats_file, 'a+')
			sfes_handle = open(self.final_stats_expanded_singletons_file, 'a+')

			clustered_bgcs = set([])
			gcf_identifier = 1
			with open(mcxdump_out_file) as omo:
				for gcf in omo:
					gcf = gcf.strip()
					gcf_mems = gcf.split()
					if len(gcf_mems) < 2: continue
					diffs = set([])
					samp_counts = defaultdict(int)
					samp_ogs = defaultdict(set)
					products = defaultdict(float)
					core_gene_cluster_counts = 0
					for a, bgc1 in enumerate(gcf_mems):
						samp1 = self.bgc_sample[bgc1]
						samp_counts[samp1] += 1
						for prod in self.bgc_product[bgc1]: products[prod] += 1.0 / len(self.bgc_product[bgc1])
						if core_gene_cluster_counts < self.bgc_core_counts[bgc1]: core_gene_cluster_counts = self.bgc_core_counts[bgc1]
						samp_ogs[samp1] = samp_ogs[samp1].union(self.bgc_hgs[bgc1])
						clustered_bgcs.add(bgc1)
						for b, bgc2 in enumerate(gcf_mems):
							if a < b:
								diffs.add(self.pairwise_relations[bgc1][bgc2])
					multi_same_sample = 0
					num_ogs = []
					for si, s in enumerate(samp_counts):
						if samp_counts[s] > 1:
							multi_same_sample += 1
						if si == 0:
							scc = samp_ogs[s]
						else:
							scc = scc.intersection(samp_ogs[s])
						num_ogs.append(len(samp_ogs[s]))
					stdev = "NA"
					mean = "NA"
					try:
						mean = statistics.mean(num_ogs)
						stdev = statistics.stdev(num_ogs)
					except:
						pass
					gcf_stats = ['GCF_' + str(gcf_identifier), len(gcf_mems), len(samp_counts.keys()), multi_same_sample,
								 len(scc),
								 mean, stdev, min(diffs),
								 max(diffs), core_gene_cluster_counts,
								 '; '.join([x[0] + ':' + str(x[1]) for x in products.items()])]
					if run_parameter_tests:
						gcf_stats = [mip, jcp] + gcf_stats
					else:
						sfes_handle.write('\t'.join([str(x) for x in gcf_stats]) + '\n')
					sf_handle.write('\t'.join([str(x) for x in gcf_stats]) + '\n')
					gcf_identifier += 1

			singleton_bgcs = []
			for bgc in self.bgc_hgs:
				if not bgc in clustered_bgcs:
					singleton_bgcs.append(bgc)
					if not run_parameter_tests:
						products = {}
						for prod in self.bgc_product[bgc]: products[prod] = 1.0 / len(self.bgc_product[bgc])
						mean_og_count = len(self.bgc_hgs[bgc])
						gcf_stats = ['GCF_' + str(gcf_identifier), 1, 1, 0, "NA", mean_og_count, 0, 0, 0, 0,  '; '.join([x[0] + ':' + str(x[1]) for x in products.items()]) ]
						sfes_handle.write('\t'.join([str(x) for x in gcf_stats]) + '\n')
						gcf_identifier += 1

			singleton_stats = ['singletons', len(singleton_bgcs)] + (['NA'] * 7)
			if run_parameter_tests:
				singleton_stats = [mip, jcp] + singleton_stats
			sf_handle.write('\t'.join([str(x) for x in singleton_stats]) + '\n')

			sf_handle.close()
			sfes_handle.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error(
					"Problem appending information on GCFs for current parameter combination to GCF statistics report file.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

		if not run_parameter_tests:
			try:
				if self.logObject:
					self.logObject.info(
					"Writing list of BGCs for each GCF, which will be used as input for downstream programs in the suite!")
				gcf_listing_dir = outdir + 'GCF_Listings/'
				if not os.path.isdir(gcf_listing_dir): os.system('mkdir %s' % gcf_listing_dir)
				gcf_identifier = 1
				with open(mcxdump_out_file) as omo:
					for gcf in omo:
						gcf = gcf.strip()
						gcf_mems = gcf.split()
						if len(gcf_mems) < 2: continue
						outf_list = open(gcf_listing_dir + 'GCF_' + str(gcf_identifier) + '.txt', 'w')
						for bgc in gcf_mems:
							outf_list.write('%s\t%s\n' % (self.bgc_sample[bgc], self.bgc_gbk[bgc]))
						outf_list.close()
						gcf_identifier += 1
				for bgc in singleton_bgcs:
					outf_list = open(gcf_listing_dir + 'GCF_' + str(gcf_identifier) + '.txt', 'w')
					outf_list.write('%s\t%s\n' % (self.bgc_sample[bgc], self.bgc_gbk[bgc]))
					gcf_identifier += 1

				if self.logObject:
					self.logObject.info("Successfully wrote lists of BGCs for each GCF.")
			except Exception as e:
				if self.logObject:
					self.logObject.error("Problem with writing GCF lists.")
					self.logObject.error(traceback.format_exc())
				raise RuntimeError(traceback.format_exc())
		else:
			try:
				if self.logObject:
					self.logObject.info(
					"Writing list of BGCs for each GCF for each clustering parameter combination into single plot.")
				btgmf_handle = open(self.bgc_to_gcf_map_file, 'a+')

				with open(mcxdump_out_file) as omo:
					for j, gcf in enumerate(omo):
						gcf = gcf.strip()
						gcf_mems = gcf.split()
						if len(gcf_mems) < 2: continue
						for bgc in gcf_mems:
							sname = self.bgc_sample[bgc]
							btgmf_handle.write('%f\t%f\tGCF_%d\t%s\t%s\n' % (mip, jcp, j + 1, sname, bgc))
				for sbgc in singleton_bgcs:
					sname = self.bgc_sample[sbgc]
					btgmf_handle.write('%f\t%f\tGCF_singletons\t%s\t%s\n' % (mip, jcp, sname, sbgc))
				btgmf_handle.close()
			except Exception as e:
				if self.logObject:
					self.logObject.error("Problem with writing BGC to GCF relationship.")
					self.logObject.error(traceback.format_exc())
				raise RuntimeError(traceback.format_exc())

	def plotResultsFromUsingDifferentParameters(self, outdir):
		"""
		Function to create an 8x11 in PDF with figures aimed to give the user an idea on how different values of the
		MCL inflation parameter and the Jaccard similarity cutoff might impact clustering to help them choose the best
		parameter values.

		:param outdir: path to workspace directory.
		"""
		try:
			plot_input_dir = outdir + 'plotting_input/'
			if not os.path.isdir(plot_input_dir): os.system('mkdir %s' % plot_input_dir)

			singleton_counts = defaultdict(int)
			clustered_counts = defaultdict(int)
			cluster_counts = defaultdict(int)
			cluster_mixedannot_counts = defaultdict(int)
			cluster_mixedannot_wohypo_counts = defaultdict(int)
			all_annotation_classes = set([])
			with open(self.final_stats_file) as ogdf:
				for i, line in enumerate(ogdf):
					line = line.strip()
					ls = line.split("\t")
					if i == 0: continue
					mcl = round(float(ls[0]), 2)
					jcp = round(float(ls[1]), 2)
					param_id = '%f_%f' % (mcl, jcp)
					gcf_id = ls[2]
					if gcf_id == 'singletons':
						singleton_counts[param_id] = int(ls[3])
					else:
						clustered_counts[param_id] += int(ls[3])
						cluster_counts[param_id] += 1
						annotations = set(
							[x.split(':')[0].replace('-', '.') for x in ls[-1].split('; ') if
							 x.split(':')[0] != 'NA'])
						if len(ls[-1].split('; ')) > 1: cluster_mixedannot_counts[param_id] += 1
						if len(annotations) > 1: cluster_mixedannot_wohypo_counts[param_id] += 1
						all_annotation_classes = all_annotation_classes.union(annotations)

			plot_overview_file = plot_input_dir + 'plotting_input_1.txt'
			plot_input_1_file = plot_input_dir + 'plotting_input_2.txt'
			plot_input_2_file = plot_input_dir + 'plotting_input_3.txt'
			plot_sankey_file = plot_input_dir + 'plotting_input_4.txt'
			plot_annot_file = plot_input_dir + 'plotting_input_5.txt'

			plot_input_1_handle = open(plot_input_1_file, 'w')
			plot_input_2_handle = open(plot_input_2_file, 'w')
			plot_annot_handle = open(plot_annot_file, 'w')
			plot_overview_handle = open(plot_overview_file, 'w')
			plot_sankey_handle = open(plot_sankey_file, 'w')

			plot_annot_handle.write('\n'.join(['Unknown'] + sorted(list(all_annotation_classes))))
			plot_annot_handle.close()

			coord_tuples = defaultdict(set)
			plot_input_1_handle.write('\t'.join(
				['GCF', 'Parameters', 'JaccardSim', 'Inflation', 'Samples', 'Samples_Offset',
				 'MultiBGCSamples_Offset',
				 'SCCExists', 'SCCSize', 'CoreGeneClusters', 'AvgGeneCount', 'StdDevGeneCount', 'Unknown'] + sorted(
					list(all_annotation_classes))) + '\n')
			plot_input_2_handle.write('\t'.join(['GCF', 'Parameters', 'Annotation', 'Count', 'Total BGCs']) + '\n')
			plot_overview_handle.write('\t'.join(
				['Parameters', 'JaccardSim', 'Inflation', 'Clustered', 'Singletons', 'SingletonToClustered',
				 'MixedAnnotations', 'MixedAnnotationsWoHypothetical', 'NumberClusters',
				 'MixedAnnotationProportion',
				 'MixedAnnotationWoHypotheticalProportion']) + '\n')
			plot_sankey_handle.write('\t'.join(
				['JaccardSim', 'Parameter_Source', 'Parameter_Dest', 'Source_Size', 'Dest_Size', 'Source_Inflation',
				 'Source_Type', 'Dest_Type', 'PassageWeight']) + '\n')

			jcp_gcfs = defaultdict(lambda: defaultdict(set))
			with open(self.bgc_to_gcf_map_file) as obgmf:
				for line in obgmf:
					line = line.strip()
					mip, jcp, gcf_id, sname, gbkpath = line.split('\t')
					jcp = float(jcp)
					mip = float(mip)
					param_mod_id = 'Inf: %.1f JacSim: %.1f' % (mip, jcp)
					gcf_name = gcf_id + ' - ' + param_mod_id
					jcp_gcfs[jcp][gcf_name].add(gbkpath)

			inflation_order = {0.8: 1, 1.4: 2, 2: 3, 2.5: 4, 3: 5, 3.5: 6, 4: 7, 5: 8}
			for jcp in jcp_gcfs:
				for i, gcf1 in enumerate(jcp_gcfs[jcp]):
					for j, gcf2 in enumerate(jcp_gcfs[jcp]):
						gcf1_inf = float(gcf1.split('Inf: ')[1].split(' JacSim:')[0].strip())
						gcf2_inf = float(gcf2.split('Inf: ')[1].split(' JacSim:')[0].strip())
						if (inflation_order[gcf1_inf] + 1) != (inflation_order[gcf2_inf]): continue
						gcf1_bgcs = jcp_gcfs[jcp][gcf1]
						gcf2_bgcs = jcp_gcfs[jcp][gcf2]
						intersect = len(gcf1_bgcs.intersection(gcf2_bgcs))
						if intersect > 0:
							gcf1_type = 'cluster';
							gcf2_type = 'cluster'
							if 'singleton' in gcf1: gcf1_type = 'singleton'
							if 'singleton' in gcf2: gcf2_type = 'singleton'
							plot_sankey_handle.write('\t'.join([str(x) for x in
																[jcp, gcf1, gcf2, len(gcf1_bgcs), len(gcf2_bgcs),
																 '%.1f vs. %.1f' % (gcf1_inf, gcf2_inf), gcf1_type,
																 gcf2_type, intersect]]) + '\n')
			plot_sankey_handle.close()

			with open(self.final_stats_file) as ogdf:
				for i, line in enumerate(ogdf):
					line = line.strip()
					ls = line.split("\t")
					if i == 0: continue
					mcl = round(float(ls[0]), 2)
					jcp = round(float(ls[1]), 2)
					param_id = '%f_%f' % (mcl, jcp)
					gcf_id = ls[2]
					if gcf_id == 'singletons': continue

					param_mod_id = 'MCL Inflation: %.1f; Jaccard Similarity Cutoff: %.1f' % (mcl, jcp)
					# get stats for plotting

					samples = int(ls[4])
					samples_with_multi_bgcs = int(ls[5])
					size_of_scc = int(ls[6])
					number_of_core_genes = ls[-2]
					avg_gene_count = ls[-6]
					stdev_gene_count = ls[-5]
					annotations = ls[-1].split('; ')

					unique_coord = False
					random_offset_1 = random.randint(-50, 51) / 100.0
					random_offset_2 = random.randint(-50, 51) / 100.0
					coord_tuple = tuple([samples + random_offset_1, samples_with_multi_bgcs + random_offset_2])
					while not unique_coord:
						if not coord_tuple in coord_tuples[param_mod_id]:
							unique_coord = True
						else:
							random_offset_1 = random.randint(-50, 51) / 100.0
							random_offset_2 = random.randint(-50, 51) / 100.0
							coord_tuple = tuple(
								[samples + random_offset_1, samples_with_multi_bgcs + random_offset_2])
					coord_tuples[param_mod_id].add(coord_tuple)
					scc_exists = "SCC DNE"
					if size_of_scc > 0:
						scc_exists = "SCC Exists"

					annot_count = defaultdict(float)
					for an in annotations:
						ans = an.split(':')
						annot_count[ans[0].replace('-', '.')] = float(ans[1])

					annotation_class_abd = [annot_count['NA']]
					for ac in sorted(all_annotation_classes):
						annotation_class_abd.append(annot_count[ac])

					plot_input_1_handle.write('\t'.join([str(x) for x in
														 [gcf_id + ' - ' + param_mod_id, param_mod_id, jcp, mcl,
														  samples, samples + random_offset_1,
														  samples_with_multi_bgcs + random_offset_2, scc_exists,
														  size_of_scc, number_of_core_genes, avg_gene_count,
														  stdev_gene_count] + annotation_class_abd]) + '\n')
					if int(gcf_id.split('_')[1]) <= 50:
						tot_count = sum(annot_count.values())
						for an in annot_count:
							an_count = annot_count[an]
							if an == 'NA':
								an = 'Unknown'
							plot_input_2_handle.write('\t'.join([str(x) for x in
																 [gcf_id + ' - ' + param_mod_id, param_mod_id, an,
																  an_count, tot_count]]) + '\n')

			for pid in clustered_counts:
				mcl, jcp = [float(x) for x in pid.split('_')]
				param_mod_id = 'MCL Inflation: %.1f; Jaccard Similarity Cutoff: %.1f' % (mcl, jcp)

				plot_overview_handle.write(
					'\t'.join([str(x) for x in [param_mod_id, jcp, mcl, clustered_counts[pid],
												singleton_counts[pid],
												singleton_counts[pid] / clustered_counts[pid],
												cluster_mixedannot_counts[pid],
												cluster_mixedannot_wohypo_counts[pid],
												cluster_counts[pid],
												cluster_mixedannot_counts[pid] / float(
													cluster_counts[pid]),
												cluster_mixedannot_wohypo_counts[pid] / float(
													cluster_counts[pid])]]) + '\n')

				for sgcf_id in range(1, singleton_counts[pid] + 1):
					sgcf_name = 'SGCF_' + str(sgcf_id)
					plot_input_1_handle.write('\t'.join(
						[sgcf_name + ' - ' + param_mod_id, param_mod_id, str(jcp), str(mcl), '1', 'NA', 'NA',
						 'Not Applicable'] + ['NA'] * (5 + len(all_annotation_classes))) + '\n')

			plot_input_1_handle.close()
			plot_input_2_handle.close()
			plot_overview_handle.close()

			plot_pdf_file = outdir + 'Plots_Depicting_Parameter_Influence_on_GCF_Clustering.pdf'
			rscript_plot_cmd = ["Rscript", RSCRIPT_FOR_CLUSTER_ASSESSMENT_PLOTTING, plot_input_1_file,
								plot_annot_file,
								plot_input_2_file, plot_overview_file, plot_sankey_file, plot_pdf_file]
			if self.logObject:
				self.logObject.info('Running R-based plotting with the following command: %s' % ' '.join(rscript_plot_cmd))
			try:
				subprocess.call(' '.join(rscript_plot_cmd), shell=True, stdout=sys.stderr,
								stderr=sys.stderr,
								executable='/bin/bash')
				if self.logObject:
					self.logObject.info('Successfully ran: %s' % ' '.join(rscript_plot_cmd))
			except Exception as e:
				if self.logObject:
					self.logObject.error('Had an issue running: %s' % ' '.join(rscript_plot_cmd))
					self.logObject.error(traceback.format_exc())
				raise RuntimeError('Had an issue running: %s' % ' '.join(rscript_plot_cmd))

			if self.logObject:
				self.logObject.info('Plotting completed!')

		except Exception as e:
			if self.logObject:
				self.logObject.error(
					"Problem creating plot(s) for assessing best parameter choices for GCF clustering.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def convertGenbanksIntoFastas(self, fasta_dir, fasta_listing_file, sample_retention_set=None):
		"""
		Function to convert Genbanks for BGC instances into FASTA format and listing files. Note, there will
		be one FASTA per sample not per BGC Genbank (in case sample's have more than one associated Genbank).

		:param listing_file: tab-delimited file with two columns: (1) sample name (2) path to BGC Genbank
		:param fasta_listing_file: tab-delimited file with two columns: (1) sample name (2) path to BGC FASTA
		"""
		try:
			sample_index = defaultdict(int)

			all_samples = set([])
			with open(self.bgc_genbanks_listing) as obsf:
				for i, line in enumerate(obsf):
					line = line.strip()
					try:
						assert (len(line.split('\t')) >= 2)
					except Exception as e:
						if self.logObject:
							self.logObject.error(
							"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (
									i + 1))
							self.logObject.error(traceback.format_exc())
						raise RuntimeError(
							"More than two columns exist at line %d in BGC specification/listing file. Exiting now ..." % (
									i + 1))
					if len(line.split('\t')) == 2:
						sample, gbk = line.split('\t')
					else:
						sample, gbk = line.split('\t')[:-1]
					sample = util.cleanUpSampleName(sample)
					if sample_retention_set != None and sample not in sample_retention_set: continue
					try:
						assert (util.is_genbank(gbk))
						bgc_id = sample
						if sample_index[sample] > 0:
							bgc_id = sample + '_' + str(sample_index[sample] + 1)
						sample_index[sample] += 1

						sample_fasta = fasta_dir + sample + '.fasta'
						sample_handle = open(sample_fasta, 'a+')
						with open(gbk) as ogbk:
							for ri, rec in enumerate(SeqIO.parse(ogbk, 'genbank')):
								sample_handle.write('>' + bgc_id + '|' + str(ri) + '\n' + str(rec.seq) + '\n')
						sample_handle.close()
						all_samples.add(sample)
						if self.logObject:
							self.logObject.info("Added Genbank %s into sample %s's GCF relevant FASTA." % (gbk, sample))
					except Exception as e:
						if self.logObject:
							self.logObject.warning("Unable to validate %s as Genbank. Skipping its incorporation into analysis.")
							self.logObject.warning(traceback.format_exc())
						raise RuntimeWarning("Unable to validate %s as Genbank. Skipping ...")

			outf = open(fasta_listing_file, 'w')
			for s in all_samples:
				outf.write(s + '\t' + fasta_dir + s + '.fasta\n')
			outf.close()
		except Exception as e:
			if self.logObject:
				self.logObject.error("Issues converting BGC Genbank listing into FASTA listing.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def readInPopulationsSpecification(self, pop_specs_file, prune_set=None):
		"""
		Read in population specifications per sample and assign population to each BGC.

		:param pop_specs_file: path to file which has tab separated sample (1st column) and corresponding population (2nd column)
		"""
		try:
			self.bgc_population = defaultdict(lambda: "NA")
			self.sample_population = defaultdict(lambda: "NA")
			with open(pop_specs_file) as opsf:
				for i, line in enumerate(opsf):
					if i == 0 and line.startswith('name'): continue
					line = line.strip()
					sample, population = line.split('\t')
					if prune_set != None and not sample in prune_set: continue
					self.sample_population[sample] = population
					for bgc in self.sample_bgcs[sample]:
						self.bgc_population[bgc] = population
			self.logObject.info("Successfully parsed population specifications file. There are %d populations." % len(population))
		except Exception as e:
			if self.logObject:
				self.logObject.error("Issue in parsing population specifications per sample and associating with each BGC.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())

	def constructHMMProfiles(self, outdir, initial_sample_prokka_data, cpus=1, quick_mode=False):
		"""
		Wrapper function to construct Hmmer3 HMMs for each of the homolog groups.

		:param outdir: The path to the workspace / output directory.
		:param initial_sample_prokka_data: Dictionary with keys being sample identifiers and values being paths to
										   genbank and proteome files from Prokka based annotation
		:param cpus: The number of cpus (will be used for parallelizing)
		"""
		outdir = os.path.abspath(outdir) + '/'
		prot_seq_dir = outdir + 'Protein_Sequences/'
		prot_alg_dir = outdir + 'Protein_Alignments/'
		prot_hmm_dir = outdir + 'Profile_HMMs/'
		search_ref_res_dir = outdir + 'Reflexive_Alignment_Results/'
		hg_differentiation_file = open(outdir + 'Homolog_Groups_Differentiation.txt', 'w')

		if not os.path.isdir(prot_seq_dir): os.system('mkdir %s' % prot_seq_dir)
		if not os.path.isdir(prot_alg_dir): os.system('mkdir %s' % prot_alg_dir)
		if not os.path.isdir(prot_hmm_dir): os.system('mkdir %s' % prot_hmm_dir)
		if not os.path.isdir(search_ref_res_dir): os.system('mkdir %s' % search_ref_res_dir)

		try:
			inputs = []
			sample_hgs = defaultdict(set)
			for hg in self.hg_genes:
				sample_sequences = {}
				for gene in self.hg_genes[hg]:
					gene_info = self.comp_gene_info[gene]
					bgc_id = gene_info['bgc_name']
					sample_id = self.bgc_sample[bgc_id]
					sample_hgs[sample_id].add(hg)
					prot_seq = gene_info['prot_seq']
					sample_sequences[sample_id] = prot_seq
				inputs.append([hg, sample_sequences, prot_seq_dir, prot_alg_dir, prot_hmm_dir, self.logObject])

			sample_hg_counts = [len(sample_hgs[x]) for x in sample_hgs]
			self.lowerbound_hg_count = math.floor(min(sample_hg_counts))

			p = multiprocessing.Pool(cpus)
			p.map(create_hmm_profiles, inputs)
			p.close()

			if self.logObject:
				self.logObject.info("Successfully created profile HMMs for each homolog group. Now beginning concatenation into single file.")

			self.concatenated_profile_HMM = outdir + 'All_GCF_Homologs.hmm'
			os.system('rm -f %s %s.h3*' % (self.concatenated_profile_HMM, self.concatenated_profile_HMM))

			for f in os.listdir(prot_hmm_dir):
				os.system('cat %s >> %s' % (prot_hmm_dir + f, self.concatenated_profile_HMM))

			if quick_mode:
				self.consensus_sequence_HMM = outdir + 'All_GCF_Homologs_ConsensusSequences.fasta'
				self.consensus_sequence_HMM_db = outdir + 'All_GCF_Homologs_ConsensusSequences'

				hmmemit_cmd = ['hmmemit', '-o', self.consensus_sequence_HMM, '-c', self.concatenated_profile_HMM]
				if self.logObject:
					self.logObject.info(
						'Running hmmemit on concatenated profiles HMM to get consensus sequences with the following command: %s' % ' '.join(hmmemit_cmd))
				try:
					subprocess.call(' '.join(hmmemit_cmd), shell=True, stdout=subprocess.DEVNULL,
									stderr=subprocess.DEVNULL,
									executable='/bin/bash')
					if self.logObject:
						self.logObject.info('Successfully ran: %s' % ' '.join(hmmemit_cmd))
				except:
					if self.logObject:
						self.logObject.error('Had an issue running: %s' % ' '.join(hmmemit_cmd))
						self.logObject.error(traceback.format_exc())
					raise RuntimeError('Had an issue running: %s' % ' '.join(hmmemit_cmd))

				makedb_cmd = ['diamond', 'makedb', '--in', self.consensus_sequence_HMM, '-d', self.consensus_sequence_HMM_db]
				if self.logObject:
					self.logObject.info(
						'Running Diamond makedb on concatenated profiles with the following command: %s' % ' '.join(makedb_cmd))
				try:
					subprocess.call(' '.join(makedb_cmd), shell=True, stdout=subprocess.DEVNULL,
									stderr=subprocess.DEVNULL,
									executable='/bin/bash')
					if self.logObject:
						self.logObject.info('Successfully ran: %s' % ' '.join(makedb_cmd))
				except:
					if self.logObject:
						self.logObject.error('Had an issue running: %s' % ' '.join(makedb_cmd))
						self.logObject.error(traceback.format_exc())
					raise RuntimeError('Had an issue running: %s' % ' '.join(makedb_cmd))

				diamond_cmds = []
				for sample in initial_sample_prokka_data:
					sample_proteome = initial_sample_prokka_data[sample]['predicted_proteome']
					result_file = search_ref_res_dir + sample + '.txt'
					diamond_cmd = ['diamond', 'blastp', '--threads', '1',  '--sensitive', '--db',
								   self.consensus_sequence_HMM_db, '--query', sample_proteome, '--outfmt', '6', '--out',
								   result_file, self.logObject]
					diamond_cmds.append(diamond_cmd)

				p = multiprocessing.Pool(cpus)
				p.map(util.multiProcess, diamond_cmds)
				p.close()

			else:
				hmmpress_cmd = ['hmmpress', self.concatenated_profile_HMM]
				if self.logObject:
					self.logObject.info(
						'Running hmmpress on concatenated profiles with the following command: %s' % ' '.join(hmmpress_cmd))
				try:
					subprocess.call(' '.join(hmmpress_cmd), shell=True, stdout=subprocess.DEVNULL,
									stderr=subprocess.DEVNULL,
									executable='/bin/bash')
					if self.logObject:
						self.logObject.info('Successfully ran: %s' % ' '.join(hmmpress_cmd))
				except:
					if self.logObject:
						self.logObject.error('Had an issue running: %s' % ' '.join(hmmpress_cmd))
						self.logObject.error(traceback.format_exc())
					raise RuntimeError('Had an issue running: %s' % ' '.join(hmmpress_cmd))

				hmmscan_cmds = []
				for sample in initial_sample_prokka_data:
					sample_proteome = initial_sample_prokka_data[sample]['predicted_proteome']
					result_file = search_ref_res_dir + sample + '.txt'
					hmmscan_cmd = ['hmmscan', '--max', '--cpu', '1', '--tblout', result_file, self.concatenated_profile_HMM,
								   sample_proteome, self.logObject]
					hmmscan_cmds.append(hmmscan_cmd)

				p = multiprocessing.Pool(cpus)
				p.map(util.multiProcess, hmmscan_cmds)
				p.close()

			best_hits = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
			for sample in initial_sample_prokka_data:
				result_file = search_ref_res_dir + sample + '.txt'
				assert (os.path.isfile(result_file))

				hits_per_gene = defaultdict(list)
				if quick_mode:
					with open(result_file) as orf:
						for line in orf:
							if line.startswith("#"): continue
							line = line.strip()
							ls = line.split()
							hg = ls[1].split('-consensus')[0]
							gene_id = ls[0]
							eval = decimal.Decimal(ls[10])
							hits_per_gene[gene_id].append([hg, eval])
				else:
					with open(result_file) as orf:
						for line in orf:
							if line.startswith("#"): continue
							line = line.strip()
							ls = line.split()
							hg = ls[0]
							gene_id = ls[2]
							eval = decimal.Decimal(ls[4])
							hits_per_gene[gene_id].append([hg, eval])

				for gene_id in hits_per_gene:
					for i, hit in enumerate(sorted(hits_per_gene[gene_id], key=itemgetter(1))):
						if i == 0:
							in_hg_flag = False
							if gene_id in self.gene_to_hg and hit[0] == self.gene_to_hg[gene_id]: in_hg_flag = True
							best_hits[hit[0]][sample][in_hg_flag].append(hit[1])

			for hg in best_hits:
				true_hits = []
				false_hits = []
				for sample in best_hits[hg]:
					true_hits += best_hits[hg][sample][True]
					false_hits += best_hits[hg][sample][False]
				worst_true_hit = 10000.0
				best_false_hit = 10000.0
				if len(true_hits) > 0:
					worst_true_hit = max(true_hits)
				if len(false_hits) > 0:
					best_false_hit = min(false_hits)
				able_to_differentiate = False
				eval_threshold = decimal.Decimal(1e-10)
				if best_false_hit > worst_true_hit:
					able_to_differentiate = True
					if worst_true_hit == 0.0:
						worst_true_hit = 1e-300
						if quick_mode:
							worst_true_hit = 1e-500
					if best_false_hit == 0.0:
						best_false_hit = 1e-300
						if quick_mode:
							best_false_hit = 1e-500
					wth_log10 = decimal.Decimal(worst_true_hit).log10()
					bfh_log10 = decimal.Decimal(best_false_hit).log10()
					eval_threshold = max([(decimal.Decimal(10.0) ** decimal.Decimal((wth_log10 + bfh_log10)/decimal.Decimal(2.0))), (decimal.Decimal(10.0) ** decimal.Decimal(bfh_log10 + decimal.Decimal(-5)))])
					eval_threshold = min([decimal.Decimal(eval_threshold), decimal.Decimal(1e-10)])
				hg_differentiation_file.write('\t'.join([hg, str(worst_true_hit), str(best_false_hit), str(able_to_differentiate), str(eval_threshold)]) + '\n')
				self.hg_differentiation_stats[hg] = {'worst_true_hit': worst_true_hit, 'best_false_hit': best_false_hit,
													 'able_to_differentiate': able_to_differentiate}
				self.hg_max_self_evalue[hg] = [eval_threshold, able_to_differentiate]
		except:
			if self.logObject:
				self.logObject.error("Issues with constructing profile HMMs.")
				self.logObject.error(traceback.format_exc())
			raise RuntimeError(traceback.format_exc())
		hg_differentiation_file.close()

	def runHMMScan(self, outdir, expanded_sample_prokka_data, cpus=1, quick_mode=False, annotation_pickle_file=None):
		"""
		Function to run hmmscan and process results as well as read in Genbanks for expanded list of samples.

		:param outdir: The path to the workspace / output directory.
		:param expanded_sample_prokka_data: Dictionary with keys being sample identifiers and values being paths to
											genbank and proteome files from Prokka based annotation
		:param cpus: The number of cpus (will be used for parallelizing)
		"""
		search_res_dir = os.path.abspath(outdir + 'Alignment_Results') + '/'
		if not os.path.isdir(search_res_dir): os.system('mkdir %s' % search_res_dir)

		if not annotation_pickle_file:
			with multiprocessing.Manager() as manager:
				sample_gbk_info = manager.dict()
				genbanks = []
				for sample in expanded_sample_prokka_data:
					sample_genbank = expanded_sample_prokka_data[sample]['genbank']
					genbanks.append([sample, sample_genbank, sample_gbk_info])

				with manager.Pool(cpus) as pool:
					pool.map(util.parseGenbankAndFindBoundaryGenes, genbanks)

				for sample in sample_gbk_info:
					gene_location, scaff_genes, bound_genes, gito, goti = sample_gbk_info[sample]
					self.gene_location[sample] = gene_location
					self.scaffold_genes[sample] = scaff_genes
					self.boundary_genes[sample] = bound_genes
					self.gene_id_to_order[sample] = gito
					self.gene_order_to_id[sample] = goti
		else:
			try:
				with open(annotation_pickle_file, "rb") as input_file:
					genbank_info = cPickle.load(input_file)
					for sample in genbank_info:
						self.gene_location[sample] = genbank_info[sample][0]
						self.scaffold_genes[sample] = genbank_info[sample][1]
						self.boundary_genes[sample] = genbank_info[sample][2]
						self.gene_id_to_order[sample] = genbank_info[sample][3]
						self.gene_order_to_id[sample] = genbank_info[sample][4]
			except Exception as e:
				if self.logObject:
					self.logObject.error("Had issues loading data from extention sample set pickle object.")
					self.logObject.error(traceback.format_exc())
				raise RuntimeError(traceback.format_exc())

		alignment_cmds = []
		for sample in expanded_sample_prokka_data:
			sample_proteome = expanded_sample_prokka_data[sample]['predicted_proteome']
			result_file = search_res_dir + sample + '.txt'
			if quick_mode:
				diamond_cmd = ['diamond', 'blastp', '--threads', '1', '--sensitive', '--db',
							   self.consensus_sequence_HMM_db, '--query', sample_proteome, '--outfmt', '6', '--out',
							   result_file, self.logObject]
				alignment_cmds.append(diamond_cmd)
			else:
				hmmscan_cmd = ['hmmscan', '--max', '--cpu', '1', '--tblout', result_file, self.concatenated_profile_HMM,
							   sample_proteome, self.logObject]
				alignment_cmds.append(hmmscan_cmd)

		p = multiprocessing.Pool(cpus)
		p.map(util.multiProcess, alignment_cmds)
		p.close()

		hg_valid_length_range = {}
		for hg in self.hg_genes:
			gene_lengths = []
			for g in self.hg_genes[hg]:
				gstart = self.comp_gene_info[g]['start']
				gend = self.comp_gene_info[g]['end']
				gene_lengths.append(abs(gstart - gend))
			min_gene_nucl_seq_len = min(gene_lengths)
			max_gene_nucl_seq_len = max(gene_lengths)
			median_gene_nucl_seq_lens = statistics.median(gene_lengths)
			mad_gene_nucl_seq_lens = max(stats.median_abs_deviation(gene_lengths, scale="normal"), 25)
			hg_valid_length_range[hg] = {'min_gene_length': min_gene_nucl_seq_len,
										 'max_gene_length': max_gene_nucl_seq_len,
										 'median_gene_length': median_gene_nucl_seq_lens,
										 'gene_length_deviation': mad_gene_nucl_seq_lens}

		for sample in expanded_sample_prokka_data:
			result_file = search_res_dir + sample + '.txt'
			assert (os.path.isfile(result_file))

			best_hit_per_gene = defaultdict(lambda: [set([]), 100000.0])
			with open(result_file) as orf:
				for line in orf:
					line = line.strip()
					ls = line.split()
					if quick_mode:
						hg = ls[1].split('-consensus')[0]
						gene_id = ls[0]
						eval = decimal.Decimal(ls[10])
						if eval < best_hit_per_gene[gene_id][1]:
							best_hit_per_gene[gene_id][1] = eval
							best_hit_per_gene[gene_id][0] = set([hg])
						elif eval == best_hit_per_gene[gene_id][1]:
							best_hit_per_gene[gene_id][0].add(hg)
					else:
						if line.startswith("#"): continue
						hg = ls[0]
						gene_id = ls[2]
						eval = decimal.Decimal(ls[4])
						if eval < best_hit_per_gene[gene_id][1]:
							best_hit_per_gene[gene_id][1] = eval
							best_hit_per_gene[gene_id][0] = set([hg])
						elif eval == best_hit_per_gene[gene_id][1]:
							best_hit_per_gene[gene_id][0].add(hg)

			with open(result_file) as orf:
				for line in orf:
					line = line.strip()
					ls = line.split()
					if quick_mode:
						hg = ls[1].split('-consensus')[0]
						gene_id = ls[0]
						gene_length = abs(self.gene_location[sample][gene_id]['start'] - self.gene_location[sample][gene_id]['end'])
						is_boundary_gene = gene_id in self.boundary_genes[sample]
						if not hg in best_hit_per_gene[gene_id][0]: continue
						scaffold = self.gene_location[sample][gene_id]['scaffold']
						eval = decimal.Decimal(ls[10])
						if eval < decimal.Decimal(1e-10):
							self.hmmscan_results_lenient[sample][gene_id] = [hg, eval]
						if (not is_boundary_gene) and \
								(not (gene_length <= hg_valid_length_range[hg]['max_gene_length'] and gene_length >= hg_valid_length_range[hg]['min_gene_length'])) and \
								(abs(gene_length - hg_valid_length_range[hg]['median_gene_length']) >= (2 * hg_valid_length_range[hg]['gene_length_deviation'])): continue
						if eval <= self.hg_max_self_evalue[hg][0]:
							self.hmmscan_results[gene_id].append([hg, eval, sample, scaffold])
						elif eval <= decimal.Decimal(1e-10) and is_boundary_gene:
							self.hmmscan_results[gene_id].append([hg, eval, sample, scaffold])
					else:
						if line.startswith("#"): continue
						hg = ls[0]
						gene_id = ls[2]
						gene_length = abs(self.gene_location[sample][gene_id]['start'] - self.gene_location[sample][gene_id]['end'])
						is_boundary_gene = gene_id in self.boundary_genes[sample]
						if not hg in best_hit_per_gene[gene_id][0]: continue
						scaffold = self.gene_location[sample][gene_id]['scaffold']
						eval = decimal.Decimal(ls[4])
						if eval < decimal.Decimal(1e-10):
							self.hmmscan_results_lenient[sample][gene_id] = [hg, eval]
						if (not is_boundary_gene) and \
								(not (gene_length <= hg_valid_length_range[hg]['max_gene_length'] and gene_length >= hg_valid_length_range[hg]['min_gene_length'])) and \
								(abs(gene_length - hg_valid_length_range[hg]['median_gene_length']) >= (2 * hg_valid_length_range[hg]['gene_length_deviation'])): continue
						if eval <= self.hg_max_self_evalue[hg][0]:
							self.hmmscan_results[gene_id].append([hg, eval, sample, scaffold])
						elif eval <= decimal.Decimal(1e-10) and is_boundary_gene:
							self.hmmscan_results[gene_id].append([hg, eval, sample, scaffold])

def create_hmm_profiles(inputs):
	"""
	Function to create MAFFT based MSAs of proteins for homolog group and then use Hmmer3 to create HMM profiles.
	"""
	hg, sample_sequences, prot_seq_dir, prot_alg_dir, prot_hmm_dir, logObject = inputs

	hg_prot_fasta = prot_seq_dir + '/' + hg + '.faa'
	hg_prot_msa = prot_alg_dir + '/' + hg + '.msa.faa'
	hg_prot_hmm = prot_hmm_dir + '/' + hg + '.hmm'

	hg_prot_handle = open(hg_prot_fasta, 'w')
	for s in sample_sequences:
		hg_prot_handle.write('>' + s + '\n' + str(sample_sequences[s]) + '\n')
	hg_prot_handle.close()

	mafft_cmd = ['mafft', '--maxiterate', '1000', '--localpair', hg_prot_fasta, '>', hg_prot_msa]
	if logObject:
		logObject.info('Running mafft with the following command: %s' % ' '.join(mafft_cmd))
	try:
		subprocess.call(' '.join(mafft_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		if logObject:
			logObject.info('Successfully ran: %s' % ' '.join(mafft_cmd))
	except:
		if logObject:
			logObject.error('Had an issue running: %s' % ' '.join(mafft_cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(mafft_cmd))

	hmmbuild_cmd = ['hmmbuild', '--amino', '-n', hg, hg_prot_hmm, hg_prot_msa]
	if logObject:
		logObject.info('Running hmmbuild (from HMMER3) with the following command: %s' % ' '.join(hmmbuild_cmd))
	try:
		subprocess.call(' '.join(hmmbuild_cmd), shell=True, stdout=subprocess.DEVNULL,
						stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		if logObject:
			logObject.info('Successfully ran: %s' % ' '.join(hmmbuild_cmd))
	except:
		if logObject:
			logObject.error('Had an issue running: %s' % ' '.join(hmmbuild_cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(hmmbuild_cmd))

	if logObject:
		logObject.info('Constructed profile HMM for homolog group %s' % hg)
