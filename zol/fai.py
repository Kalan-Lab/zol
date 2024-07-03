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
from pomegranate.hmm import DenseHMM
from pomegranate.distributions import Categorical
import numpy
import traceback
from scipy.stats import pearsonr
import shutil
from ete3 import Tree
import pickle
import statistics
import pandas as pd
import torch

# code for setup and finding location of programs based on conda vs. bioconda installation
zol_exec_directory = str(os.getenv("ZOL_EXEC_PATH")).strip()
conda_setup_success = None
plot_prog = None
tinyaai_prog = None
phylo_plot_prog = None
split_diamond_results_prog = None
if zol_exec_directory != 'None':
	try:
		zol_exec_directory = os.path.abspath(zol_exec_directory) + '/'
		plot_prog = zol_exec_directory + 'plotSegments.R'
		phylo_plot_prog = zol_exec_directory + 'phyloHeatmap.R'
		tinyaai_prog = zol_exec_directory + 'plotTinyAAI.R'
		split_diamond_results_prog = zol_exec_directory + 'splitDiamondResultsForFai'
		conda_setup_success = True
	except:
		conda_setup_success = False
if zol_exec_directory == 'None' or conda_setup_success == False:
	zol_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
	plot_prog = zol_main_directory + 'zol/plotSegments.R'
	tinyaai_prog = zol_main_directory + 'zol/plotTinyAAI.R'
	phylo_plot_prog = zol_main_directory + 'zol/phyloHeatmap.R'
	split_diamond_results_prog = os.path.abspath(os.path.dirname(__file__) + '/') + '/splitDiamondResultsForFai'
if plot_prog == None or phylo_plot_prog == None or split_diamond_results_prog == None or not os.path.isfile(plot_prog) or not os.path.isfile(phylo_plot_prog) or not os.path.isfile(split_diamond_results_prog) or not os.path.isfile(phylo_plot_prog):
	sys.stderr.write('Issues in setup of the zol-suite (in fai.py) - please describe your installation process and post an issue on GitHub!\n')
	sys.exit(1)

def subsetGenBankForQueryLocus(full_gw_genbank, locus_genbank, locus_proteins, reference_contig, reference_start, reference_end, logObject):
	"""
	Description:
	This function extracts a locus specific GenBank from a reference genome based on user provided coordinates.
	********************************************************************************************************************
	Parameters:
	- full_gw_genbank: Reference genome GenBank file - should have CDS features.
	- locus_genbank: The locus-specific GenBank file to create.
	- locus_proteins: The locus-specific protein FASTA file to create.
	- reference_contig: The identifier of the reference scaffold/contig with the locus.
	- reference_start: The starting position of the locus in the reference genome.
	- reference_end: The ending position of the locus in the reference genome.
	- logObject: A logging object.
	********************************************************************************************************************
	"""
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
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
def parseHomologGroupMatrix(orthogroup_matrix_file, logObject):
	"""
	Description:
	This function parses an ortholog group matrix and returns a dictionary mapping protein names/locus_tags to ortholog
	groups.
	********************************************************************************************************************
	Parameters:
	- orthogroup_matrix_file: The ortholog group vs sample matrix file, where cells correspond to locus tag identifiers.
	- logObject: A logging object.
	********************************************************************************************************************
	Return:
	- protein_to_hg: A dictionary mapping proteins names/locus_tags to ortholog groups.
	********************************************************************************************************************
	"""
	protein_to_hg = {}
	try:
		with open(orthogroup_matrix_file) as oog:
			for i, line in enumerate(oog):
				line = line.strip()
				ls = line.split('\t')
				if i == 0: continue
				hg = ls[0]
				for lts in ls[1:]:
					for lt in lts.split(', '):
						protein_to_hg[lt] = hg
	except Exception as e:
		sys.stderr.write('Issues mapping proteins to ortholog groups.\n')
		logObject.error('Issues mapping proteins to ortholog groups.')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return protein_to_hg

def collapseProteinsUsingCDHit(protein_fasta, nr_protein_fasta, logObject):
	"""
	Description:
	This function collapses query proteins using CD-HIT.
	********************************************************************************************************************
	Parameters:
	- protein_fasta: A file with the original set of protein queries in FASTA form.
	- nr_protein_fasta: The path to the non-redudnant (nr) protein queries after CD-HIT representative selection in
	                    FASTA format.
	- logObject: A logging object.
	********************************************************************************************************************
	Return:
	- protein_to_hg: A dictionary mapping proteins names/locus_tags to representative proteins.
	********************************************************************************************************************
	"""
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
			sys.stderr.write(traceback.format_exc())
			sys.exit(1)

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
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return protein_to_hg

def parseCoordsFromGenbank(genbanks, logObject):
	"""
	Description:
	This function gathers the coordinates for CDS features from multiple GenBank files.
	********************************************************************************************************************
	Parameters:
	- genbanks: A set of GenBank files.
	- logObject: A logging object.
	********************************************************************************************************************
	Return:
	- comp_gene_info: A 3-tier dicitonary with the first key corresponding to the gene cluster prefix/name, the second
	                  key corresponding to the locus tag identifier, and the third is information pertaining to CDS
	                  locations (scaffold, start position, end position).
	********************************************************************************************************************
	"""
	# get coords for individual genes
	comp_gene_info = defaultdict(dict)
	try:
		for gbk in genbanks:
			prefix = '.'.join(gbk.split('/')[-1].split('.')[:-1])
			gene_locations = util.parseGbk(gbk, prefix, logObject, use_either_lt_or_pi=True)
			for g in gene_locations:
				comp_gene_info[prefix][g] = gene_locations[g]
	except Exception as e:
		sys.stderr.write('Issues with parsing coordinates of genes in GenBanks!\n')
		logObject.error('Issues with parsing coordinates of genes in GenBanks!')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return dict(comp_gene_info)

def genConsensusSequences(genbanks, outdir, logObject, cpus=1, use_super5=False):
	"""
	Description:
	This function performs zol-based ortholog group determination across CDS protein sequences from query GenBank files
	provided and then determines consensus protein sequences for each ortholog group identified.
	********************************************************************************************************************
	Parameters:
	- genbanks: A set of query GenBank files.
	- outdir: The output directory/workspace where to write intermediate files/analyses.
	- logObject: A logging object.
	- cpus: The number of CPUs to use.
	- use_super5: Whether to use the SUPER5 algorithm for MUSCLE based alignment.
	********************************************************************************************************************
	Return:
	- orthogroup_matrix_file: An ortholog group vs sample matrix file, where cells correspond to locus tag identifiers.
	- consensus_prot_seqs_faa: Path to a FASTA file containing consensus protein sequences for ortholog groups
	                           determined.
	********************************************************************************************************************
	"""
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
				sys.stderr.write(traceback.format_exc())
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
			sys.stderr.write(traceback.format_exc())
			sys.exit(1)

	except Exception as e:
		sys.stderr.write('Issues with determining ortholog/ortholog groups!\n')
		logObject.error('Issues with determining ortholog/ortholog groups!')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

	# some of these folders are needed for the function (as it is used in zol) but just
	# created here to apease the function requirements.
	proc_dir = outdir + 'Homolog_Group_Processing/'
	hg_prot_dir = proc_dir + 'OG_Protein_Sequences/'
	hg_nucl_dir = proc_dir + 'OG_Nucleotide_Sequences/'
	prot_algn_dir = proc_dir + 'OG_Protein_Alignments/'
	phmm_dir = proc_dir + 'OG_Profile_HMMs/'
	cons_dir = proc_dir + 'OG_Consensus_Sequences/'
	consensus_prot_seqs_faa = outdir + 'OG_Consensus_Seqs.faa'
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

def loadTargetGenomeInfo(target_annotation_information, target_genomes_pkl_dir, diamond_reuslts, valid_tg_samples,
						 logObject, lowmem_mode=True,  min_genes_per_scaffold=1):
	"""
	Description:
	Load information pertaining to CDS locations from target genomes from sample pickle files created during prepTG.
	********************************************************************************************************************
	Parameters:
	- target_annotation_information: A set/dictionary of target genome sample identifiers.
	- target_genomes_pkl_dir: The directory where pickle files on target genome CDS location information are stored.
	- diamond_results: 
	- valid_tg_samples:
	- logObject: A logging object.
	- lowmem_mode: Whether to only regard CDS found on scaffolds with hits from the gene cluster queries.
	- min_genes_per_scaffold: The minimum number of ORFs hit by query proteins for the scaffold to be considered and
	                          CDS information for it stored.
	********************************************************************************************************************
	Return:
	- target_genome_info: A dictionary of dictionaries:
		- gene_locations: A multi-tiered dictionary with location information for CDS features.
		- scaffold_genes: A double dictionary with keys corresponding to samples and secondary keys to scaffolds and
		                  values to CDS features found on them.
		- boundary_genes: A dictionary with keys corresponding to samples and values to sets of CDS features nearby
		                  scaffold edges.
		- gene_id_to_order: A multi-tiered dictionary mapping CDS identifiers to their order along scaffolds.
		- gene_order_to_id: A multi-tiered dictionary mapping gene order indices along scaffolds to CDS identifiers.
	********************************************************************************************************************
	"""
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
				sys.stderr.write(traceback.format_exc())
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
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)


	target_genome_gene_info = {'gene_locations': gene_locations, 'scaffold_genes': scaffold_genes,
							   'boundary_genes': boundary_genes, 'gene_id_to_order': gene_id_to_order,
							   'gene_order_to_id': gene_order_to_id}
	return(target_genome_gene_info)

def runDiamondBlastp(target_concat_genome_db, query_fasta, diamond_work_dir, logObject,
					 diamond_sensitivity='very-sensitive', evalue_cutoff=1e-10, cpus=1, compute_query_coverage=False):
	"""
	Description:
	This functions runs DIAMOND blastp analysis in fai for assessing homology of target genome proteins to query
	proteins.
	********************************************************************************************************************
	Parameters:
	- target_concat_genome_db: The path to the concatenated FASTA file comprising proteomes from all target genomes.
	- query_fasta: The path to the query/reference FASTA file.
	- diamond_work_dir: The workspace directory where DIAMOND blastp alignment and intermediate files should be written.
	- logObject: A logging object.
	- diamond_sensitivity: The sensitivity mode to run DIAMOND blastp with.
	- evalue_cutoff: The maximum E-value cutoff to regard an alignment to a target genome protein as homologous to a
	                 query protein.
	- cpus: The number of CPUs to use.
	- compute_query_coverage: Whether to compute the query coverage - used for the simple BLASTp approach of abon,
	                          atpoc, and apos.
	********************************************************************************************************************
	Return:
	- diamond_results_file: Path to the DIAMOND blastp resulting alignment file.
	********************************************************************************************************************
	"""

	diamond_results_file = diamond_work_dir + 'DIAMOND_Results.txt'
	try:
		diamond_blastp_cmd = ['diamond', 'blastp', '--ignore-warnings', '--threads', str(cpus), '--' + diamond_sensitivity,
					   '--query', query_fasta, '--db', target_concat_genome_db, '--outfmt', '6', 'qseqid', 'sseqid',
					   'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
					   'bitscore', 'qlen', 'slen', '-k0', '--out', diamond_results_file, '--evalue', str(evalue_cutoff)]
		if compute_query_coverage:
			diamond_blastp_cmd = ['diamond', 'blastp', '--ignore-warnings', '--threads', str(cpus), '--' + diamond_sensitivity,
				'--query', query_fasta, '--db', target_concat_genome_db, '--outfmt', '6', 'qseqid', 'sseqid',
				'pident', 'evalue', 'bitscore', 'qlen', 'slen', 'qcovhsp', '-k0', '--out', diamond_results_file, 
				'--evalue', str(evalue_cutoff)]
		try:
			subprocess.call(' '.join(diamond_blastp_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			logObject.info('Successfully ran: %s' % ' '.join(diamond_blastp_cmd))
		except Exception as e:
			logObject.error('Had an issue running: %s' % ' '.join(diamond_blastp_cmd))
			sys.stderr.write('Had an issue running: %s' % ' '.join(diamond_blastp_cmd))
			logObject.error(e)
			sys.stderr.write(traceback.format_exc())
			sys.exit(1)
		assert(os.path.isfile(diamond_results_file) and os.path.getsize(diamond_results_file) > 0)

	except Exception as e:
		logObject.error('Issues with running DIAMOND blastp.')
		sys.stderr.write('Issues with running DIAMOND blastp.\n')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return(diamond_results_file)

def processDiamondBlastp(target_annot_information, diamond_results_file, work_dir, logObject):
	"""
	Description:
	This functions processes DIAMOND blastp results.
	********************************************************************************************************************
	Parameters:
	- target_annot_information: A multi-tiered dictionary with
	- diamond_results_file: The path to the query/reference FASTA file.
	- work_dir: The workspace directory where DIAMOND blastp alignment and intermediate files should be written.
	- logObject: A logging object.
	********************************************************************************************************************
	Return:
	- diamond_results: A dictionary with locus tags as keys and values as a list of 5 items:
		- hg: The best matching query/reference ortholog group.
		- best_hit_per_lt_bitscore: The bitscore of the best hit for the CDS.
		- best_hit_per_lt_eval: The E-value of the best hit for the CDS.
		- sample: The sample identifier.
		- best_hit_per_lt_identity: The identity of the best hit for the CDS.
		- best_hit_per_lt_sqlratio: The identity of the subject to query length ratio for the best hit for the CDS.
	********************************************************************************************************************
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
			sys.stderr.write(traceback.format_exc())
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
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

	return(dict(diamond_results))

def mapKeyProteinsToHomologGroups(query_fasta, key_protein_queries_fasta, work_dir, logObject, cpus=1):
	"""
	Description:
	Function to align key protein queries FASTA to general query ortholog FASTA file and determine which general
	query ortholog groups / non-redundant protein sequences correspond to the key ones provided by the user. The best
	hit for key proteins to general query sequences will be taken and assigned by bitscore (assuming an E-value
	theshold of 1e-10 is met).
	********************************************************************************************************************
	Parameters:
	- query_fasta: A FASTA file of query ortholog group (consensus) sequences.
	- key_protein_queries_fasta: A FASTA file of key sequences.
	- work_dir: The workspace directory where to perform mapping alignment and write intermediate files.
	- logObject: A logging object.
	- cpus: The number of CPUs to use.
	********************************************************************************************************************
	Return:
	- key_hgs: A set of "key" designated query ortholog groups.
	********************************************************************************************************************
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
			sys.stderr.write(traceback.format_exc())
			sys.exit(1)

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
					sys.stderr.write(traceback.format_exc())
					sys.exit(1)

	except Exception as e:
		logObject.error('Issues with mapping key proteins to general proteins or consensus sequence of ortholog groups.')
		sys.stderr.write('Issues with mapping key proteins to general proteins or consensus sequence of ortholog groups.\n')
		sys.stderr.write(str(e) + '\n')
		raise RuntimeError(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return key_hgs


gc_genbanks_dir, gc_info_dir, query_gene_info, lt_to_hg, model, model_labels, target_annotation_info, boundary_genes = [None]*8
sample_lt_to_evalue, sample_lt_to_identity, sample_lt_to_sqlratio, sample_lt_to_bitscore = [None]*4
lts_ordered_dict, hgs_ordered_dict, gene_locations, gene_id_to_order, gene_order_to_id = [None]*5

def identifyGCInstances(query_information, target_information, diamond_results, work_dir, logObject, min_hits=5,
						 min_key_hits=3, draft_mode=False, gc_to_gc_transition_prob=0.9, bg_to_bg_transition_prob=0.9,
						 gc_emission_prob_with_hit=0.95,  gc_emission_prob_without_hit=0.2,
						 syntenic_correlation_threshold=0.8, max_int_genes_for_merge=0,  kq_evalue_threshold=1e-20,
						 flanking_context=1000, cpus=1, block_size=3000, gc_delineation_mode="GENE-CLUMPER"):
	"""
	Description:
	This function sets up identification of homologous/orthologous instances of gene clusters in target genomes based
	on DIAMOND blastp results and location information of genes in the query and target genomes. It is void and produces
	two types of files in subdirectories of work_dir, one being an extracted instance of homologous gene-clusters from
	target genomes and the other a TSV with information for individual CDS genes on how well they are matching query
	sequence.
	********************************************************************************************************************
	Parameters:
	- query_information: A multi-tiered dictionary with query sequence information:
		- protein_to_hg: Dictionary mapping individual query sequences to ortholog groups.
		- query_fasta: The path to a FASTA of query sequences.
	    - comp_gene_info: Dictionary linking to query CDS location/syntenic information.
	- target_information: A multi-tiered dictionary with information on CDS feature locations across target genomes.
	- diamond_results: A dictionary with locus tags as keys and values as a list of 5 items:
		- hg: The best matching query/reference ortholog group.
		- best_hit_per_lt_bitscore: The bitscore of the best hit for the CDS.
		- best_hit_per_lt_eval: The E-value of the best hit for the CDS.
		- sample: The sample identifier.
		- best_hit_per_lt_identity: The identity of the best hit for the CDS.
		- best_hit_per_lt_sqlratio: The identity of the subject to query length ratio for the best hit for the CDS.
	- work_dir: The workspace where to write files to for the analysis/function.
	- logObject: A logging object.
	- min_hits: The minimum number of non-redundant or ortholog group query sequences needed to be hit for gene cluster
	            detection.
	- min_key_hits: The minimum number of key non-redundant or ortholog group query sequences needed to be hit for gene
	                cluster detection.
	- draft_mode: Whether to use draft mode and regard cutoffs in aggregate across candidate gene cluster instances
	              found near scaffold edges.
	- gc_to_gc_transition_prob: The gene cluster to gene cluster transition probability to set for HMM.
	- bg_to_bg_transition_prob: The background to background transition probability to set for HMM.
	- gc_emission_prob_with_hit: The gene cluster emission probability for the gene cluster state when a DIAMOND hit is
	                             present at a CDS/ORF.
	- gc_emission_prob_without_hit: The gene cluster emission probability for the gene cluster state  when no DIAMOND
	                                hit is present at a CDS/ORF.
	- syntenic_correlation_threshold: The threshold for the absolute syntenic correlation between the candidate gene
									  cluster homologous instance and the query gene cluster.
	- max_int_genes_for_merge: The maximum number of intermediate genes between those aligning to the query proteins to
	                           be regarded as co-located enough to correspond to a potential single instance of the
	                           gene cluster.
	- kq_evalue_threshold: The E-value threshold for "key" designated ortholog groups to be regarded as actually being
	                       hit.
	- flanking_context: The number of basepairs to take around identified gene clusters when extracting the gene cluster
	                    specific GenBank from the full target genome GenBank.
	- cpus: The number of CPUs to use.
	- block_size: The maximum number of samples to consider at a time. This might be good to make an adjustable option
	              for those interested in applying fai to metagenomes, where setting it to a much lower value likely
	              makes sense.
	- gc_delineation_mode: The method to use for gene-cluster delineation, can either be "GENE-CLUMPER" or "HMM"
	********************************************************************************************************************
	"""

	try:
		global gc_genbanks_dir, gc_info_dir, query_gene_info, lt_to_hg, model, model_labels, target_annotation_info, boundary_genes
		global single_query_mode, sample_lt_to_evalue, sample_lt_to_identity, sample_lt_to_sqlratio, sample_lt_to_bitscore
		global lts_ordered_dict, hgs_ordered_dict, gene_locations, gene_id_to_order, gene_order_to_id

		gc_genbanks_dir = os.path.abspath(work_dir + 'GeneCluster_Genbanks') + '/'
		gc_info_dir = os.path.abspath(work_dir + 'GeneCluster_Info') + '/'
		util.setupReadyDirectory([gc_genbanks_dir, gc_info_dir])

		# unpack information in dictionaries
		query_gene_info = query_information['comp_gene_info']
		single_query_mode = query_information['single_query_mode']
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

		# Create DenseHMM - using the pomegrenate library
		gc_hg_probs = [gc_emission_prob_without_hit]
		bg_hg_probs = [1.0-gc_emission_prob_without_hit]
		model_labels = ['background']
		for hg in all_hgs:
			model_labels.append(hg)
			gc_hg_probs.append(gc_emission_prob_with_hit)
			bg_hg_probs.append(1.0-gc_emission_prob_with_hit)

		gc_cat = Categorical([gc_hg_probs])
		bg_cat = Categorical([bg_hg_probs])

		model = DenseHMM()
		model.add_distributions([gc_cat, bg_cat])

		gc_to_gc = gc_to_gc_transition_prob
		gc_to_bg = 1.0 - gc_to_gc_transition_prob
		bg_to_bg = bg_to_bg_transition_prob
		bg_to_gc = 1.0 - bg_to_bg_transition_prob

		start_to_gc = 0.5
		start_to_bg = 0.5
		gc_to_end = 0.5
		bg_to_end = 0.5

		model.add_edge(model.start, gc_cat, start_to_gc)
		model.add_edge(model.start, bg_cat, start_to_bg)
		model.add_edge(gc_cat, model.end, gc_to_end)
		model.add_edge(bg_cat, model.end, bg_to_end)
		model.add_edge(gc_cat, gc_cat, gc_to_gc)
		model.add_edge(gc_cat, bg_cat, gc_to_bg)
		model.add_edge(bg_cat, gc_cat, bg_to_gc)
		model.add_edge(bg_cat, bg_cat, bg_to_bg)

		"""
		# Previous code for older versions of pomegrenate

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
		"""

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
				#if len(sample_hgs[sample]) < 3: continue
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
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def identify_gc_instances(input_args):
	"""
	Description:
	This function is the actual one which identifies gene-clusters and is separated from identifyGCInstances() to allow
	multiple processing. It is sample specific. It uses global variable set by identifyGCInstances().
	********************************************************************************************************************
	Parameters:
	- sample: The sample identifier.
	- min_hits: The minimum number of non-redundant or ortholog group query sequences needed to be hit for gene cluster
	            detection.
	- min_key_hits: The minimum number of key non-redundant or ortholog group query sequences needed to be hit for gene
	                cluster detection.
	- key_hgs: The set of key query ortholog group identifiers.
	- kq_evalue_threshold: The E-value threshold for "key" designated ortholog groups to be regarded as actually being
	                       hit.
	- syntenic_correlation_threshold: The threshold for the absolute syntenic correlation between the candidate gene
									  cluster homologous instance and the query gene cluster.
	- max_int_genes_for_merge: The maximum number of intermediate genes between those aligning to the query proteins to
	                           be regarded as co-located enough to correspond to a potential single instance of the
	                           gene cluster.
	- flanking_context: The number of basepairs to take around identified gene clusters when extracting the gene cluster
	                    specific GenBank from the full target genome GenBank.
	- draft_mode: Whether to use draft mode and regard cutoffs in aggregate across candidate gene cluster instances
	              found near scaffold edges.
	- gc_delineation_mode: The method to use for gene-cluster delineation, can either be "GENE-CLUMPER" or "HMM"
	********************************************************************************************************************
	"""
	sample, min_hits, min_key_hits, key_hgs, kq_evalue_threshold, syntenic_correlation_threshold, max_int_genes_for_merge, flanking_context, draft_mode, gc_delineation_mode = input_args
	bg_set = set(['background'])

	if single_query_mode:
		gc_sample_listing_handle = open(gc_info_dir + sample + '.bgcs.txt', 'w')
		gc_hg_evalue_handle = open(gc_info_dir + sample + '.hg_evalues.txt', 'w')

		sample_gc_id = 1
		for scaffold in hgs_ordered_dict[sample]:
			hgs_ordered = hgs_ordered_dict[sample][scaffold]
			lts_ordered = lts_ordered_dict[sample][scaffold]
			for i, hg in enumerate(hgs_ordered):
				if hg != 'background':
					lt = lts_ordered[i]
					gc_genbank_file = gc_genbanks_dir + sample + '_fai-gene-cluster-' + str(sample_gc_id) + '.gbk'
					sample_gc_id += 1

					min_gc_pos = gene_locations[sample][lt]['start'] - flanking_context
					max_gc_pos = gene_locations[sample][lt]['end'] + flanking_context

					util.createGenbank(target_annotation_info[sample]['genbank'], gc_genbank_file, scaffold,
									   min_gc_pos, max_gc_pos)
					gc_sample_listing_handle.write('\t'.join([sample, gc_genbank_file]) + '\n')

					identity = 0.0
					sqlratio = 0.0
					bitscore = 0.0
					if lt in sample_lt_to_identity[sample]: identity = sample_lt_to_identity[sample][lt]
					if lt in sample_lt_to_sqlratio[sample]: sqlratio = sample_lt_to_sqlratio[sample][lt]
					if lt in sample_lt_to_bitscore[sample]: bitscore = sample_lt_to_bitscore[sample][lt]
					gc_hg_evalue_handle.write('\t'.join([gc_genbank_file, sample, lt, hg, str(bitscore), str(identity),
								                         str(sqlratio), str(hg in key_hgs)]) + '\n')
					sample_gc_id += 1

		gc_hg_evalue_handle.close()
		gc_sample_listing_handle.close()
		return

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
						key_hgs_detected = set([])
						if len(key_hgs.intersection(set(gc_state_hgs).difference(bg_set))) > 0:
							for j, lt in enumerate(gc_state_lts):
								curr_hg = gc_state_hgs[j]
								if curr_hg in key_hgs and lt in sample_lt_to_evalue[sample] and sample_lt_to_evalue[sample][lt] <= kq_evalue_threshold:
									key_hgs_detected.add(curr_hg)
									features_key_hg = True
						if len(boundary_genes[sample].intersection(set(gc_state_lts).difference(bg_set))) > 0:
							boundary_lt_featured = True
						sample_gc_predictions.append([gc_state_lts, gc_state_hgs, len(gc_state_lts),
													   len(set(gc_state_hgs).difference(bg_set)),
													   len(key_hgs_detected),
													   scaffold, boundary_lt_featured, features_key_hg, gcs_id, key_hgs_detected])
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
				key_hgs_detected = set([])
				if len(key_hgs.intersection(set(gc_state_hgs).difference(bg_set))) > 0:
					for j, lt in enumerate(gc_state_lts):
						curr_hg = gc_state_hgs[j]
						if curr_hg in key_hgs and lt in sample_lt_to_evalue[sample] and sample_lt_to_evalue[sample][lt] <= kq_evalue_threshold:
							features_key_hg = True
							key_hgs_detected.add(curr_hg)
				if len(boundary_genes[sample].intersection(set(gc_state_lts).difference(bg_set))) > 0:
					boundary_lt_featured = True
				sample_gc_predictions.append([gc_state_lts, gc_state_hgs, len(gc_state_lts),
											  len(set(gc_state_hgs).difference(bg_set)),
											  len(key_hgs_detected),
											  scaffold, boundary_lt_featured, features_key_hg, gcs_id, key_hgs_detected])
				gcs_id += 1
	else:
		gcs_id = 1
		for scaffold in hgs_ordered_dict[sample]:
			hgs_ordered = hgs_ordered_dict[sample][scaffold]
			lts_ordered = lts_ordered_dict[sample][scaffold]

			hg_seq = numpy.array([[[model_labels.index(hg)] for hg in list(hgs_ordered)]])
			hmm_predictions = model.predict(hg_seq)[0]

			"""
			hg_seq = numpy.array(list(hgs_ordered))
			hmm_predictions = model.predict(hg_seq)
			"""

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
					if len(set(gc_state_hgs).difference(bg_set)) >= 3:
						boundary_lt_featured = False
						features_key_hg = False
						key_hgs_detected = set([])
						if len(key_hgs.intersection(set(gc_state_hgs).difference(bg_set))) > 0:
							for j, lt in enumerate(gc_state_lts):
								curr_hg = gc_state_hgs[j]
								if curr_hg in key_hgs and lt in sample_lt_to_evalue[sample] and sample_lt_to_evalue[sample][lt] <= kq_evalue_threshold:
									features_key_hg = True
									key_hgs_detected.add(curr_hg)
						if len(boundary_genes[sample].intersection(set(gc_state_lts).difference(bg_set))) > 0:
							boundary_lt_featured = True
						sample_gc_predictions.append([gc_state_lts, gc_state_hgs, len(gc_state_lts),
													   len(set(gc_state_hgs).difference(bg_set)),
													   len(key_hgs_detected),
													   scaffold, boundary_lt_featured, features_key_hg, gcs_id, key_hgs_detected])
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
					if len(set(gc_state_hgs).difference(bg_set)) >= 3:
						if '1' in set(gc_states):
							gc_state_lts = []
							gc_state_hgs = []
							gc_states = []
							continue
						boundary_lt_featured = False
						features_key_hg = False
						key_hgs_detected = set([])
						if len(key_hgs.intersection(set(gc_state_hgs).difference(bg_set))) > 0:
							for j, lt in enumerate(gc_state_lts):
								curr_hg = gc_state_hgs[j]
								if curr_hg in key_hgs and lt in sample_lt_to_evalue[sample] and sample_lt_to_evalue[sample][lt] <= kq_evalue_threshold:
									features_key_hg = True
									key_hgs_detected.add(curr_hg)
						if len(boundary_genes[sample].intersection(set(gc_state_lts).difference(bg_set))) > 0:
							boundary_lt_featured = True
						sample_gc_predictions.append([gc_state_lts, gc_state_hgs, len(gc_state_lts),
													   len(set(gc_state_hgs).difference(bg_set)),
													   len(key_hgs_detected),
													   scaffold, boundary_lt_featured, features_key_hg, gcs_id, key_hgs_detected])
						gcs_id += 1
					gc_state_lts = []
					gc_state_hgs = []
					gc_states = []
				if i == (len(hmm_predictions) - 1):
					if len(set(gc_state_hgs).difference(bg_set)) >= 3:
						if '1' in set(gc_states):
							gc_state_lts = []
							gc_state_hgs = []
							gc_states = []
							continue
						boundary_lt_featured = False
						features_key_hg = False
						key_hgs_detected = set([])
						if len(key_hgs.intersection(set(gc_state_hgs).difference(bg_set))) > 0:
							for j, lt in enumerate(gc_state_lts):
								curr_hg = gc_state_hgs[j]
								if curr_hg in key_hgs and lt in sample_lt_to_evalue[sample] and sample_lt_to_evalue[sample][lt] <= kq_evalue_threshold:
									features_key_hg = True
									key_hgs_detected.add(curr_hg)
						if len(boundary_genes[sample].intersection(set(gc_state_lts).difference(bg_set))) > 0:
							boundary_lt_featured = True
						sample_gc_predictions.append([gc_state_lts, gc_state_hgs, len(gc_state_lts),
													   len(set(gc_state_hgs).difference(bg_set)),
													   len(key_hgs_detected),
													   scaffold, boundary_lt_featured, features_key_hg, gcs_id, key_hgs_detected])
						gcs_id += 1
					gc_state_lts = []
					gc_state_hgs = []
					gc_states = []

	if len(sample_gc_predictions) == 0: return

	sorted_sample_gc_predictions = [x for x in sorted(sample_gc_predictions, key=itemgetter(3), reverse=True)]

	cumulative_edge_hgs = set([])
	cumulative_edge_key_hgs = set([])
	visited_scaffolds_with_edge_gc_segment = set([])
	sample_gc_predictions_filtered = []
	sample_edge_gc_predictions_filtered = []

	for gc_segment in sorted_sample_gc_predictions:
		if (gc_segment[3] >= min_hits and gc_segment[4] >= min_key_hits) or (gc_segment[3] >= 3 and gc_segment[6] and not gc_segment[5] in visited_scaffolds_with_edge_gc_segment):
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
						sys.stderr.write(traceback.format_exc())
						sys.exit(1)
						#pass
			if best_corr != 'irrelevant' and (best_corr == None or best_corr < syntenic_correlation_threshold): continue
			if (gc_segment[3] >= min_hits and gc_segment[4] >= min_key_hits):
				sample_gc_predictions_filtered.append(gc_segment + [best_corr])
				if gc_segment[6]:
					cumulative_edge_hgs = cumulative_edge_hgs.union(set(gc_segment[1]))
					cumulative_edge_key_hgs = cumulative_edge_key_hgs.union(set(gc_segment[9]))
					visited_scaffolds_with_edge_gc_segment.add(gc_segment[5])
			elif gc_segment[3] >= 3 and gc_segment[6] and not gc_segment[5] in visited_scaffolds_with_edge_gc_segment:
				sample_edge_gc_predictions_filtered.append(gc_segment + [best_corr])
				visited_scaffolds_with_edge_gc_segment.add(gc_segment[5])
				cumulative_edge_hgs = cumulative_edge_hgs.union(set(gc_segment[1]))
				cumulative_edge_key_hgs = cumulative_edge_key_hgs.union(set(gc_segment[9]))
					
	if len(sample_edge_gc_predictions_filtered) >= 1 and draft_mode:
		if len(cumulative_edge_hgs) >= min_hits and len(cumulative_edge_key_hgs) >= min_key_hits:
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
	gc_corr_info_handle = open(gc_info_dir + sample + '.corr_info.txt', 'w')

	sample_gc_id = 1
	for gc_segment in dereplicated_sample_gc_predictions_filtered:
		gc_genbank_file = gc_genbanks_dir + sample + '_fai-gene-cluster-' + str(sample_gc_id) + '.gbk'
		sample_gc_id += 1

		gc_segment_scaff = gc_segment[5]
		min_gc_pos = min([gene_locations[sample][g]['start'] for g in gc_segment[0]])-flanking_context
		max_gc_pos = max([gene_locations[sample][g]['end'] for g in gc_segment[0]])+flanking_context

		util.createGenbank(target_annotation_info[sample]['genbank'], gc_genbank_file, gc_segment_scaff,
							  min_gc_pos, max_gc_pos)
		gc_sample_listing_handle.write('\t'.join([sample, gc_genbank_file]) + '\n')
		gc_corr_info_handle.write('\t'.join([sample, gc_genbank_file, str(gc_segment[10])]) + '\n')

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
	gc_corr_info_handle.close()

def filterParalogousSegmentsAndConcatenateIntoMultiRecordGenBanks(hmm_work_dir, homologous_gbk_dir, single_query_mode, logObject):
	"""
	Description:
	This function allows for filtering paralogous gene-cluster segments identified in target genomes where the segments
	have overlap in the query non-redundant/ortholog groups they are mapping to. Gene cluster instances are considered
	paralogous if at least two overlapping query hits are found. Gene-clusters retained are prioritized based on
	cumulative bitscores to query genes.
	********************************************************************************************************************
	Parameters:
	- hmm_work_dir: The workspace where identifyGCInstances() saved extracted GenBanks for gene cluster instances from
	                target genomes along with homology information to the query gene cluster.
	- homologous_gbk_dir: The directory where to save the final extracted gene-cluster files.
	- logObject: A logging object.
	********************************************************************************************************************
	"""
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
					# consider segments paralogous if more than 2 reference proteins/ortholog groups are overlapping
					# suggesting paralogy beyond fragmentation that might have split a gene in two.
					# This is not the case if single_query_mode is True.
					intersection_hgs = gcs1_hg.intersection(gcs2_hg)
					if len(intersection_hgs) >= 2 or (single_query_mode and len(intersection_hgs) >= 1):
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
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def plotOverviews(target_annotation_info, hmm_work_dir, protein_to_hg, plot_work_dir, plot_result_pdf, plot_naming_pdf,
				  logObject, height=10, width=15):
	"""
	Description:
	This function serves to create plots from fai analysis showing the similarity of homologous / orthologous
	gene-cluster instances identified to the query gene cluster proteins.
	********************************************************************************************************************
	Parameters:
	- target_annotation_info: Dictionary linking sample names to their respective full genome GenBank files.
	- hmm_work_dir: The workspace where identifyGCInstances() saved extracted GenBanks for gene cluster instances from
	     			target genomes along with homology information to the query gene cluster.
	- protein_to_hg: A dictionary linking query proteins to non-redundant clusters / ortholog groups.
	- plot_work_dir: The workspace where to write intermediate files used for creating plots.
	- plot_result_pdf: A final plot (PDF) showcasing homology of target genome gene clusters to the query
	                   gene cluster(s).
	- plot_naming_pdf: A plot (PDF) based table which shows the mapping of query protein identifiers to ortholog group
	                   or non-redudnant cluster identifiers.
	- logObject: A logging object.
	- height: The height of the plot in inches.
	- width: The width of the plot in inches.
	********************************************************************************************************************
	"""
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
						lt = None
						try:
							lt = feature.qualifiers.get('locus_tag')[0]
						except:
							try:
								lt = feature.qualifiers.get('protein_id')[0]
							except:
								sys.stderr.write('The GenBank %s, cataloging a homologous instance to the query gene cluster had at least one CDS without either a locus_tag or protein_id feature.\n' % target_annotation_info[sample]['genbank'])
								sys.exit(1)
						if lt in relevant_lts:
							scaff_is_relevant = True
							break
					if not scaff_is_relevant: continue
					for feature in rec.features:
						if feature.type != 'CDS': continue
						lt = None
						try:
							lt = feature.qualifiers.get('locus_tag')[0]
						except:
							try:
								lt = feature.qualifiers.get('protein_id')[0]
							except:
								sys.stderr.write('The GenBank %s, cataloging a homologous instance to the query gene cluster had at least one CDS without either a locus_tag or protein_id feature.\n' % target_annotation_info[sample]['genbank'])
								sys.exit(1)
						all_coords = []
						if not 'join' in str(feature.location):
							start = min([int(x.strip('>').strip('<')) for x in
										 str(feature.location)[1:].split(']')[0].split(':')]) + 1
							end = max([int(x.strip('>').strip('<')) for x in
									   str(feature.location)[1:].split(']')[0].split(':')])
							direction = str(feature.location).split('(')[1].split(')')[0]
							all_coords.append([start, end, direction])
						elif 'order' in str(feature.location):
							for exon_coord in str(feature.location)[6:-1].split(', '):
								start = min(
									[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
								end = max(
									[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
								direction = exon_coord.split('(')[1].split(')')[0]
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
			sys.stderr.write(traceback.format_exc())
			sys.exit(1)

	except Exception as e:
		logObject.error('Issues with plotting overviews of homologous gene-cluster segments identified.')
		sys.stderr.write('Issues with plotting overviews of homologous gene-cluster segments identified.\n')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def plotTreeHeatmap(homologous_gbk_dir, hmm_work_dir, species_tree, plot_phylo_dir, plot_result_pdf, logObject, height=10, width=10):
	"""
	Description:
	This function serves to create a phylo-heatmap plot when/if a species tree is provided to fai. Takes 
	the best bitscore available for a query protein/ortholog-group per sample in a locus deemed homologous.
	********************************************************************************************************************
	Parameters:
	- homologous_gbk_dir: Results directory with homologous gene cluster GenBank files.
	- hmm_work_dir: Directory with itermediate files and information from HMM search for homologous gene cluster 
	                instances.
	- species_tree: Species tree in Newick format.
	- plot_phylo_dir: Directory with intermediate files for plotting phylo-heatmap.
	- plot_pdf: The final PDF result of the phylo-heatmap.
	- logObject: A logging object.
	- height: The height of the plot in inches.
	- width: The width of the plot in inches.
	********************************************************************************************************************
	"""
	try:
		heatmap_info_file = plot_phylo_dir + 'Heatmap_Track.txt'
		heatmap_info_file_handle = open(heatmap_info_file, 'w')
		heatmap_info_file_handle.write('label\tquery_prot_id\tbitscore\n')

		sample_final_lts = defaultdict(set) 
		for gbk in os.listdir(homologous_gbk_dir):
			sample = '.gbk'.join(gbk.split('.gbk')[:-1]).split('_fai-gene-cluster-')[0]
			with open(homologous_gbk_dir + gbk) as ogf:
				for rec in SeqIO.parse(ogf, 'genbank'):
					for feature in rec.features:
						if feature.type != 'CDS': continue
						lt = None
						try:
							lt = feature.qualifiers.get('locus_tag')[0]
						except:
							sys.stderr.write('The GenBank %s, cataloging a homologous instance to the query gene cluster had at least one CDS without a locus_tag feature.\n' % target_annotation_info[sample]['genbank'])
							sys.exit(1)
						sample_final_lts[sample].add(lt)
		if len(sample_final_lts) == 0:
			logObject.warning('Unable to generate phylogenetic-heatmap because no gene cluster instances were detected in target genomes.')
			sys.stderr.write('Warning: unable to generate phylogenetic-heatmap because no gene cluster instances were detected in target genomes.\n')
			return
		
		t = Tree(species_tree)
		all_samples_in_tree = set([])	
		samples_accounted = set([])
		for node in t.traverse('postorder'):
			if node.is_leaf and node.name != "":
				all_samples_in_tree.add(node.name)
			
		if len(all_samples_in_tree.intersection(set(sample_final_lts.keys()))) == 0:
			logObject.warning('Unable to generate phylogenetic-heatmap because species tree provided doesn\'t match target genomes searched against.')
			sys.stderr.write('Warning: Unable to generate phylogenetic-heatmap because species tree provided doesn\'t match target genomes searched against.\n')
			return
		
		gbk_info_dir = hmm_work_dir + 'GeneCluster_Info/'
		example_og = None
		if os.path.isdir(gbk_info_dir):
			for gcs_info in os.listdir(gbk_info_dir):
				if not gcs_info.endswith('.hg_evalues.txt') or os.path.getsize(gbk_info_dir + gcs_info) == 0: continue
				sample = gcs_info.split('.hg_evalues.txt')[0]
				if not sample in all_samples_in_tree:
					sys.stderr.write('Warning: sample %s not in phylogeny or has an unexpected name alteration.' % sample)
					logObject.warning('Sample %s not in phylogeny or has an unexpected name alteration.' % sample)
					continue
				samples_accounted.add(sample)
				with open(gbk_info_dir + gcs_info) as ogif:
					for line in ogif:
						line = line.strip()
						ls = line.split('\t')
						lt = ls[2]
						if not lt in sample_final_lts[sample]: continue
						og = ls[3]
						if og == 'background': continue
						example_og = og
						bitscore = float(ls[4])
						heatmap_info_file_handle.write(sample + '\t' + og + '\t' + str(bitscore) + '\n')

		if example_og == None:
			logObject.warning('Unable to generate phylogenetic-heatmap because no gene cluster instances were detected in target genomes.')
			sys.stderr.write('Warning: unable to generate phylogenetic-heatmap because no gene cluster instances were detected in target genomes.\n')
			return
			
		for sample in all_samples_in_tree:
			if not sample in samples_accounted:
				heatmap_info_file_handle.write(sample + '\t' + example_og + '\tNA\n')
		heatmap_info_file_handle.close()

		plot_cmd = ['Rscript', phylo_plot_prog, species_tree, heatmap_info_file, plot_result_pdf, str(height), str(width)]
		try:
			subprocess.call(' '.join(plot_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(plot_result_pdf))
			logObject.info('Successfully ran: %s' % ' '.join(plot_cmd))
		except Exception as e:
			logObject.error('Had an issue running R based phylo-heatmap plotting: %s' % ' '.join(plot_cmd))
			sys.stderr.write('Had an issue running R based phylo-heatmap plotting: %s\n' % ' '.join(plot_cmd))
			logObject.error(e)
			sys.stderr.write(traceback.format_exc())
			sys.exit(1)
				
	except Exception as e:
		logObject.error('Issues with creating phylo-heatmap.')
		sys.stderr.write('Issues with creating phylo-heatmap.\n')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def createOverviewSpreadsheetAndTinyAAIPlot(hmm_work_dir, protein_to_hg, key_hgs, spreadsheet_result_file, spreadsheet_result_tsv_dir, homologous_gbk_dir, plot_work_dir, tiny_aai_plot_file, logObject):
	"""
	Description:
	This function serves to create an overview spreadsheet detailing the candidate homologous gene cluster segments 
	identified. It will also generate a tiny AAI vs. genes-shared plot.
	********************************************************************************************************************
	Parameters:
	- hmm_work_dir: The workspace where identifyGCInstances() saved extracted GenBanks for gene cluster instances from
	     			target genomes along with homology information to the query gene cluster.
	- protein_to_hg: A dictionary linking query proteins to non-redundant clusters / ortholog groups.
	- key_hgs: Key homolog groups as defined by user.
	- spreadsheet_result_file: The path to the spreadsheet XLSX file with the overview info. 
	- spreadsheet_result_tsv_dir: Directory where to write the TSV files that ultimately go into the final spreadsheet.
	- plot_work_dir: The workspace where to write intermediate files used for creating plots.
	- tiny_aai_plot_file: A final plot (PDF) showcasing average amino-acid identity of homologous gene clusters identified
	                      to proportion of gene clusters represented.
	- logObject: A logging object.
	********************************************************************************************************************
	"""
	try:
		gbk_info_dir = hmm_work_dir + 'GeneCluster_Info/'

		hg_prots = defaultdict(set)	
		hg_list = []
		hg_headers = []
		for p in protein_to_hg:
			hg_prots[protein_to_hg[p]].add(p)
			hg_list.append(protein_to_hg[p])
			hg = protein_to_hg[p]
			hg_name = hg
			if hg in key_hgs:
				hg_name = 'KEY-PROTEIN_' + protein_to_hg[p]
			hg_headers.append(hg_name + ' - PID')
			hg_headers.append(hg_name + ' - SQR')

		gcgbk_to_corr = {}
		for f in os.listdir(gbk_info_dir):
			info_file = gbk_info_dir + f
			if f.endswith('.corr_info.txt') and os.path.getsize(info_file) > 0:
				with open(info_file) as oif:
					for line in oif:
						line = line.strip()
						sample, gc_genbank_file, corr_value = line.split('\t')
						if corr_value == 'irrelevant':
							gcgbk_to_corr[gc_genbank_file] = float('nan')
						else:
							gcgbk_to_corr[gc_genbank_file] = round(float(corr_value), 3)

		tiny_aai_plot_input = plot_work_dir + 'Tiny_AAI_Plot_Data.txt'
		tiny_aai_plot_handle = open(tiny_aai_plot_input, 'w')
		tiny_aai_plot_handle.write('Sample\tAAI\tProp_Genes_Found\tMean_Syntenic_Correlation\n')

		tot_tsv_file = spreadsheet_result_tsv_dir + 'total_gcs.tsv'
		igc_tsv_file = spreadsheet_result_tsv_dir + 'individual_gcs.tsv'

		tot_tsv_handle = open(tot_tsv_file, 'w')
		igc_tsv_handle = open(igc_tsv_file, 'w')

		data_tot_header = ['sample', 'aggregate-bitscore', 'aai-to-query', 'mean-sequence-to-query-ratio', 'proportion-query-genes-found', 
					 	   'avg-syntenic-correlation', 'number-background-genes', 'number-gene-clusters', 'copy-counts'] + hg_headers
		data_igc_header = ['sample', 'gene-cluster-id', 'aggregate-bitscore', 'aai-to-query', 'mean-sequence-to-query-ratio', 'proportion-query-genes-found', 
					 	   'syntenic-correlation', 'number-background-genes', 'copy-counts'] + hg_headers

		tot_tsv_handle.write('\t'.join(data_tot_header) + '\n')
		igc_tsv_handle.write('\t'.join(data_igc_header) + '\n')

		tot_row_count = 0
		igc_row_count = 0

		igc_aggregate_bitscores = []
		aggregate_bitscores = []
		for f in os.listdir(gbk_info_dir):
			info_file = gbk_info_dir + f
			if f.endswith('.hg_evalues.txt') and os.path.getsize(info_file) > 0:
				tot_hg_best_hits = defaultdict(lambda: [0.0, [], []])
				tot_hg_hits_copy_counts = defaultdict(int)
				tot_bg_genes = 0
				igc_hg_best_hits = defaultdict(lambda: defaultdict(lambda: [0.0, [], []]))
				igc_hg_hits_copy_counts = defaultdict(lambda: defaultdict(int))
				igc_bg_genes = defaultdict(int)
				tot_corr_values = []
				sample = None
				igcs = set([])
				with open(info_file) as oif:
					for line in oif:
						line = line.strip('\n')
						igc, sample, lt, hg, bitscore, identity, sqlratio, is_key_hg = line.split('\t')
						bitscore = float(bitscore)
						igcs.add(igc)
						tot_corr_values.append(gcgbk_to_corr[igc])
						if hg == 'background':
							tot_bg_genes += 1
							igc_bg_genes[igc] += 1
						else:
							tot_hg_hits_copy_counts[hg] += 1
							igc_hg_hits_copy_counts[igc][hg] += 1

							if tot_hg_best_hits[hg][0] < bitscore:
								identity = float(identity)
								sqlratio = float(sqlratio)
								tot_hg_best_hits[hg] = [bitscore, [identity], [sqlratio]]
							elif tot_hg_best_hits[hg][0] == bitscore:
								identity = float(identity)
								sqlratio = float(sqlratio)
								tot_hg_best_hits[hg][1].append(identity)
								tot_hg_best_hits[hg][2].append(sqlratio)

							if igc_hg_best_hits[igc][hg][0] < bitscore:
								identity = float(identity)
								sqlratio = float(sqlratio)
								igc_hg_best_hits[igc][hg] = [bitscore, [identity], [sqlratio]]
							elif igc_hg_best_hits[igc][hg][0] == bitscore:
								identity = float(identity)
								sqlratio = float(sqlratio)
								igc_hg_best_hits[igc][hg][1].append(identity)
								igc_hg_best_hits[igc][hg][2].append(sqlratio)

				copy_counts = []
				hg_pid_info = []
				hg_sql_info = []
				hg_pid_sql_pair = []
				aggregate_bitscore = 0.0
				for hg in hg_list:
					max_pid = 0.0
					max_sql = 0.0
					aggregate_bitscore += tot_hg_best_hits[hg][0]
					for i, pid in enumerate(tot_hg_best_hits[hg][1]):
						if pid >= max_pid:
							max_pid = pid
							max_sql = tot_hg_best_hits[hg][2][i]
					hg_pid_info.append(max_pid)
					hg_sql_info.append(max_sql)
					hg_pid_sql_pair.append(str(round(max_pid,3)) + '\t' + str(round(max_sql,3)))
					copy_counts.append(str(tot_hg_hits_copy_counts[hg]))
			
				mean_aai = statistics.mean([x for x in hg_pid_info if x > 0.0])
				mean_sql = statistics.mean([x for x in hg_sql_info if x > 0.0])
				prop_hg_found = len([x for x in hg_pid_info if x > 0.0])/float(len(hg_list))
				copy_count_string = ','.join(copy_counts)
				avg_corr_value = statistics.mean(tot_corr_values)
				number_gcs = len(igcs)
				
				tot_row = [str(x) for x in [sample, round(aggregate_bitscore, 3), round(mean_aai, 3), round(mean_sql, 3), round(prop_hg_found, 3), round(avg_corr_value, 3), tot_bg_genes, number_gcs, copy_count_string]] + hg_pid_sql_pair
				aggregate_bitscores.append(aggregate_bitscore)
				tot_tsv_handle.write('\t'.join(tot_row) + '\n')
				tiny_aai_plot_handle.write('\t'.join([str(x) for x in [sample, mean_aai, prop_hg_found, avg_corr_value]]) + '\n')
				tot_row_count += 1

				for igc in igcs:
					igc_copy_counts = []
					igc_hg_pid_info = []
					igc_hg_sql_info = []
					igc_hg_pid_sql_pair = []
					igc_aggregate_bitscore = 0.0
					for hg in hg_list:
						max_pid = 0.0
						max_sql = 0.0
						igc_aggregate_bitscore += igc_hg_best_hits[igc][hg][0]
						for i, pid in enumerate(igc_hg_best_hits[igc][hg][1]):
							if pid >= max_pid:
								max_pid = pid
								max_sql = igc_hg_best_hits[igc][hg][2][i]
						igc_hg_pid_info.append(max_pid)
						igc_hg_sql_info.append(max_sql)
						igc_hg_pid_sql_pair.append(str(round(max_pid, 3)) + '\t' + str(round(max_sql, 3)))
						igc_copy_counts.append(str(igc_hg_hits_copy_counts[igc][hg]))
			
					igc_mean_aai = statistics.mean([x for x in igc_hg_pid_info if x > 0.0])
					igc_mean_sql = statistics.mean([x for x in igc_hg_sql_info if x > 0.0])
					igc_prop_hg_found = len([x for x in igc_hg_pid_info if x > 0.0])/float(len(hg_list))
					copy_count_string = ','.join(igc_copy_counts)

					igc_path = igc
					if not os.path.isfile(igc):
						igc_path = homologous_gbk_dir + igc.split('/')[-1]
						try:
							assert(os.path.isfile(igc_path))
						except:
							sys.stderr.write('Warning: Gene cluster GenBank file not found/confirmed for: %s\n' % igc)
							logObject.warning('Gene cluster GenBank file not found/confirmed for: %s' % igc)

					igc_row = [str(x) for x in [sample, igc_path, round(igc_aggregate_bitscore, 3), round(igc_mean_aai, 3), round(igc_mean_sql, 3), 
								igc_prop_hg_found, gcgbk_to_corr[igc], igc_bg_genes[igc], copy_count_string]] + igc_hg_pid_sql_pair
					igc_aggregate_bitscores.append(igc_aggregate_bitscore)
					igc_tsv_handle.write('\t'.join(igc_row) + '\n')		
					igc_row_count += 1

		tot_tsv_handle.close()
		igc_tsv_handle.close()

		# Generate Excel spreadsheet
		writer = pd.ExcelWriter(spreadsheet_result_file, engine='xlsxwriter')
		workbook = writer.book
		dd_sheet = workbook.add_worksheet('Data Dictionary')
		dd_sheet.write(0, 0, 'Data Dictionary describing columns of "Overview" spreadsheets can be found on fai\'s Wiki page at:')
		dd_sheet.write(1, 0, 'https://github.com/Kalan-Lab/zol/wiki/3.-more-info-on-fai#explanation-of-report')

		numeric_columns = set(['aggregate-bitscore', 'aai-to-query', 'mean-sequence-to-query-ratio', 'proportion-query-genes-found', 
					 	   		'syntenic-correlation', 'avg-syntenic-correlation', 'number-background-genes', 'number-gene-clusters'] + hg_headers)

		header_format = workbook.add_format({'bold': True, 'text_wrap': True, 'valign': 'top', 'fg_color': '#D7E4BC', 'border': 1})

		tot_results_df = util.loadTableInPandaDataFrame(tot_tsv_file, numeric_columns)
		tot_results_df.to_excel(writer, sheet_name='Genome Wide - Report', index=False, na_rep="NA")

		igc_results_df = util.loadTableInPandaDataFrame(igc_tsv_file, numeric_columns)
		igc_results_df.to_excel(writer, sheet_name='Gene Cluster Instance - Report', index=False, na_rep="NA")

		tot_worksheet = writer.sheets['Genome Wide - Report']
		igc_worksheet = writer.sheets['Gene Cluster Instance - Report']

		tot_worksheet.conditional_format('A1:HA1', {'type': 'cell', 'criteria': '!=', 'value': 'NA', 'format': header_format})
		igc_worksheet.conditional_format('A1:HA1', {'type': 'cell', 'criteria': '!=', 'value': 'NA', 'format': header_format})

		#data_tot_header = ['sample', 'aggregate-bitscore', 'mean-aai-to-query', 'mean-sequence-to-query-ratio', 'proportion-query-genes-found', 
		#    			 	'avg-syntenic-correlation', 'number-background-genes', 'number-gene-clusters', 'copy-counts'] + hg_headers

		#data_igc_header = ['sample', 'gene-cluster-id', 'aggregate-bitscore', 'mean-aai-to-query', 'mean-sequence-to-query-ratio', 'proportion-query-genes-found', 
		#			 	    'syntenic-correlation', 'number-background-genes', 'copy-counts'] + hg_headers

		tot_row_count += 1
		igc_row_count += 1

		if len(aggregate_bitscores) > 0:
			# aggregate bitscore coloring
			tot_worksheet.conditional_format('B2:B' + str(tot_row_count), {'type': '2_color_scale', 'min_color': "#d9d9d9", 'max_color': "#949494", 'min_type': 'num', 'max_type': 'num', "min_value": min(aggregate_bitscores), "max_value": max(aggregate_bitscores)})
			igc_worksheet.conditional_format('C2:C' + str(igc_row_count), {'type': '2_color_scale', 'min_color': "#d9d9d9", 'max_color': "#949494", 'min_type': 'num', 'max_type': 'num', "min_value": min(igc_aggregate_bitscores), "max_value": max(igc_aggregate_bitscores)})

			# mean aai 
			tot_worksheet.conditional_format('C2:C' + str(tot_row_count), {'type': '2_color_scale', 'min_color': "#d8eaf0", 'max_color': "#83c1d4", "min_value": 0.0, "max_value": 100.0, 'min_type': 'num', 'max_type': 'num'})
			igc_worksheet.conditional_format('D2:D' + str(igc_row_count), {'type': '2_color_scale', 'min_color': "#d8eaf0", 'max_color': "#83c1d4", "min_value": 0.0, "max_value": 100.0, 'min_type': 'num', 'max_type': 'num'})

			# mean sequence-to-query ratio 
			tot_worksheet.conditional_format('D2:D' + str(tot_row_count), {'type': '2_color_scale', 'min_color': "#dbd5e8", 'max_color': "#afa1cf", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})
			igc_worksheet.conditional_format('E2:E' + str(igc_row_count), {'type': '2_color_scale', 'min_color': "#dbd5e8", 'max_color': "#afa1cf", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})

			# proportion query genes found
			tot_worksheet.conditional_format('E2:E' + str(tot_row_count), {'type': '2_color_scale', 'min_color': "#fce8ee", 'max_color': "#eb8da9", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})
			igc_worksheet.conditional_format('F2:F' + str(igc_row_count), {'type': '2_color_scale', 'min_color': "#fce8ee", 'max_color': "#eb8da9", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})

			# syntenic-correlation
			tot_worksheet.conditional_format('F2:F' + str(tot_row_count), {'type': '2_color_scale', 'min_color': "#fcf6eb", 'max_color': "#d1c0a1", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})
			igc_worksheet.conditional_format('G2:G' + str(igc_row_count), {'type': '2_color_scale', 'min_color': "#fcf6eb", 'max_color': "#d1c0a1", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})

			hg_index = 9
			for hg in hg_list:
				columnid_pid = util.determineColumnNameBasedOnIndex(hg_index)
				columnid_sql = util.determineColumnNameBasedOnIndex(hg_index+1)
				hg_index += 2

				# percent identity
				tot_worksheet.conditional_format(columnid_pid + '2:' + columnid_pid + str(tot_row_count), {'type': '2_color_scale', 'min_color': "#d8eaf0", 'max_color': "#83c1d4", "min_value": 0.0, "max_value": 100.0, 'min_type': 'num', 'max_type': 'num'})
				igc_worksheet.conditional_format(columnid_pid + '2:' + columnid_pid + str(tot_row_count), {'type': '2_color_scale', 'min_color': "#d8eaf0", 'max_color': "#83c1d4", "min_value": 0.0, "max_value": 100.0, 'min_type': 'num', 'max_type': 'num'})

				# sequence-to-query ratio
				tot_worksheet.conditional_format(columnid_sql + '2:' + columnid_sql + str(igc_row_count), {'type': '2_color_scale', 'min_color': "#dbd5e8", 'max_color': "#afa1cf", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})
				igc_worksheet.conditional_format(columnid_sql + '2:' + columnid_sql + str(igc_row_count), {'type': '2_color_scale', 'min_color': "#dbd5e8", 'max_color': "#afa1cf", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})

		else:
			tot_worksheet.write(1, 0, 'No hits to query gene cluster in target genomes found at requested thresholds! Try loosening parameters perhaps.')
			igc_worksheet.write(1, 0, 'No hits to query gene cluster in target genomes found at requested thresholds! Try loosening parameters perhaps.')

		# Freeze the first row of both sheets
		tot_worksheet.freeze_panes(1, 0)
		igc_worksheet.freeze_panes(1, 0)

		# close workbook
		workbook.close()
		
		tiny_aai_plot_handle.close()

		# plot tiny AAI figure
		plot_cmd = ['Rscript', tinyaai_prog, tiny_aai_plot_input, tiny_aai_plot_file]
		try:
			subprocess.call(' '.join(plot_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
			assert (os.path.isfile(tiny_aai_plot_file))
			logObject.info('Successfully ran: %s' % ' '.join(plot_cmd))
		except Exception as e:
			logObject.error('Had an issue running R based plotting: %s' % ' '.join(plot_cmd))
			sys.stderr.write('Had an issue running R based plotting: %s\n' % ' '.join(plot_cmd))
			logObject.error(e)
			sys.stderr.write(traceback.format_exc())
			sys.exit(1)

	except Exception as e:
		logObject.error('Issues with producing spreadsheet/AAI plot overviews for homologous gene-cluster segments identified.')
		sys.stderr.write('Issues with producing spreadsheet/AAI plot overviews of homologous gene-cluster segments identified.\n')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
