import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import logging
import subprocess
from operator import itemgetter
from collections import defaultdict
import traceback
from scipy import stats
from ete3 import Tree
import itertools
import math
import numpy as np
import gzip
import pathlib

valid_alleles = set(['A', 'C', 'G', 'T'])
curr_dir = os.path.abspath(pathlib.Path(__file__).parent.resolve()) + '/'
main_dir = '/'.join(curr_dir.split('/')[:-2]) + '/'


def multiProcess(input):
	"""
	Genralizable function to be used with multiprocessing to parallelize list of commands. Inputs should correspond
	to space separated command (as list), with last item in list corresponding to a logging object handle for logging
	progress.
	"""
	input_cmd = input[:-1]
	logObject = input[-1]
	logObject.info('Running the following command: %s' % ' '.join(input_cmd))
	try:
		subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(input_cmd))
	except Exception as e:
		logObject.error('Had an issue running: %s' % ' '.join(input_cmd))
		sys.stderr.write('Had an issue running: %s' % ' '.join(input_cmd))
		logObject.error(e)
		sys.exit(1)

def setupReadyDirectory(directories):
	try:
		assert (type(directories) is list)
		for d in directories:
			if os.path.isdir(d):
				os.system('rm -rf %s' % d)
			os.system('mkdir %s' % d)
	except Exception as e:
		raise RuntimeError(traceback.format_exc())


def p_adjust_bh(p):
	"""
	Benjamini-Hochberg p-value correction for multiple hypothesis testing.
	"""
	p = np.asfarray(p)
	by_descend = p.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(p)) / np.arange(len(p), 0, -1)
	q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
	return q[by_orig]


def is_newick(newick):
	"""
	Function to validate if Newick phylogeny file is correctly formatted.
	"""
	try:
		t = Tree(newick)
		return True
	except:
		return False


def is_fastq(fastq):
	"""
	Function to validate if FASTA file is correctly formatted.
	"""
	try:
		with open(fastq) as of:
			SeqIO.parse(of, 'fastq')
		return True
	except:
		return False


def is_fasta(fasta):
	"""
	Function to validate if FASTA file is correctly formatted.
	"""
	try:
		if fasta.endswith('.gz'):
			with gzip.open(fasta, 'rt') as ogf:
				SeqIO.parse(ogf, 'fasta')
		else:
			with open(fasta) as of:
				SeqIO.parse(of, 'fasta')
		return True
	except:
		return False


def is_genbank(gbk):
	"""
	Function to check in Genbank file is correctly formatted.
	"""
	try:
		assert (gbk.endswith('.gbk') or gbk.endswith('.gbff') or gbk.endswith('.gbk.gz') or gbk.endswith('.gbff.gz'))
		if gbk.endswith('.gz'):
			with gzip.open(gbk, 'rt') as ogf:
				SeqIO.parse(ogf, 'genbank')
		else:
			with open(gbk) as of:
				SeqIO.parse(of, 'genbank')
		return True
	except:
		return False

def checkValidGenBank(genbank_file):
	try:
		number_of_cds = 0
		lt_has_comma = False
		with open(genbank_file) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type == 'CDS':
						number_of_cds += 1
						lt = feature.qualifiers.get('locus_tag')[0]
						if ',' in lt:
							lt_has_comma = True

		if number_of_cds > 0 and not lt_has_comma:
			return True
		else:
			return False
	except:
		return False

def parseGenbankForCDSProteinsAndDNA(gbk_path, logObject):
	try:
		proteins = {}
		nucleotides = {}
		with open(gbk_path) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				full_sequence = str(rec.seq).upper()
				for feature in rec.features:
					if feature.type != 'CDS': continue
					lt = feature.qualifiers.get('locus_tag')[0]
					prot_seq = feature.qualifiers.get('translation')[0]
					all_coords = []
					if not 'join' in str(feature.location):
						start = min([int(x.strip('>').strip('<')) for x in
									 str(feature.location)[1:].split(']')[0].split(':')]) + 1
						end = max([int(x.strip('>').strip('<')) for x in
								   str(feature.location)[1:].split(']')[0].split(':')])
						direction = str(feature.location).split('(')[1].split(')')[0]
						all_coords.append([start, end, direction])
					else:
						all_starts = []
						all_ends = []
						all_directions = []
						for exon_coord in str(feature.location)[5:-1].split(', '):
							start = min([int(x.strip('>').strip('<')) for x in
										 exon_coord[1:].split(']')[0].split(':')]) + 1
							end = max([int(x.strip('>').strip('<')) for x in
									   exon_coord[1:].split(']')[0].split(':')])
							direction = exon_coord.split('(')[1].split(')')[0]
							all_starts.append(start);
							all_ends.append(end);
							all_directions.append(direction)
							all_coords.append([start, end, direction])
						assert (len(set(all_directions)) == 1)
					nucl_seq = ''
					for sc, ec, dc in sorted(all_coords, key=itemgetter(0), reverse=False):
						if ec >= len(full_sequence):
							nucl_seq += full_sequence[sc - 1:]
						else:
							nucl_seq += full_sequence[sc - 1:ec]
					if direction == '-':
						nucl_seq = str(Seq(nucl_seq).reverse_complement())
					proteins[lt] = prot_seq
					nucleotides[lt] = nucl_seq
		return([proteins, nucleotides])
	except Exception as e:
		sys.stderr.write('Issues with parsing the GenBank %s\n' % gbk_path)
		logObject.error('Issues with parsing the GenBank %s' % gbk_path)
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)


def parseVersionFromSetupPy():
	"""
	Parses version from setup.py program.
	"""
	setup_py_prog = main_dir + 'setup.py'
	version = 'NA'
	with open(setup_py_prog) as osppf:
		for line in osppf:
			line = line.strip()
			if line.startswith('version='):
				version = line.split('version=')[1][:-1]
	return version

def createLoggerObject(log_file):
	"""
	Function which creates logger object.
	:param log_file: path to log file.
	:return: logging logger object.
	"""
	logger = logging.getLogger('task_logger')
	logger.setLevel(logging.DEBUG)
	# create file handler which logs even debug messages
	fh = logging.FileHandler(log_file)
	fh.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M")
	fh.setFormatter(formatter)
	logger.addHandler(fh)
	return logger


def closeLoggerObject(logObject):
	"""
	Function which closes/terminates loggerObject.
	:param logObject: logging logger object to close
	"""
	handlers = logObject.handlers[:]
	for handler in handlers:
		handler.close()
		logObject.removeHandler(handler)


def logParameters(parameter_names, parameter_values):
	"""
	Function to log parameters of executable program to std.stderr
	"""
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		sys.stderr.write(pn + ': ' + str(pv) + '\n')


def logParametersToFile(parameter_file, parameter_names, parameter_values):
	"""
	Function to log parameters of executable program to text file.
	"""
	parameter_handle = open(parameter_file, 'w')
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		parameter_handle.write(pn + ': ' + str(pv) + '\n')
	parameter_handle.close()


def logParametersToObject(logObject, parameter_names, parameter_values):
	"""
	Function to log parameters of executable program to text file.
	"""
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		logObject.info(pn + ': ' + str(pv))


def is_numeric(x):
	try:
		x = float(x)
		return True
	except:
		return False


def castToNumeric(x):
	try:
		x = float(x)
		return (x)
	except:
		return float('nan')

""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

numeric_columns = set(['GCF Count', 'hg order index', 'hg consensus direction', 'median gene length',
					   'proportion of samples with hg', 'proportion of total populations with hg',
					   'hg median copy count', 'num of hg instances', 'samples with hg', 'ambiguous sites proporition',
					   'Tajimas D', 'proportion variable sites', 'proportion nondominant major allele',
					   'median beta rd',
					   'median dn ds', 'mad dn ds', 'populations with hg', 'proportion of total populations with hg',
					   'most significant Fisher exact pvalues presence absence', 'median Tajimas D per population',
					   'mad Tajimas D per population'])


def loadSampleToGCFIntoPandaDataFrame(gcf_listing_dir):
	import pandas as pd
	panda_df = None
	try:
		data = []
		data.append(['GCF', 'Sample', 'BGC Instances'])
		for f in os.listdir(gcf_listing_dir):
			gcf = f.split('.txt')[0]
			sample_counts = defaultdict(int)
			with open(gcf_listing_dir + f) as ogldf:
				for line in ogldf:
					line = line.strip()
					sample, bgc_path = line.split('\t')
					sample_counts[sample] += 1
			for s in sample_counts:
				data.append([gcf, s, sample_counts[s]])

		panda_dict = {}
		for ls in zip(*data):
			key = ' '.join(ls[0].split('_'))
			vals = ls[1:]
			panda_dict[key] = vals
		panda_df = pd.DataFrame(panda_dict)

	except Exception as e:
		raise RuntimeError(traceback.format_exc())
	return panda_df

def calculateSelectDistances(newick_file, selected_pairs):
	try:
		t = Tree(newick_file)
		leafs = set([])
		for node in t.traverse('postorder'):
			if node.is_leaf():
				leafs.add(node.name)
		pw_info = defaultdict(lambda: "nan")
		for i, n1 in enumerate(sorted(leafs)):
			for j, n2 in enumerate(sorted(leafs)):
				if i >= j: continue
				if n1 == n2: continue
				gn1 = n1.split('|')[0]
				gn2 = n2.split('|')[0]
				pw_key = tuple([gn1, gn2])
				pw_dist = t.get_distance(n1, n2)
				if pw_key in selected_pairs:
					if pw_key in pw_info and pw_dist < pw_info[pw_key]:
						pw_info[pw_key] = pw_dist
					else:
						pw_info[pw_key] = pw_dist
		return (pw_info)
	except Exception as e:
		sys.stderr.write('Issues with calculating pairwise distances for tree: %s.\n' % newick_file)
		sys.stderr.write(str(e) + '\n')
		raise RuntimeError(traceback.format_exc())
		sys.exit(1)

def computeCongruence(hg, gene_tree, gc_pw_info, selected_pairs, outf, logObject):
	try:
		hg_pw_info = calculateSelectDistances(gene_tree, selected_pairs)
		hg_pw_dists_filt = []
		gc_pw_dists_filt = []
		for pair in selected_pairs:
			if hg_pw_info[pair] != 'nan' and gc_pw_info[pair] != 'nan':
				hg_pw_dists_filt.append(hg_pw_info[pair])
				gc_pw_dists_filt.append(gc_pw_info[pair])

		congruence_stat = 'NA'
		if len(hg_pw_dists_filt) >= 3:
			stat, pval = stats.pearsonr(hg_pw_dists_filt, gc_pw_dists_filt)
			if stat != 'nan':
				congruence_stat = stat
		out_handle = open(outf, 'w')
		out_handle.write(str(congruence_stat) + '\n')
		out_handle.close()

	except Exception as e:
		sys.stderr.write('Issues with computing congruence of gene tree for homolog group %s to gene-cluster consensus tree.\n' % hg)
		logObject.error('Issues with computing congruence of gene tree for homolog group %s to gene-cluster consensus tree.' % hg)
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def checkCoreHomologGroupsExist(ortho_matrix_file):
	try:
		core_hgs = set([])
		with open(ortho_matrix_file) as omf:
			for i, line in enumerate(omf):
				if i == 0: continue
				line = line.rstrip('\n')
				ls = line.split('\t')
				hg = ls[0]
				sample_count = 0
				for lts in ls[1:]:
					if not lts.strip() == '':
						sample_count += 1
				if sample_count / float(len(ls[1:])) == 1.0:
					core_hgs.add(hg)
		assert(len(core_hgs) != 0)
		return True
	except:
		return False

def loadSamplesIntoPandaDataFrame(sample_annot_file, pop_spec_file=None):
	import pandas as pd
	panda_df = None
	try:
		if pop_spec_file != None:
			sample_pops = {}
			with open(pop_spec_file) as opsf:
				for line in opsf:
					line = line.strip()
					ls = line.split('\t')
					sample_pops[ls[0]] = ls[1]

		data = []
		if pop_spec_file != None:
			data.append(['Sample', 'Population/Clade'])
		else:
			data.append(['Sample'])
		with open(sample_annot_file) as saf:
			for line in saf:
				line = line.strip('\n')
				ls = line.split('\t')
				if pop_spec_file != None:
					pop = sample_pops[ls[0]]
					data.append([ls[0], pop])
				else:
					data.append([ls[0]])

		panda_dict = {}
		for ls in zip(*data):
			key = ' '.join(ls[0].split('_'))
			vals = ls[1:]
			panda_dict[key] = vals
		panda_df = pd.DataFrame(panda_dict)

	except Exception as e:
		raise RuntimeError(traceback.format_exc())
	return panda_df

def parseGbk(gbk_path, prefix, logObject):
	try:
		gc_gene_locations = {}
		with open(gbk_path) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type != 'CDS': continue
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
					dir = all_coords[0][2]
					for sc, ec, dc in sorted(all_coords, key=itemgetter(0), reverse=False):
						if sc < start:
							start = sc
						if ec > end:
							end = ec
					location = {'start': start, 'end': end, 'direction': dir}
					gc_gene_locations[prefix + '|' + lt] = location
		return gc_gene_locations
	except Exception as e:
		sys.stderr.write('Issue parsing GenBank %s\n' % gbk_path)
		logObject.error('Issue parsing GenBank %s' % gbk_path)
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def gatherAnnotationFromDictForHomoloGroup(hg, dict):
	try:
		annot_set_filt = set([x for x in dict[hg][0] if x.strip() != ''])
		assert(len(annot_set_filt) > 0)
		return('; '.join(annot_set_filt) + ' (' + str(dict[hg][1]) + ')')
	except:
		return('NA')

def gatherValueFromDictForHomologGroup(hg, dict):
	try:
		return (dict[hg])
	except:
		return ("NA")


def loadTableInPandaDataFrame(input_file, numeric_columns):
	import pandas as pd
	panda_df = None
	try:
		data = []
		with open(input_file) as oif:
			for line in oif:
				line = line.strip('\n')
				ls = line.split('\t')
				data.append(ls)

		panda_dict = {}
		for ls in zip(*data):
			key = ls[0]
			cast_vals = ls[1:]
			if key in numeric_columns:
				cast_vals = []
				for val in ls[1:]:
					cast_vals.append(castToNumeric(val))
			panda_dict[key] = cast_vals
		panda_df = pd.DataFrame(panda_dict)

	except Exception as e:
		raise RuntimeError(traceback.format_exc())
	return panda_df