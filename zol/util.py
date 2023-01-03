import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import logging
import subprocess
import statistics
from operator import itemgetter
from collections import defaultdict
import traceback
import multiprocessing
import copy
from scipy import stats
from ete3 import Tree
import itertools
import math
import numpy as np
import gzip
import pathlib
import operator

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


def calculateTajimasD(sequences):
	"""
	The code for this functionality was largely taken from Tom Whalley's Tajima's D implementation in Python and further
	modified/corrected based on Wikipedia's page for Tajima's D (Mathematical details).
	"""

	"""Calculate pi"""
	numseqs = len(sequences)
	divisor = math.comb(numseqs, 2)
	combos = itertools.combinations(sequences, 2)
	differences = 0
	for pair in combos:
		seqA = pair[0]
		seqB = pair[1]
		for p, a in enumerate(seqA):
			b = seqB[p]
			if a != b and a != '-' and b != '-':
				differences += 1
	pi = float(differences) / divisor

	"""Calculate s, number of segregation sites)."""
	# Assume if we're in here seqs have already been checked
	combos = itertools.combinations(sequences, 2)
	indexes = set([])
	for pair in combos:
		seqA = pair[0]
		seqB = pair[1]
		for idx, (i, j) in enumerate(zip(seqA, seqB)):
			if i != j and i != '-' and j != '-':
				indexes.add(idx)

	indexes = list(indexes)
	S = len(indexes)

	"""
	Now we have pi (pairwise differences) and s (number
	of segregating sites). This gives us 'little d', so
	now we need to divide it by sqrt of variance.
	"""
	l = len(sequences)

	# calculate D
	a1 = sum([(1.0 / float(i)) for i in range(1, l)])
	a2 = sum([(1.0 / (i ** 2)) for i in range(1, l)])

	b1 = float(l + 1) / (3 * (l - 1))
	b2 = float(2 * ((l ** 2) + l + 3)) / (9 * l * (l - 1))

	c1 = b1 - (1.0 / a1)
	c2 = b2 - (float(l + 2) / (a1 * l)) + (float(a2) / (a1 ** 2.0))

	e1 = float(c1) / a1
	e2 = float(c2) / ((a1 ** 2) + a2)
	if S >= 3:
		D = (float(pi - (float(S) / a1)) / math.sqrt((e1 * S) + ((e2 * S) * (S - 1))))
		return (D)
	else:
		return ("NA")


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


def loadTableInPandaDataFrame(input_file):
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
			key = ' '.join(ls[0].split('_'))
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


def loadCustomPopGeneTableInPandaDataFrame(input_file):
	import pandas as pd
	panda_df = None
	try:
		ignore_data_cats = {'hg_median_copy_count', 'proportion_of_total_populations_with_hg',
							'proportion_variable_sites', 'proportion_nondominant_major_allele', 'median_dn_ds',
							'mad_dn_ds', 'all_domains', 'most_significant_Fisher_exact_pvalues_presence_absence',
							'median_Tajimas_D_per_population', 'mad_Tajimas_D_per_population',
							'most_negative_population_Tajimas_D', 'most_positive_population_Tajimas_D',
							'population_entropy', 'median_fst_like_estimate'}
		data = []
		with open(input_file) as oif:
			for line in oif:
				line = line.strip('\n')
				ls = line.split('\t')
				data.append(ls)

		panda_dict = {}
		for ls in zip(*data):
			if ls[0] in ignore_data_cats: continue
			if ls[0] == 'gcf_annotation':
				updated_ans = []
				for ans in ls[1:]:
					upans = []
					for an in ans.split('; '):
						if not an.startswith('NA:'):
							upans.append(an)
					updated_ans.append('|'.join(upans))
				panda_dict[' '.join(ls[0].split('_'))] = updated_ans
			elif ls[0] == 'annotation':
				key = ' '.join(ls[0].split('_'))
				cleaned_annots = []
				for val in ls[1:]:
					ca = []
					for a in val.split('; '):
						if a != 'hypothetical protein':
							ca.append(a)
					if len(ca) == 0:
						cleaned_annots.append('hypothetical protein')
					else:
						cleaned_annots.append('; '.join(ca))
				panda_dict[key] = cleaned_annots
			else:
				key = ' '.join(ls[0].split('_'))
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


