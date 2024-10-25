import os
import sys
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import decimal
from Bio.Seq import Seq
import logging
import subprocess
from operator import itemgetter
from collections import defaultdict
import traceback
import numpy as np
import gzip
import copy
import itertools
import multiprocessing
import pickle
import resource
import pkg_resources  # part of setuptools
import shutil
from ete3 import Tree
import math
from zol import fai, zol
import statistics
import pyhmmer
import pandas as pd
from scipy import stats

version = pkg_resources.require("zol")[0].version

valid_alleles = set(['A', 'C', 'G', 'T'])

def memory_limit(mem):
	"""
	Description:
	Experimental function to limit memory.
	********************************************************************************************************************
	Parameters:
	- mem: The memory limit in GB.
	********************************************************************************************************************
	"""
	max_virtual_memory = mem*1000000000
	soft, hard = resource.getrlimit(resource.RLIMIT_AS)
	resource.setrlimit(resource.RLIMIT_AS, (max_virtual_memory, hard))
	print(resource.getrlimit(resource.RLIMIT_AS))


def runCmdViaSubprocess(cmd, logObject, check_files=[], check_directories=[], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL):
	"""
	Description:
	Function to run some command via subprocess.
	********************************************************************************************************************
	Parameters:
	- cmd: The command as a list.
	- logObject: A logging object.
	- check_files: Files to check the existence of assuming successful run of the command.
	- check_directories: Directories to check the existence of assuming successful run of the command.
	- stdout: Where to have subprocess direct standard output.
	- stderr: Where to have subprocess direct standard errorr.
	********************************************************************************************************************
	"""
	logObject.info('Running %s' % ' '.join(cmd))
	try:
		subprocess.call(' '.join(cmd), shell=True, stdout=stdout, stderr=stderr,
						executable='/bin/bash')
		for cf in check_files:
			assert (os.path.isfile(cf))
		for cd in check_directories:
			assert (os.path.isdir(cd))
		logObject.info('Successfully ran: %s' % ' '.join(cmd))
	except:
		logObject.error('Had an issue running: %s' % ' '.join(cmd))
		logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(cmd))


def cleanUpSampleName(original_name):
	"""
	Description:
	Function to clean up sample names for troublesome characters that makes unix based file creation tricky.
	********************************************************************************************************************
	Parameters:
	- original_name: The original name of the sample.
	********************************************************************************************************************
	"""
	return original_name.replace('#', '').replace('*', '_').replace(':', '_').replace(';', '_').replace(' ',
																										'_').replace(
		':', '_').replace('|', '_').replace('"', '_').replace("'", '_').replace("=", "_").replace('-', '_').replace('(',
																													'').replace(
		')', '').replace('/', '').replace('\\', '').replace('[', '').replace(']', '').replace(',', '')

def readInAnnotationFilesForExpandedSampleSet(expansion_listing_file, full_dir, logObject=None):
	"""
	Description:
	Function to read in GenBank paths from expansion listing file and load into dictionary with keys corresponding to
	sample IDs.
	********************************************************************************************************************
	Parameters:
	- expansion_listing_file: A tab-delimited file with two columns: (1) sample ID (2) GenBank file name.
	- full_dir: The path to where target genome GenBanks are stored.
	- logObject: A logging object.
	********************************************************************************************************************
	Returns:
	- sample_annotation_data: A dictionary of dictionaries with primary keys as sample names and the secondary key as
	                          "genbank" with final values being paths to the corresponding GenBank file for a sample
	                          target genome.
	********************************************************************************************************************
	"""
	sample_annotation_data = defaultdict(dict)
	try:
		with open(expansion_listing_file) as oalf:
			for line in oalf:
				line = line.strip()
				sample, genbank = line.split('\t')
				genbank = full_dir + genbank.split('/')[-1]
				sample = cleanUpSampleName(sample)
				try:
					assert(os.path.isfile(genbank) and os.path.isfile(genbank))
					sample_annotation_data[sample]['genbank'] = genbank
				except Exception as e:
					if logObject:
						logObject.warning('Ignoring sample %s, because at least one of two annotation files does not seem to exist.' % sample)
					else:
						sys.stderr.write('Ignoring sample %s, because at least one of two annotation files does not seem to exist.\n' % sample)
		assert (len(sample_annotation_data) >= 1)
		return (sample_annotation_data)
	except Exception as e:
		if logObject:
			logObject.error("Input file listing the location of annotation files for samples leads to incorrect paths or something else went wrong with processing of it. Exiting now ...")
			logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def createGenbank(full_genbank_file, new_genbank_file, scaffold, start_coord, end_coord):
	"""
	Description:
	Function to prune full genome-sized GenBank for only features in BGC of interest.
	********************************************************************************************************************
	Parameters:
	- full_genbank_file: GenBank file for full genome.
	- new_genbank_file: Path to gene cluster specific GenBank to be created.
	- scaffold: Scaffold identifier.
	- start_coord: Start coordinate.
	- end_coord: End coordinate.
	********************************************************************************************************************
	"""
	try:
		ngf_handle = open(new_genbank_file, 'w')
		pruned_coords = set(range(start_coord, end_coord + 1))
		full_genbank_index = SeqIO.index(full_genbank_file, 'genbank')
		rec = full_genbank_index[scaffold]
		original_seq = str(rec.seq)
		filtered_seq = ""
		start_coord = max(start_coord, 1)
		if end_coord >= len(original_seq):
			filtered_seq = original_seq[start_coord - 1:]
		else:
			filtered_seq = original_seq[start_coord - 1:end_coord]

		new_seq_object = Seq(filtered_seq)
		updated_rec = copy.deepcopy(rec)
		updated_rec.seq = new_seq_object

		updated_features = []
		#print(start_coord)
		#print(end_coord)
		#print('--------')
		for feature in rec.features:
			start = None
			end = None
			direction = None
			all_coords = []
			#print(str(feature.location))

			if not 'join' in str(feature.location) and not 'order' in str(feature.location):
				start = min([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
				end = max([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
				direction = str(feature.location).split('(')[1].split(')')[0]
				all_coords.append([start, end, direction])
			elif 'order' in str(feature.location):
				all_starts = []
				all_ends = []
				all_directions = []
				for exon_coord in str(feature.location)[6:-1].split(', '):
					start = min(
						[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
					end = max([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
					direction = exon_coord.split('(')[1].split(')')[0]
					all_starts.append(start)
					all_ends.append(end)
					all_directions.append(direction)
					all_coords.append([start, end, direction])
				start = min(all_starts)
				end = max(all_ends)
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
					all_coords.append([start, end, direction])
				start = min(all_starts)
				end = max(all_ends)

			feature_coords = set(range(start, end + 1))
			edgy_feat = 'False'
			if (start <= 2000) or (end + 1 >= len(original_seq)-2000):
				edgy_feat = 'True'
			part_of_cds_hanging = False
			if len(feature_coords.intersection(pruned_coords)) > 0:
				fls = []
				for sc, ec, dc in all_coords:
					exon_coord = set(range(sc, ec+1))
					if len(exon_coord.intersection(pruned_coords)) == 0: continue
					updated_start = sc - start_coord + 1
					updated_end = ec - start_coord + 1
					if ec > end_coord:
						# note overlapping genes in prokaryotes are possible so avoid proteins that overlap
						# with boundary proteins found by the HMM.
						if feature.type == 'CDS':
							part_of_cds_hanging = True
							continue
						else:
							updated_end = end_coord - start_coord + 1  # ; flag1 = True
					if sc < start_coord:
						if feature.type == 'CDS':
							part_of_cds_hanging = True
							continue
						else:
							updated_start = 1  # ; flag2 = True
					strand = 1
					if dc == '-':
						strand = -1
					fls.append(FeatureLocation(updated_start - 1, updated_end, strand=strand))
				if len(fls) > 0 and not part_of_cds_hanging:
					updated_location = fls[0]
					if len(fls) > 1:
						updated_location = sum(fls)
					feature.location = updated_location
					feature.qualifiers['near_scaffold_edge'] = edgy_feat
					updated_features.append(feature)
		updated_rec.features = updated_features
		SeqIO.write(updated_rec, ngf_handle, 'genbank')
		ngf_handle.close()
	except Exception as e:
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def multiProcess(input):
	"""
	Description:
	This is a generalizable function to be used with multiprocessing to parallelize list of commands. Inputs should
	correspond to space separated command (as list), with last item in list corresponding to a logging object handle for
	logging progress.
	********************************************************************************************************************
	Parameters:
	- input: A list corresponding to a command to run with the last item in the list corresponding to a logging object
	         for the function.
	********************************************************************************************************************
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
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def setupReadyDirectory(directories, delete_if_exist=False):
	"""
	Description:
	This is a generalizable function to create directories.
	********************************************************************************************************************
	Parameters:
	- dictionaries: A list of paths to directories to create or recreate (after removing).
	********************************************************************************************************************
	"""
	try:
		assert (type(directories) is list)
		for d in directories:
			if os.path.isdir(d):
				if delete_if_exist:
					os.system('rm -rfi %s' % d)
			else:
				os.system('mkdir %s' % d)
	except Exception as e:
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def is_fasta(fasta):
	"""
	Description:
	Function to validate if a file is actually a FASTA file.
	********************************************************************************************************************
	Parameters:
	- fasta: A file that should be in FASTA format.
	********************************************************************************************************************
	Returns:
	- True or False statement depending on whether file is in FASTA format.
	********************************************************************************************************************
	"""
	try:
		recs = 0
		if fasta.endswith('.gz'):
			with gzip.open(fasta, 'rt') as ogf:
				for rec in SeqIO.parse(ogf, 'fasta'):
					recs += 1
					break
		else:
			with open(fasta) as of:
				for rec in SeqIO.parse(of, 'fasta'):
					recs += 1
					break
		if recs > 0:
			return True
		else:
			return False
	except:
		return False

def is_genbank(gbk, check_for_cds=False):
	"""
	Description:
	Function to validate if a file is actually a GenBank file.
	********************************************************************************************************************
	Parameters:
	- gbk: A file that should be in GenBank format.
	- check_for_cds: Whether to also check that the GenBank contains CDS features.
	********************************************************************************************************************
	Returns:
	- True or False statement depending on whether file is in GenBank format.
	********************************************************************************************************************
	"""
	try:
		recs = 0
		cds_flag = False
		assert (gbk.endswith('.gbk') or gbk.endswith('.gbff') or gbk.endswith('.gbk.gz') or gbk.endswith('.gbff.gz') or gbk.endswith('.gb') or gbk.endswith('.gb.gz') or gbk.endswith('genbank') or gbk.endswith('.genbank.gz'))
		if gbk.endswith('.gz'):
			with gzip.open(gbk, 'rt') as ogf:
				for rec in SeqIO.parse(ogf, 'genbank'):
					if check_for_cds:
						for feature in rec.features:
							if feature.type == 'CDS':
								cds_flag = True
					if not check_for_cds or cds_flag:
						recs += 1
						break
		else:
			with open(gbk) as ogf:
				for rec in SeqIO.parse(ogf, 'genbank'):
					if check_for_cds:
						for feature in rec.features:
							if feature.type == 'CDS':
								cds_flag = True
					if not check_for_cds or cds_flag:
						recs += 1
						break
		if recs > 0:
			return True
		else:
			return False
	except:
		return False

def checkValidGenBank(gbk, quality_assessment=False, draft_assessment=False, use_either_lt_or_pi=False):
	"""
	Description:
	Function to check whether gene cluster GenBanks provided to zol as input meets the criteria requested by the user.
	********************************************************************************************************************
	Parameters:
	- gbk: The path to the GenBank file.
	- quality_assessment: Whether to check that most bases in the gene-cluster are non-ambiguous.
	- draft_assessment: Whether to check that gene-cluster does not lie on edge of the scaffold.
	- use_either_lt_or_pi: Whether protein_id is acceptable to use if locus_tag unavailable (currently for usage in fai,
	                       not zol).
	********************************************************************************************************************
	Returns:
	- True or False statement depending on whether file meets criteria for inclusion in zol analysis.
	********************************************************************************************************************
	"""
	try:
		number_of_cds = 0
		lt_has_comma = False
		lt_is_none = False
		seqs = ''
		recs = 0
		edgy_cds = False
		with open(gbk) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type == 'CDS':
						number_of_cds += 1
						lt = None
						pi = None
						try:
							lt = feature.qualifiers.get('locus_tag')[0]
						except:
							pass
						try:
							pi = feature.qualifiers.get('protein_id')[0]
						except:
							pass
						if use_either_lt_or_pi:
							if lt == None and pi != None:
								lt = pi
						if lt == None:
							lt_is_none = True
						if ',' in lt:
							lt_has_comma = True
						try:
							edgy_cds = feature.qualifiers.get('near_scaffold_edge')[0] == 'True'
							#print(edgy_cds)
						except:
							pass
				seqs += str(rec.seq)
				recs += 1
		prop_missing = sum([1 for bp in seqs if not bp in set(['A', 'C', 'G', 'T'])])/len(seqs)
		if number_of_cds > 0 and not lt_is_none and not lt_has_comma:
			if (quality_assessment and prop_missing >= 0.1) or (draft_assessment and (recs > 1 or edgy_cds)):
					return False
			else:
				return True
		else:
			return False
	except:
		return False

def convertGenbankToCDSProtsFasta(gbk, protein, logObject, use_either_lt_or_pi=False):
	"""
	Description:
	This function extracts protein sequences for CDS features from a GenBank.
	********************************************************************************************************************
	Parameters:
	- gbk: The path to the GenBank file.
	- protein: The path to the file to which to write the protein sequences to in FASTA format.
	- logObject: A logging object.
	- use_either_lt_or_pi: Whether protein_id is acceptable to use if locus_tag unavailable.
	********************************************************************************************************************
	"""
	try:
		prot_handle = open(protein, 'w')
		with open(gbk) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type != 'CDS': continue
					lt = None
					pi = None
					try:
						lt = feature.qualifiers.get('locus_tag')[0]
					except:
						pass
					try:
						pi = feature.qualifiers.get('protein_id')[0]
					except:
						pass
					if use_either_lt_or_pi:
						if lt == None and pi != None:
							lt = pi
					prot_seq = feature.qualifiers.get('translation')[0]
					prot_handle.write('>' + lt + '\n' + str(prot_seq) + '\n')
		prot_handle.close()
	except:
		sys.stderr.write('Difficulties in processing input GenBank and converting to protein fasta. Please make sure translation and locus_tag fields are available for GenBank!\n')
		logObject.error('Difficulties in processing input GenBank and converting to protein fasta. Please make sure translation and locus_tag fields are available for GenBank!')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def parseGenbankForCDSProteinsAndDNA(gbk, logObject, allow_edge_cds=True, feature_type='CDS'):
	"""
	Description:
	This function parses GenBank for CDS protein and nucleotide sequences.
	********************************************************************************************************************
	Parameters:
	- gbk: Path to the GenBank file.
	- logObject: A logging object.
	- allow_edge_cds: Whether to regard CDS features near scaffold edges.
	- feature_type: The focal type of feature. Default is CDS.
	********************************************************************************************************************
	Returns:
	- A list which can be expanded to the following dictionaries:
		- proteins: A dictionary mapping locus tag identifiers to protein sequences.
		- nucleotides: A dictionary mapping locus tag identifiers to nucleotide sequences.
		- upstream_regions: A dictionary mapping locus tag identifiers to upstream region sequences.
	********************************************************************************************************************
	"""
	try:
		proteins = {}
		nucleotides = {}
		upstream_regions = {}
		with open(gbk) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				full_sequence = str(rec.seq).upper()
				for feature in rec.features:
					if feature.type != feature_type: continue
					lt = feature.qualifiers.get('locus_tag')[0]

					all_coords = []
					all_starts = []
					all_ends = []
					if not 'join' in str(feature.location):
						start = min([int(x.strip('>').strip('<')) for x in
									 str(feature.location)[1:].split(']')[0].split(':')]) + 1
						end = max([int(x.strip('>').strip('<')) for x in
								   str(feature.location)[1:].split(']')[0].split(':')])
						direction = str(feature.location).split('(')[1].split(')')[0]
						all_coords.append([start, end, direction])
						all_starts.append(start)
						all_ends.append(end)
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


					max_ec = max(all_ends)
					min_sc = min(all_starts)
					nucl_seq = ''
					for sc, ec, dc in sorted(all_coords, key=itemgetter(0), reverse=False):
						if ec >= len(full_sequence):
							nucl_seq += full_sequence[sc - 1:]
						else:
							nucl_seq += full_sequence[sc - 1:ec]

					upstream_region = None
					if direction == '-':
						nucl_seq = str(Seq(nucl_seq).reverse_complement())
						if ec + 100 >= len(full_sequence):
							upstream_region = str(Seq(full_sequence[max_ec:max_ec+100]).reverse_complement())
					else:
						if sc - 100 >= 0:
							upstream_region = str(Seq(full_sequence[min_sc-101:min_sc-1]))

					final_prot_seq = None
					final_nucl_seq = None
					final_upstream_region = None
					edgy_cds = False

					try:
						final_upstream_region = feature.qualifiers.get('orf_upstream')[0]
					except:
						final_upstream_region = upstream_region

					try:
						final_nucl_seq = feature.qualifiers.get('open_reading_frame')[0]
					except:
						final_nucl_seq = nucl_seq
					try:
						final_prot_seq = feature.qualifiers.get('translation')[0]
					except:
						final_prot_seq = str(Seq(nucl_seq).translate())
					try:
						edgy_cds = feature.qualifiers.get('near_scaffold_edge')[0] == 'True'
					except:
						edgy_cds = False
					if allow_edge_cds or not edgy_cds:
						proteins[lt] = final_prot_seq
						nucleotides[lt] = final_nucl_seq
						if upstream_region != None:
							upstream_regions[lt] = final_upstream_region

		return([proteins, nucleotides, upstream_regions])
	except Exception as e:
		sys.stderr.write('Issues with parsing the GenBank %s\n' % gbk)
		logObject.error('Issues with parsing the GenBank %s' % gbk)
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
def getVersion():
	"""
	Description:
	Parses the version of the zol suite from the setup.py program.
	********************************************************************************************************************
	"""
	return(str(version))

def default_to_regular(d):
	"""
	Description:
	Function to convert defaultdict of defaultdict to dict of dict
	Taken from Martijn Pieters response in StackOverflow:
	https://stackoverflow.com/questions/26496831/how-to-convert-defaultdict-of-defaultdicts-of-defaultdicts-to-dict-of-dicts-o
	********************************************************************************************************************
	Parameters:
	- d: Defaultdict dictionary.
	********************************************************************************************************************
	Returns:
	- A regular old dictionary.
	******************************************************************************************************************
	"""
	if isinstance(d, defaultdict):
		d = {k: default_to_regular(v) for k, v in d.items()}
	return d

def createLoggerObject(log_file):
	"""
	Description:
	This function creates a logging object.
	********************************************************************************************************************
	Parameters:
	- log_file: Path to file to which to write logging.
	********************************************************************************************************************
	Returns:
	- logger: A logging object.
	********************************************************************************************************************
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
	Description:
	This function closes a logging object.
	********************************************************************************************************************
	Parameters:
	- logObject: A logging object.
	********************************************************************************************************************
	"""

	handlers = logObject.handlers[:]
	for handler in handlers:
		handler.close()
		logObject.removeHandler(handler)


def logParametersToFile(parameter_file, parameter_names, parameter_values):
	"""
	Description:
	This function serves to create a parameters input file for major programs, e.g. fai and zol.
	********************************************************************************************************************
	Parameters:
	- parameter_file: The path to the file where to write parameter information to. Will overwrite each time.
	- parameter_names: A list containing the parameter names.
	- parameter_values: A list in the same order as parameter_names which contains the respective arguments provided.
	********************************************************************************************************************
	"""
	parameter_handle = open(parameter_file, 'w')
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		parameter_handle.write(pn + ': ' + str(pv) + '\n')
	parameter_handle.close()

def is_integer(x):
	"""
	Description:
	This function checks whether the input variable corresponds to an integer.
	********************************************************************************************************************
	Parameters:
	- x: Input variable.
	********************************************************************************************************************
	Returns:
	- True or False statement depending on whether input variable is an integer.
	********************************************************************************************************************
	"""
	try:
		x = int(x)
		return True
	except:
		return False

def is_numeric(x):
	"""
	Description:
	This function checks whether the input variable is numeric.
	********************************************************************************************************************
	Parameters:
	- x: Input variable.
	********************************************************************************************************************
	Returns:
	- True or False statement depending on whether input variable is numeric.
	********************************************************************************************************************
	"""
	try:
		x = float(x)
		return True
	except:
		return False

def castToNumeric(x):
	"""
	Description:
	This function attempts to cast a variable into a float. A special exception is whether "< 3 segregating sites!" is
	the value of the variable, which will simply be retained as a string.
	********************************************************************************************************************
	Parameters:
	- x: Input variable.
	********************************************************************************************************************
	Returns:
	- A float casting of the variable's value if numeric or "nan" if not.
	********************************************************************************************************************
	"""
	try:
		if x == '< 3 segregating sites!':
			return(x)
		else:
			x = float(x)
			return (x)
	except:
		return float('nan')

def checkCoreHomologGroupsExist(ortho_matrix_file):
	"""
	Description:
	This function checks whether at least one core ortholog group exists within an ortholog group by sample matrix.
	********************************************************************************************************************
	Parameters:
	- orthogroup_matrix_file: The ortholog group vs sample matrix file, where cells correspond to locus tag identifiers.
	********************************************************************************************************************
	Returns:
	- True or False depending on whether a core ortholog group is found.
	********************************************************************************************************************
	"""

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

def processGenomesUsingMiniprot(reference_proteome, sample_genomes, additional_miniprot_outdir,
								additional_proteomes_directory, additional_genbanks_directory, logObject, threads=1,
								locus_tag_length=3):
	"""
	Description:
	This function oversees processing of input genomes to create proteome and GenBank files using miniprot.
	********************************************************************************************************************
	Parameters:
	- reference_proteome: The reference proteome (in FASTA format) to use to map CDS features onto target sample
	                      genomes.
	- sample_genomes: A dictionary mapping sample identifiers to the path of their genomes in FASTA format.
	- additional_miniprot_outdir: Workspace where miniprot (intermediate) results should be written to directly.
	- additional_proteomes_directory: Directory where final proteome files (in FASTA format) for target genomes will be
	                                  saved.
	- additional_genbanks_directory: Directory where final GenBank files for target genomes will be saved.
	- logObject: A logging object.
	- threads: The number of threads to use.
	- locus_tag_length: The length of the locus tags to generate.
	********************************************************************************************************************
	"""
	try:
		alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		possible_locustags = sorted(list([''.join(list(x)) for x in list(itertools.product(alphabet, repeat=locus_tag_length))]))

		miniprot_cmds = []
		for i, sample in enumerate(sample_genomes):
			sample_assembly = sample_genomes[sample]
			sample_locus_tag = ''.join(list(possible_locustags[i]))

			sample_mp_db = additional_miniprot_outdir + sample + '.mpi'
			sample_mp_gff = additional_miniprot_outdir + sample + '.gff'
			sample_mp_gbk = additional_genbanks_directory + sample + '.gbk'
			sample_mp_faa = additional_proteomes_directory + sample + '.faa'

			miniprot_index_cmd = ['miniprot', '-t1', '-d', sample_mp_db, sample_assembly]
			miniprot_run_cmd = ['miniprot', '--gff', '-t1', sample_mp_db, reference_proteome, '>', sample_mp_gff]
			miniprot_process_cmd = ['convertMiniprotGffToGbkAndProt.py', '-g', sample_mp_gff, '-f', sample_assembly,
									'-l', sample_locus_tag, '-og', sample_mp_gbk, '-op', sample_mp_faa]
			miniprot_cmds.append(miniprot_index_cmd + [';'] + miniprot_run_cmd + [';'] + miniprot_process_cmd + [logObject])

		p = multiprocessing.Pool(threads)
		p.map(multiProcess, miniprot_cmds)
		p.close()

		for sample in sample_genomes:
			try:
				sample_mp_gbk = additional_genbanks_directory + sample + '.gbk'
				sample_mp_faa = additional_proteomes_directory + sample + '.faa'
				assert (os.path.isfile(sample_mp_gbk) and os.path.isfile(sample_mp_faa))
			except:
				sys.stderr.write("Unable to validate successful genbank/predicted-proteome creation for sample %s" % sample)
				sys.stderr.write(traceback.format_exc())
				sys.exit(1)
	except Exception as e:
		logObject.error("Problem with creating commands for running miniprot or convertMiniprotGffToGbkAndProt.py. Exiting now ...")
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def processGenomesUsingProdigal(sample_genomes, prodigal_outdir, prodigal_proteomes, prodigal_genbanks, logObject,
								threads=1, locus_tag_length=3, gene_calling_method="pyrodigal", meta_mode=False,
								avoid_locus_tags=set([])):
	"""
	Description:
	This function oversees processing of input genomes to create proteome and GenBank files using p(y)rodigal.
	********************************************************************************************************************
	Parameters:
	- sample_genomes: A dictionary mapping sample identifiers to the path of their genomes in FASTA format.
	- prodigal_outdir: Workspace where prodigal (intermediate) results should be written to directly.
	- prodigal_proteomes: Directory where final proteome files (in FASTA format) for target genomes will be saved.
	- prodigal_genbanks_directory: Directory where final GenBank files for target genomes will be saved.
	- logObject: A logging object.
	- threads: The number of threads to use.
	- locus_tag_length: The length of the locus tags to generate.
	- gene_calling_method: Whether to use pyrodigal (default), prodigal, or prodigal-gv.
	- meta_mode: Whether to run pyrodigal/prodigal in metagenomics mode.
	- avoid_locus_tags: Whether to avoid using certain locus tags.
	********************************************************************************************************************
	"""
	try:
		alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		possible_locustags = sorted(list(
			set([''.join(list(x)) for x in list(itertools.product(alphabet, repeat=locus_tag_length))]).difference(
				avoid_locus_tags)))

		prodigal_cmds = []
		for i, sample in enumerate(sample_genomes):
			sample_assembly = sample_genomes[sample]
			sample_locus_tag = ''.join(list(possible_locustags[i]))

			prodigal_cmd = ['runProdigalAndMakeProperGenbank.py', '-i', sample_assembly, '-s', sample, '-gcm', gene_calling_method,
							'-l', sample_locus_tag, '-o', prodigal_outdir]
			if meta_mode:
				prodigal_cmd += ['-m']
			prodigal_cmds.append(prodigal_cmd + [logObject])

		p = multiprocessing.Pool(threads)
		p.map(multiProcess, prodigal_cmds)
		p.close()

		for sample in sample_genomes:
			try:
				assert (os.path.isfile(prodigal_outdir + sample + '.faa') and os.path.isfile(
					prodigal_outdir + sample + '.gbk'))
				os.system('mv %s %s' % (prodigal_outdir + sample + '.gbk', prodigal_genbanks))
				os.system('mv %s %s' % (prodigal_outdir + sample + '.faa', prodigal_proteomes))
			except:
				sys.stderr.write("Unable to validate successful genbank/predicted-proteome creation for sample %s\n" % sample)
				sys.stderr.write(traceback.format_exc())
				sys.exit(1)
	except Exception as e:
		logObject.error(
			"Problem with creating commands for running prodigal via script runProdigalAndMakeProperGenbank.py. Exiting now ...")
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
		
def processGenomesAsGenbanks(sample_genomes, proteomes_directory, genbanks_directory, gene_name_mapping_outdir,
							 logObject, threads=1, locus_tag_length=3, avoid_locus_tags=set([]),
							 rename_locus_tags=False):
	"""
	Description:
	This function oversees processing of input genomes as GenBanks with CDS features already available.
	********************************************************************************************************************
	Parameters:
	- sample_genomes: A dictionary mapping sample identifiers to the path of their genomes in GenBank format with CDS
	                  features available.
	- proteomes_directory: Directory where final proteome files (in FASTA format) for target genomes will be saved.
	- genbanks_directory: Directory where final GenBank files for target genomes will be saved.
	- gene_name_mapping_outdir: Directory where mapping files for original locus tags to new locus tags will be saved.
	- logObject: A logging object.
	- threads: The number of threads to use.
	- locus_tag_length: The length of the locus tags to generate.
	- avoid_locus_tags: Whether to avoid using certain locus tags.
	- rename_locus_tags: Whether to rename locus tags.
	********************************************************************************************************************
	Returns:
	- sample_genomes_updated: Dictionary mapping sample names to paths of final/processed sample GenBanks.
	********************************************************************************************************************
	"""

	sample_genomes_updated = {}
	process_cmds = []
	try:
		alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		possible_locustags = sorted(list(
			set([''.join(list(x)) for x in list(itertools.product(alphabet, repeat=locus_tag_length))]).difference(
				avoid_locus_tags)))
		lacking_cds_gbks = set([])

		for i, sample in enumerate(sample_genomes):
			sample_locus_tag = possible_locustags[i]
			sample_genbank = sample_genomes[sample]
			process_cmd = ['processNCBIGenBank.py', '-i', sample_genbank, '-s', sample, 
						   '-g', genbanks_directory, '-p', proteomes_directory, '-n', gene_name_mapping_outdir]
			if rename_locus_tags:
				process_cmd += ['-l', sample_locus_tag]
			process_cmds.append(process_cmd + [logObject])

		p = multiprocessing.Pool(threads)
		p.map(multiProcess, process_cmds)
		p.close()

		for sample in sample_genomes:
			if sample in lacking_cds_gbks:
				continue
			try:
				assert (os.path.isfile(proteomes_directory + sample + '.faa') and
						os.path.isfile(genbanks_directory + sample + '.gbk') and
						os.path.isfile(gene_name_mapping_outdir + sample + '.txt'))
				sample_genomes_updated[sample] = genbanks_directory + sample + '.gbk'
			except:
				sys.stderr.write("Unable to validate successful genbank/predicted-proteome creation for sample %s" % sample)
				sys.stderr.write(traceback.format_exc())
				sys.exit(1)
	except Exception as e:
		logObject.error("Problem with processing existing Genbanks to (re)create genbanks/proteomes. Exiting now ...")
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return sample_genomes_updated

def determineGenomeFormat(inputs):
	"""
	Description:
	This function determines whether a target sample genome is provided in GenBank or FASTA format.
	********************************************************************************************************************
	Parameters:
	- input a list which can be expanded to the following:
		- sample: The sample identifier / name.
		- genome_File: The path to the genome file.
		- format_assess_dir: The directory where the genome type information for the sample will be written.
		- logObject: A logging object.
	********************************************************************************************************************
	"""

	sample, genome_file, format_assess_dir, logObject = inputs
	try:
		gtype = 'unknown'
		if is_fasta(genome_file):
			gtype = 'fasta'
		if is_genbank(genome_file, check_for_cds=True):
			if gtype == 'fasta':
				gtype = 'unknown'
			else:
				gtype = 'genbank'
		sample_res_handle = open(format_assess_dir + sample + '.txt', 'w')
		sample_res_handle.write(sample + '\t' + str(gtype) + '\n')
		sample_res_handle.close()
	except:
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def parseSampleGenomes(genome_listing_file, format_assess_dir, format_predictions_file, logObject, threads=1):
	"""
	Description:
	This function parses the input sample target genomes and determines whether they are all provided in the same format
	and whether everything aligns with expectations.
	********************************************************************************************************************
	Parameters:
	- genome_listing_file: A tab separated file with two columns: (1) sample name, (2) path to genome file.
	- format_assess_dir: The directory/workspace where genome format information will be saved.
	- format_predictions_file: The file where to concatenate genome format information.
	- logObject: A logging object.
	- threads: The number of threads to use.
	********************************************************************************************************************
	Returns:
	- sample_genomes: A dictionary which maps sample names to genome file paths (note, unknown format files will be
	                  dropped).
	- format_prediction: The format prediction for genome files.
	********************************************************************************************************************
	"""
	try:
		sample_genomes = {}
		assess_inputs = []
		with open(genome_listing_file) as oglf:
			for line in oglf:
				line = line.strip()
				ls = line.split('\t')
				sample, genome_file = ls
				assess_inputs.append([sample, genome_file, format_assess_dir, logObject])
				try:
					assert (os.path.isfile(genome_file))
				except:
					logObject.warning(
						"Problem with finding genome file %s for sample %s, skipping" % (genome_file, sample))
					continue
				if sample in sample_genomes:
					logObject.warning('Skipping genome %s for sample %s because a genome file was already provided for this sample' % (genome_file, sample))
					continue
				sample_genomes[sample] = genome_file

		p = multiprocessing.Pool(threads)
		p.map(determineGenomeFormat, assess_inputs)
		p.close()

		os.system('find %s -maxdepth 1 -type f | xargs cat >> %s' % (format_assess_dir, format_predictions_file))

		format_prediction = 'mixed'
		gtypes = set([])
		with open(format_predictions_file) as ofpf:
			for line in ofpf:
				line = line.strip()
				sample, gtype = line.split('\t')
				if gtype == 'unknown':
					sys.stderr.write('unsure about format for genome %s for sample %s, skipping inclusion...\n' % (sample_genomes[sample], sample))
					logObject.warning('unsure about format for genome %s for sample %s, skipping inclusion...' % (sample_genomes[sample], sample))
					del sample_genomes[sample]
				else:
					gtypes.add(gtype)

		if len(gtypes) == 1:
			format_prediction = list(gtypes)[0]

		return ([sample_genomes, format_prediction])

	except Exception as e:
		logObject.error("Problem with creating commands for running Prodigal. Exiting now ...")
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def filterRecordsNearScaffoldEdge(gbk, filt_genbank_file, logObject, quality_assessment=False):
	"""
	Description:
	This function filters specific records in a GenBank if they are near a scaffold edge.
	********************************************************************************************************************
	Parameters:
	- gbk: The GenBank file.
	- filt_genbank_file: The filtered GenBank file.
	- logObject: A logging object.
	- quality_assessment: Whether to perform quality assessment and drop the GenBank if >10% of nucleotides are
	                      ambiguous.
	********************************************************************************************************************
	"""
	try:
		number_of_cds = 0
		seqs = ""
		recs = 0
		recs_with_edgy_cds = set([])
		with open(gbk) as ogbk:
			for rec_it, rec in enumerate(SeqIO.parse(ogbk, 'genbank')):
				edgy_cds = False
				for feature in rec.features:
					if feature.type == 'CDS':
						number_of_cds += 1
						try:
							if feature.qualifiers.get('near_scaffold_edge')[0] == 'True':
								recs_with_edgy_cds.add(rec_it)
								edgy_cds = True
						except:
							pass
				if not edgy_cds:
					seqs += str(rec.seq)
				recs += 1

		if len(seqs) == 0: 
			return
		prop_missing = sum([1 for bp in seqs if not bp in set(['A', 'C', 'G', 'T'])]) / len(seqs)
		recs_without_edgy_cds = recs-len(recs_with_edgy_cds)
		if number_of_cds > 0 and (prop_missing <= 0.1 or not quality_assessment) and recs_without_edgy_cds > 0:
			out_handle = open(filt_genbank_file, 'w')
			with open(gbk) as ogbk:
				for rec_it, rec in enumerate(SeqIO.parse(ogbk, 'genbank')):
					if rec_it in recs_with_edgy_cds: continue
					SeqIO.write(rec, out_handle, 'genbank')
			out_handle.close()

	except Exception as e:
		sys.stderr.write('Issue parsing GenBank %s and CDS locus tag renaming.\n' % gbk)
		logObject.error('Issue parsing GenBank %s and CDS locus tag renaming.' % gbk)
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def renameCDSLocusTag(gbk, lt, rn_genbank_file, logObject, quality_assessment=False, draft_assessment=False):
	"""
	Description:
	This function renames or creates locus tags in a GenBank.
	********************************************************************************************************************
	Parameters:
	- gbk: The GenBank file.
	- lt: The new locus tag prefix.
	- rn_genbank_file: The path to the GenBank file to be created with the new locus tag names.
	- logObject: A logging object.
	- quality_assessment: Whether to perform quality assessment and drop the GenBank if >10% of nucleotides are
	                      ambiguous.
	- draft_assessment: Whether to perform draft quality assessment and drop the GenBank if there are no records not
	                    nearby scaffold edges.
	********************************************************************************************************************
	"""
	try:
		number_of_cds = 0
		seqs = ""
		recs = 0
		recs_with_edgy_cds = set([])
		with open(gbk) as ogbk:
			for rec_it, rec in enumerate(SeqIO.parse(ogbk, 'genbank')):
				edgy_cds = False
				for feature in rec.features:
					if feature.type == 'CDS':
						number_of_cds += 1
						try:
							if feature.qualifiers.get('near_scaffold_edge')[0] == 'True':
								recs_with_edgy_cds.add(rec_it)
								edgy_cds = True
						except:
							pass
				if not edgy_cds:
					seqs += str(rec.seq)
				recs += 1
		if len(seqs) == 0:
			return
		prop_missing = sum([1 for bp in seqs if not bp in set(['A', 'C', 'G', 'T'])]) / len(seqs)
		recs_without_edgy_cds = recs-len(recs_with_edgy_cds)
		if number_of_cds > 0 and (prop_missing <= 0.1 or not quality_assessment) and (recs_without_edgy_cds > 0 or not draft_assessment):
			out_handle = open(rn_genbank_file, 'w')
			locus_tag_iterator = 1
			with open(gbk) as ogbk:
				for rec_it, rec in enumerate(SeqIO.parse(ogbk, 'genbank')):
					if rec_it in recs_with_edgy_cds and draft_assessment: continue
					for feature in rec.features:
						if feature.type != 'CDS': continue
						new_locus_tag = lt + '_'
						if locus_tag_iterator < 10:
							new_locus_tag += '00000' + str(locus_tag_iterator)
						elif locus_tag_iterator < 100:
							new_locus_tag += '0000' + str(locus_tag_iterator)
						elif locus_tag_iterator < 1000:
							new_locus_tag += '000' + str(locus_tag_iterator)
						elif locus_tag_iterator < 10000:
							new_locus_tag += '00' + str(locus_tag_iterator)
						elif locus_tag_iterator < 100000:
							new_locus_tag += '0' + str(locus_tag_iterator)
						else:
							new_locus_tag += str(locus_tag_iterator)
						feature.qualifiers['locus_tag'] = new_locus_tag
						locus_tag_iterator += 1
					SeqIO.write(rec, out_handle, 'genbank')
			out_handle.close()
	except Exception as e:
		sys.stderr.write('Issue parsing GenBank %s and CDS locus tag renaming.\n' % gbk)
		logObject.error('Issue parsing GenBank %s and CDS locus tag renaming.' % gbk)
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def parseGbk(gbk, prefix, logObject, use_either_lt_or_pi=False, feature_type='CDS'):
	"""
	Description:
	This function parses CDS coordinate information from a GenBank.
	********************************************************************************************************************
	Parameters:
	- gbk: The GenBank file.
	- prefix: The prefix to append to locus tags (often gene cluster name) to make them use.
	- logObject: A logging object.
	- use_either_lt_or_pi: Use protein_id qualifier if locus_tag is unavailable for CDS feature.
	********************************************************************************************************************
	Returns:
	- gc_gene_locations: A dictionary for CDS locations where keys correspond to "prefix|locus_tag" and the values are
	                     another dictionary with the keys scaffold, start position, end position.
	********************************************************************************************************************
	"""
	try:
		gc_gene_locations = {}
		with open(gbk) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type != feature_type: continue
					lt = None
					pi = None
					try:
						lt = feature.qualifiers.get('locus_tag')[0]
					except:
						pass
					try:
						pi = feature.qualifiers.get('protein_id')[0]
					except:
						pass
					if use_either_lt_or_pi:
						if lt == None and pi != None:
							lt = pi

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
					dir = all_coords[0][2]
					for sc, ec, dc in sorted(all_coords, key=itemgetter(0), reverse=False):
						if sc < start:
							start = sc
						if ec > end:
							end = ec
					location = {'scaffold': rec.id, 'start': start, 'end': end, 'direction': dir}
					gc_gene_locations[prefix + '|' + lt] = location
		return gc_gene_locations
	except Exception as e:
		sys.stderr.write('Issue parsing GenBank %s\n' % gbk)
		logObject.error('Issue parsing GenBank %s' % gbk)
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def determinePossibleLTs():
	"""
	Description:
	This function creates a sorted list of possible locus tag prefices of length 4 each.
	********************************************************************************************************************
	Returns:
	- possible_locustags: A sorted list of possible locus tag prefices of length 4 each.
	********************************************************************************************************************
	"""
	alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
	possible_locustags = sorted(list([''.join(list(lt)) for lt in itertools.product(alphabet, repeat=4)]))
	return possible_locustags

def gatherAnnotationFromDictForHomoloGroup(hg, db, annot_dict):
	"""
	Description:
	This function formats the annotation information for the final zol TSV and XLSX for an ortholog group for a
	particular database.
	********************************************************************************************************************
	Parameters:
	- hg: The ortholog group identifier.
	- db: The database identifier.
	- annot_dict: The dictionary with annotation information.
	********************************************************************************************************************
	Returns:
	- A string with the annotation and the E-value in parentheses or "NA" if an annotation is not available.
	********************************************************************************************************************
	"""
	try:
		assert(db in annot_dict)
		annot_set_filt = set([x for x in annot_dict[db][hg][0] if x.strip() != ''])
		assert(len(annot_set_filt) > 0)
		return('; '.join(annot_set_filt) + ' (' + str(max(annot_dict[db][hg][1])) + ')')
	except:
		return('NA')

def gatherValueFromDictForHomologGroup(hg, info_dict):
	"""
	Description:
	This function formats/gathers information for an ortholog group from a dictionary.
	********************************************************************************************************************
	Parameters:
	- hg: The ortholog group identifier.
	- info_dict: The dictionary with information to be extracted for the ortholog group.
	********************************************************************************************************************
	Returns:
	- A string with the information for the ortholog group or "NA" if information is not available.
	********************************************************************************************************************
	"""
	try:
		return (info_dict[hg])
	except:
		return ("NA")

def loadTableInPandaDataFrame(input_file, numeric_columns):
	"""
	Description:
	This function formats reads a TSV file and stores it as a pandas dataframe.
	********************************************************************************************************************
	Parameters:
	- input_file: The input TSV file, with first row corresponding to the header.
	- numeric_columns: Set of column names which should have numeric data.
	********************************************************************************************************************
	Returns:
	- panda_df: A pandas DataFrame object reprsentation of the input TSV file.
	********************************************************************************************************************
	"""
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
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return panda_df

def chunks(lst, n):
	"""
	Description:
	Function to yield successive n-sized chunks from lst.
    Solution taken from: https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
	********************************************************************************************************************
	Parameters:
	- lst: A list.
	- n: The chunk size.
	********************************************************************************************************************
	Yields:
	- chunks of size n.
	********************************************************************************************************************
    """
	for i in range(0, len(lst), n):
		yield lst[i:i + n]

def parseGenbankAndFindBoundaryGenes(inputs):
	"""
	Description:
	Function to parse a full genome GenBank and writes a dictionary of genes per scaffold, gene to scaffold, and a
	set of CDS features which lie on the boundary of scaffolds (within 2000 bp).
	********************************************************************************************************************
	Parameters:
	- input: A list which can be expanded to the following:
		- sample: The sample / target genome identifier.
		- sample_genbank: The sample / target genome GenBank file.
		- pkl_result_file: The path to the pickle file to write with the sample / target genome information mentioned
		                   in the description.
	********************************************************************************************************************
	"""

	distance_to_scaffold_boundary = 2000
	gene_location = {}
	scaffold_genes = defaultdict(set)
	boundary_genes = set([])
	gene_id_to_order = defaultdict(dict)
	gene_order_to_id = defaultdict(dict)

	sample, sample_genbank, pkl_result_file = inputs
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
			elif 'order' in str(feature.location):
				all_starts = []
				all_ends = []
				all_directions = []
				for exon_coord in str(feature.location)[6:-1].split(', '):
					ec_start = min(
						[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
					ec_end = max(
						[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
					ec_direction = exon_coord.split('(')[1].split(')')[0]
					all_starts.append(ec_start)
					all_ends.append(ec_end)
					all_directions.append(ec_direction)
				assert (len(set(all_directions)) == 1)
				start = min(all_starts)
				end = max(all_ends)
				direction = all_directions[0]
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
				assert (len(set(all_directions)) == 1)
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


	sample_data = [gene_location, dict(scaffold_genes), boundary_genes, dict(gene_id_to_order), dict(gene_order_to_id)]

	with open(pkl_result_file, 'wb') as pickle_file:
		pickle.dump(sample_data, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)

def convertGenomeGenBankToFasta(inputs):
	input_gbk, output_fasta = inputs
	output_fasta_handle = open(output_fasta, 'w')
	with open(input_gbk) as oig:
		for rec in SeqIO.parse(oig, 'genbank'):
			output_fasta_handle.write('>' + rec.id + '\n' + str(rec.seq) + '\n')
	output_fasta_handle.close()

def createNJTree(additional_genbanks_directory, species_tree, workspace_dir, logObject, threads=1):
	"""
	Description:
	Function to create a species tree using ANI estimates + a neighbor joining approach.
	********************************************************************************************************************
	Parameters:
	- additional_genbanks_directory: Directory with full genomes in GenBank format.
	- workspace_dir: Workspace directory.
	- logObject: A logging object.
	- threads: The number of threads to use.
	********************************************************************************************************************
	"""

	try:
		tmp_genome_fasta_dir = workspace_dir + 'Genome_FASTAs/'
		setupReadyDirectory([tmp_genome_fasta_dir])

		all_genomes_listing_file = workspace_dir + 'Genome_FASTAs_Listing.txt'
		all_genomes_listing_handle = open(all_genomes_listing_file, 'w')
		conversion_inputs = []
		all_samples = set([])
		for f in os.listdir(additional_genbanks_directory):
			if not f.endswith('.gbk'): continue
			input_gbk = additional_genbanks_directory + f
			sample = '.gbk'.join(f.split('.gbk')[:-1])
			output_fasta = tmp_genome_fasta_dir + sample + '.fasta'
			all_samples.add(sample)
			conversion_inputs.append([input_gbk, output_fasta])
			all_genomes_listing_handle.write(output_fasta + '\n')
		all_genomes_listing_handle.close()

		# parallelize conversion
		p = multiprocessing.Pool(threads)
		p.map(convertGenomeGenBankToFasta, conversion_inputs)
		p.close()

		# run skani triangle
		skani_result_file = workspace_dir + 'Skani_Triangle_Edge_Output.txt'
		skani_triangle_cmd = ['skani', 'triangle', '-E', '-l', all_genomes_listing_file,'-t', str(threads), '-o', skani_result_file]
		try:
			subprocess.call(' '.join(skani_triangle_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
			assert (os.path.isfile(skani_result_file))
			logObject.info('Successfully ran: %s' % ' '.join(skani_triangle_cmd))
		except Exception as e:
			logObject.error('Had an issue with running skani: %s' % ' '.join(skani_triangle_cmd))
			sys.stderr.write('Had an issue with running skani: %s\n' % ' '.join(skani_triangle_cmd))
			logObject.error(e)
			sys.stderr.write(traceback.format_exc())
			sys.exit(1)

		shutil.rmtree(tmp_genome_fasta_dir)

		dist_matrix_file = workspace_dir + 'Skani_Based_Distance_Matrix.txt'
		dist_matrix_handle = open(dist_matrix_file, 'w')
		stos_dists = defaultdict(lambda: defaultdict(lambda: 1.0))
		with open(skani_result_file) as osrf:
			for i, line in enumerate(osrf):
				if i == 0: continue
				line = line.strip()
				ls = line.split('\t')
				s1 = '.fasta'.join(ls[0].split('/')[-1].split('.fasta')[:-1])
				s2 = '.fasta'.join(ls[1].split('/')[-1].split('.fasta')[:-1])
				if s1 != s2:
					ani = float(ls[2])/100.0
					dist_ani = 1.0 - ani
					stos_dists[s1][s2] = dist_ani
					stos_dists[s2][s1] = dist_ani
				else:
					stos_dists[s1][s2] = 0.0

		dist_matrix_handle.write('sample\t' + '\t'.join(sorted(all_samples)) + '\n')
		for s1 in sorted(all_samples):
			printlist = [s1]
			for s2 in sorted(all_samples):
				printlist.append(str(stos_dists[s1][s2]))
			dist_matrix_handle.write('\t'.join(printlist) + '\n')
		dist_matrix_handle.close()

		rscript_path = workspace_dir + 'generateNjTree.R'
		unrooted_tree_file = workspace_dir + 'Unrooted_Species_Tree.nwk'
		generateNjTree(rscript_path, dist_matrix_file, unrooted_tree_file, logObject)
		nj_tree_cmd = ['Rscript', rscript_path]
		try:
			subprocess.call(' '.join(nj_tree_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(unrooted_tree_file))
			logObject.info('Successfully ran: %s' % ' '.join(nj_tree_cmd))
		except Exception as e:
			logObject.error('Had an issue with generating neigbhor-joining tree: %s' % ' '.join(nj_tree_cmd))
			sys.stderr.write('Had an issue with generating neighbor-joining tree: %s\n' % ' '.join(nj_tree_cmd))
			logObject.error(e)
			sys.stderr.write(traceback.format_exc())
			sys.exit(1)

		# Midpoint the tree using ete3
		t = Tree(unrooted_tree_file)
		R = t.get_midpoint_outgroup()
		t.set_outgroup(R)
		t.write(outfile=species_tree,format=1)

	except Exception as e:
		sys.stderr.write('Issues with creating species tree.\n')
		sys.stderr.write(traceback.format_exc())
		logObject.error('Issues with creating species tree.')
		logObject.error(traceback.format_exc())		
		sys.exit(1)

def determineColumnNameBasedOnIndex(index):
	"""
	Function to determine spreadsheet column name for a given index
	"""
	# offset at 0 
	num_to_char = {}
	alphabet = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
	alphabet_first_spot = [''] + alphabet
	for i, c in enumerate(alphabet):
		num_to_char[i] = c
	level = math.floor(index / 26)
	remainder = index % 26
	columname = alphabet_first_spot[level]
	columname += alphabet[remainder]
	return columname

def cleanUp(clean_up_dirs_and_files, logObject):
	"""
	Basic function to clean-up disk heavy files/directories that are intermediate.
	"""
	try:
		for df in clean_up_dirs_and_files:
			if os.path.isfile(df):
				logObject.warning('Deleting the file %s' % df)
				os.system('rm -fi %s' % df)
			elif os.path.isdir(df):
				logObject.warning('Deleting the file %s' % df)
				shutil.rmtree(df)
			else:
				logObject.error('Couldn\'t find %s to delete!' % df)
	except:
		sys.stderr.write('Issues with cleaning up files/directories.\n')
		sys.stderr.write(traceback.format_exc())
		logObject.error('Issues with cleaning up files/directories.\n')
		logObject.error(traceback.format_exc())
		sys.exit(1)
		
def diamondBlastAndGetBestHits(cluster_id, query_protein_fasta, key_protein_fasta, target_genomes_db, workspace_dir, logObject,
							   identity_cutoff=40.0, coverage_cutoff=70.0, evalue_cutoff=1e-3, blastp_mode="very-sensitive", 
							   prop_key_prots_needed=0.0, threads=1):
	"""
	Description:

	********************************************************************************************************************
	Parameters:
	- cluster_id: Identifier of the query gene cluster (e.g. single BGC or phage) in consideration.
	- query_protein_fasta: FASTA file with all protein sequences for the query gene cluster.
	- key_protein_fasta: FASTA file with all the key (subset) protein sequences for the query gene cluster.
	- target_genomes_db: prepTG results directory of target genomes.
	- workspace_dir: Workspace directory.
	- logObject: A logging object.
	- identity_cutoff: Minimum percent identity for homologs to query gene cluster proteins.
	- coverage_cutoff: Minimum query coverage for homologs to query gene cluster proteins.
	- evalue_cutoff: Maximum E-value for homologs to query gene cluster proteins.
	- blastp_mode: Sensitivity mode for DIAMOND blastp.
	- prop_key_prots_needed: The proportion of key proteins needed for the query gene cluster to be deemed present in a 
	                         target genome.
	- threads: The number of threads to use.
	********************************************************************************************************************
	"""
	try:
		genome_wide_tsv_result_file = workspace_dir + 'total_gcs.tsv'
		target_genome_dmnd_db = target_genomes_db + 'Target_Genomes_DB.dmnd'

		# run diamond (use's fai function: runDiamondBlastp)
		diamond_results_file = workspace_dir + 'DIAMOND_Results.txt'
		fai.runDiamondBlastp(target_genome_dmnd_db, query_protein_fasta, workspace_dir, logObject,
					 diamond_sensitivity=blastp_mode, evalue_cutoff=evalue_cutoff, compute_query_coverage=True, 
					 threads=threads)

		all_hgs = set([])
		key_hgs = set([])

		with open(query_protein_fasta) as oqpf:
			for rec in SeqIO.parse(oqpf, 'fasta'):
				all_hgs.add(rec.id)
		
		with open(key_protein_fasta) as okpf:
			for rec in SeqIO.parse(okpf, 'fasta'):
				key_hgs.add(rec.id)

		# parse results (based closely on fai function: processDiamondBlastp)
		best_hit_per_lt = defaultdict(lambda: defaultdict(lambda: [0.0, [], [], []]))
		with open(diamond_results_file) as orf:
			for line in orf:
				line = line.strip()
				ls = line.split()
				hg = ls[0]
				sample = ls[1].split('|')[0]
				lt = ls[1].split('|')[1]
				identity = float(ls[2])
				qcovhsp = float(ls[7])
				if qcovhsp < coverage_cutoff or identity < identity_cutoff: continue
				bitscore = float(ls[4])
				qlen = float(ls[5])
				slen = float(ls[6])
				sql_ratio = float(slen)/float(qlen)
				if bitscore > best_hit_per_lt[sample][lt][0]:
					best_hit_per_lt[sample][lt][0] = bitscore
					best_hit_per_lt[sample][lt][1] = [hg]
					best_hit_per_lt[sample][lt][2] = [identity]
					best_hit_per_lt[sample][lt][3] = [sql_ratio]
				elif bitscore == best_hit_per_lt[sample][lt][0]:
					best_hit_per_lt[sample][lt][1].append(hg)
					best_hit_per_lt[sample][lt][2].append(identity)
					best_hit_per_lt[sample][lt][3].append(sql_ratio)

		sample_best_hit_per_hg = defaultdict(lambda: defaultdict(list))
		for sample in best_hit_per_lt:
			for lt in best_hit_per_lt[sample]:
				for hgi, hg in enumerate(best_hit_per_lt[sample][lt][1]):
					sample_best_hit_per_hg[sample][hg].append([best_hit_per_lt[sample][lt][0],
												               best_hit_per_lt[sample][lt][2][hgi], 
															   best_hit_per_lt[sample][lt][3][hgi], lt])

		outf_handle = open(genome_wide_tsv_result_file, 'w')
		outf_handle.write('genome\tblank1\taai\tblank2\tshared_gene_prop\n')
		for sample in sample_best_hit_per_hg:
			top_identities = []
			hgs_found = set([])
			key_hgs_found = set([])
			for hg in sample_best_hit_per_hg[sample]:
				for i, hginfo in enumerate(sorted(sample_best_hit_per_hg[sample][hg], key=itemgetter(0,1,2), reverse=False)):
					if i == 0:
						top_identities.append(hginfo[1])
						hgs_found.add(hginfo[3])
						if hg in key_hgs:
							key_hgs_found.add(hg)

			key_hgs_found_prop = float(len(key_hgs_found))/len(key_hgs)
			if key_hgs_found_prop < prop_key_prots_needed: continue
			aai = statistics.mean(top_identities)
			sgp = float(len(hgs_found))/len(all_hgs)
			printrow = [sample, 'NA', str(aai), 'NA', str(sgp)]
			outf_handle.write('\t'.join(printrow) + '\n')
		outf_handle.close()
		
	except:
		sys.stderr.write('Issues with running simple BLASTp for abon, atpoc, or apos.\n')
		sys.stderr.write(traceback.format_exc())
		logObject.error('Issues with running simple BLASTp for abon, atpoc, or apos.\n')
		logObject.error(traceback.format_exc())
		sys.exit(1)

def diamondBlast(inputs):
	og_prot_file, og_prot_dmnd_db, og_blast_file, logObject = inputs

	makedb_cmd = ['diamond', 'makedb', '--ignore-warnings', '--in', og_prot_file, '-d', og_prot_dmnd_db]
	search_cmd = ['diamond', 'blastp', '--ignore-warnings', '-p', '1', '-d', og_prot_dmnd_db, '-q', og_prot_file, '-o', og_blast_file]
	try:
		subprocess.call(' '.join(makedb_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		assert (os.path.isfile(og_prot_dmnd_db))
		logObject.info('Successfully ran: %s' % ' '.join(makedb_cmd))
	except Exception as e:
		logObject.error('Had an issue running DIAMOND makedb: %s' % ' '.join(makedb_cmd))
		sys.stderr.write('Had an issue running DIAMOND makedb: %s\n' % ' '.join(makedb_cmd))
		logObject.error(e)
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

	try:
		subprocess.call(' '.join(search_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		assert (os.path.isfile(og_blast_file))
		logObject.info('Successfully ran: %s' % ' '.join(search_cmd))	
	except Exception as e:
		logObject.error('Had an issue running DIAMOND blastp: %s' % ' '.join(search_cmd))
		sys.stderr.write('Had an issue running DIAMOND blastp: %s\n' % ' '.join(search_cmd))
		logObject.error(e)
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def determineFaiParamRecommendataions(genbanks, ortho_matrix_file, hg_prot_dir, outdir, logObject, threads=1):
	"""
	Description:
	Function to determine parameter recommendations for running fai from known instances based on zol orthology 
	inference.
	********************************************************************************************************************
	Parameters:
	- genbanks: list of gene cluster GenBank files.
	- ortho_matrix_file: zol orthology inference.
	- hg_prot_dir: directory of ortholog group protein sequences, each ortholog group is its own FASTA file. 
	- outdir: workspace directory
	- logObject: A logging object.
	********************************************************************************************************************
	"""
	try:
		lt_to_og = {}
		og_gbks = defaultdict(set)
		og_gbk_lts = defaultdict(lambda: defaultdict(set))
		og_lts = defaultdict(set)
		gbks = []
		with open(ortho_matrix_file) as oomf:
			for i, line in enumerate(oomf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					gbks = ls[1:]
				else:
					og = ls[0]
					for j, lts in enumerate(ls[1:]):
						gbk = gbks[j]
						for lt in lts.split(', '):
							if lt.strip() != '':
								lt = '|'.join(lt.split('|')[1:])
								lt_to_og[lt] = og
								og_lts[og].add(lt)
								og_gbk_lts[og][gbk].add(lt)
								og_gbks[og].add(gbk)
		
		self_blast_dir = outdir + 'Self_Blast_OG_Proteins/' 
		setupReadyDirectory([self_blast_dir])
		diamond_self_blasting_inputs = []
		for f in os.listdir(hg_prot_dir):
			og = f.split('.faa')[0]
			og_prot_file = hg_prot_dir + f 
			og_blast_file = self_blast_dir + og + '.txt'
			og_prot_dmnd_db = hg_prot_dir + f
			diamond_self_blasting_inputs.append([og_prot_file, og_prot_dmnd_db, og_blast_file, logObject])
		
		p = multiprocessing.Pool(threads)
		p.map(diamondBlast, diamond_self_blasting_inputs)
		p.close()

		og_min_eval_params = outdir + 'OG_Information.txt'
		omep_handle = open(og_min_eval_params, 'w')
		maximum_evalues = []
		near_core_ogs_maximum_evalues = []
		omep_handle.write('\t'.join(['og', 'maximum_evalue', 'conservation', 'is_near_core', 'is_single_copy']) + '\n')
		near_core_ogs = set([])
		for f in os.listdir(self_blast_dir):
			self_blast_file = self_blast_dir + f 
			og = f.split('.txt')[0]
			maximum_evalue = -1.0
			with open(self_blast_file) as osbf:
				for line in osbf:
					line = line.strip()
					ls = line.split('\t')
					evalue = float(ls[-2])
					if evalue > maximum_evalue:
						maximum_evalue = evalue
			maximum_evalues.append(maximum_evalue)

			og_gbk_count = 0			
			singlecopy = True
			for gbk in og_gbk_lts[og]:
				og_gbk_count += 1
				lt_count = 0
				for lt in og_gbk_lts[og][gbk]:
					lt_count += 1
				if lt_count > 1:
					singlecopy = False
			
			conservation = og_gbk_count/float(len(genbanks))
			nearcore = False
			if conservation >= 0.8:
				nearcore = True
				near_core_ogs.add(og)
				near_core_ogs_maximum_evalues.append(maximum_evalue)

			omep_handle.write(og + '\t' + str(maximum_evalue) + '\t' + str(conservation) + '\t' + str(nearcore) + '\t' + str(singlecopy) + '\n')
		omep_handle.close()

		og_list = defaultdict(list)
		cds_counts = []
		prop_genes_nc = []
		gbk_og_counts = defaultdict(lambda: defaultdict(int))
		for gbk in genbanks:
			cds_count = 0
			nc_cds_count = 0
			with open(gbk) as ogbk:
				for i, rec in enumerate(SeqIO.parse(ogbk, 'genbank')):
					for feature in rec.features:
						if feature.type != 'CDS': continue
						lt = feature.qualifiers.get('locus_tag')[0]
						loc_str = str(feature.location)
						all_coords, start, end, direction, is_multi_part = parseFeatureCoord(loc_str)
						og = lt_to_og[lt]
						gbk_og_counts[gbk][og] += 1
						nc = False
						if og in near_core_ogs:
							nc = True
							nc_cds_count += 1
						og_list[gbk + '|' + str(i)].append([lt, og, start, nc])
						cds_count += 1
			prop_genes_nc.append(nc_cds_count/float(cds_count))
			cds_counts.append(cds_count)

		cds_between_ncs = []
		for rec in og_list:
			last_nc = None
			for i, lt_info in enumerate(sorted(og_list[rec], key=itemgetter(2))):
				is_nc = lt_info[3]
				if is_nc:
					if last_nc != None:
						cbn = i - last_nc
						cds_between_ncs.append(cbn)
					last_nc = i

		gbk_scores = {}
		for i, gbk1 in enumerate(genbanks):
			go1cs = gbk_og_counts[gbk1]
			g1ogs = set(go1cs.keys())
			sum_jaccard = 0.0
			for j, gbk2 in enumerate(genbanks):
				go2cs = gbk_og_counts[gbk2]
				g2ogs = set(go2cs.keys())
				ogs_intersect = g1ogs.intersection(g2ogs)
				ogs_union = g1ogs.union(g2ogs)
				ogs_jaccard = len(ogs_intersect)/float(len(ogs_union))
				sum_jaccard += ogs_jaccard
			gbk_scores[gbk1] = sum_jaccard
	
		ref_gbk = None
		for i, gbk in enumerate(sorted(gbk_scores.items(), key=itemgetter(1), reverse=True)):
			if i == 0:
				ref_gbk = gbk[0]

		# extract near-core OG sequences from ref GBK
		near_core_prots_faa_file = outdir + 'NearCore_Proteins_from_Representative.faa'
		ncpff_handle = open(near_core_prots_faa_file, 'w')
		ref_cds_nc = 0
		ref_cds = 0
		with open(ref_gbk) as orgf:
			for rec in SeqIO.parse(orgf, 'genbank'): 
				for feature in rec.features:
					if not feature.type == 'CDS': continue
					lt = feature.qualifiers.get('locus_tag')[0]
					seq = feature.qualifiers.get('translation')[0]
					ref_cds += 1
					og = lt_to_og[lt]
					if lt in near_core_ogs:
						ref_cds_nc += 1
						ncpff_handle.write('>' + lt + '\n' + seq + '\n')
		ncpff_handle.close()
 	
		prop_ref_cds_nc = ref_cds_nc/float(ref_cds)
		max_distance_between_ncs = max(cds_between_ncs)
		median_cds_count = statistics.median(cds_counts)
		maximum_of_maximum_evalues = max(maximum_evalues)
		maximum_of_near_core_ogs_maximum_evalues = max(near_core_ogs_maximum_evalues)
		median_prop_cds_nc = statistics.median(prop_genes_nc)

		parameter_recommendations_file = outdir + 'Parameter_Recommendations_for_fai.txt'
		prf_handle = open(parameter_recommendations_file, 'w')

		prf_handle.write('=============================================================\n')
		prf_handle.write('Recommendations for running fai to find additional instances of gene cluster:\n')
		prf_handle.write('-------------------------------------------------------------\n')
		prf_handle.write('Note, this functionality assumes that the known instances of the gene cluster\nare representative of the gene cluster/taxonomic diversity you will be searching.\n')
		prf_handle.write('=============================================================\n')
		prf_handle.write('General statistics:\n')
		prf_handle.write('=============================================================\n')
		prf_handle.write('Maximum of maximum E-values observed for any OG\t%s\n' % str(maximum_of_maximum_evalues))
		prf_handle.write('Maximum of near-core OG E-values observed:\t%s\n' % str(maximum_of_near_core_ogs_maximum_evalues))
		prf_handle.write('Maximum distance between near-core OGs:\t%s\n' % str(max_distance_between_ncs))
		prf_handle.write('Median CDS count:\t%s\n' % str(median_cds_count))
		prf_handle.write('Median proportion of CDS which are near-core (conserved in 80 percent of gene-clusters):\t%s\n' % str(median_prop_cds_nc))
		prf_handle.write('Best representative query gene-cluster instance to use:\t%s\n' % ref_gbk)
		prf_handle.write('=============================================================\n')
		prf_handle.write('Parameter recommendations - threads set to 4 by default\n')
		prf_handle.write('please provide the path to the prepTG database yourself!\n')
		prf_handle.write('=============================================================\n')
		prf_handle.write('Lenient / Sensitive Recommendations for Exploratory Analysis:\n')
		fai_cmd = ['fai', '--threads', '4', '--output_dir', 'fai_Search_Results/', '--draft_mode', 
			       '--evalue_cutoff', str(max(maximum_of_maximum_evalues, 1e-10)), '--min_prop', str(max(prop_ref_cds_nc-0.25, 0.1)), 
				   '--syntenic_correlation_threshold', '0.0', '--max_genes_disconnect', str(max_distance_between_ncs+3)]
		prf_handle.write(' '.join(fai_cmd) + '\n')
		prf_handle.write('-------------------------------------------------------------\n')
		prf_handle.write('Strict / Specific Recommendations:\n')
		fai_cmd = ['fai', '--threads', '4', '--output_dir', 'fai_Search_Results/', '--draft_mode', '--filter_paralogs',
			       '--evalue_cutoff', str(maximum_of_maximum_evalues), '--min_prop', str(max(prop_ref_cds_nc, 0.25)), 
				   '--syntenic_correlation_threshold', '0.0', '--max_genes_disconnect', str(max_distance_between_ncs),
				   'key_protein_queries', near_core_prots_faa_file, '--key_protein_min_prop', '0.5', 
				   '--key_protein_evalue_cutoff', str(maximum_of_near_core_ogs_maximum_evalues)]
		prf_handle.write(' '.join(fai_cmd) + '\n')
		prf_handle.write('-------------------------------------------------------------\n')
		prf_handle.close()

		os.system('cat %s' % parameter_recommendations_file)

	except Exception as e:
		sys.stderr.write('Issue with determining parameter recommendations for running fai based on quick zol analysis of known gene cluster instances.\n')
		sys.stderr.write(traceback.format_exc())
		logObject.error('Issue with determining parameter recommendations for running fai based on quick zol analysis of known gene cluster instances.\n')
		logObject.error(traceback.format_exc())
		sys.exit(1)


def parseFeatureCoord(str_gbk_loc):
	try:
		start = None
		end = None
		direction = None
		all_coords = []
		is_multi_part = False
		if not 'join' in str(str_gbk_loc) and not 'order' in str(str_gbk_loc):
			start = min([int(x.strip('>').strip('<')) for x in
						 str(str_gbk_loc)[1:].split(']')[0].split(':')]) + 1
			end = max([int(x.strip('>').strip('<')) for x in
					   str(str_gbk_loc)[1:].split(']')[0].split(':')])
			direction = str(str_gbk_loc).split('(')[1].split(')')[0]
			all_coords.append([start, end, direction])
		elif 'order' in str(str_gbk_loc):
			is_multi_part = True
			all_starts = []
			all_ends = []
			all_directions = []
			for exon_coord in str(str_gbk_loc)[6:-1].split(', '):
				ec_start = min(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
				ec_end = max(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
				ec_direction = exon_coord.split('(')[1].split(')')[0]
				all_starts.append(ec_start)
				all_ends.append(ec_end)
				all_directions.append(ec_direction)
				all_coords.append([ec_start, ec_end, ec_direction])
			assert (len(set(all_directions)) == 1)
			start = min(all_starts)
			end = max(all_ends)
			direction = all_directions[0]
		else:
			is_multi_part = True
			all_starts = []
			all_ends = []
			all_directions = []
			for exon_coord in str(str_gbk_loc)[5:-1].split(', '):
				ec_start = min(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
				ec_end = max(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
				ec_direction = exon_coord.split('(')[1].split(')')[0]
				all_starts.append(ec_start)
				all_ends.append(ec_end)
				all_directions.append(ec_direction)
				all_coords.append([ec_start, ec_end, ec_direction])
			assert (len(set(all_directions)) == 1)
			start = min(all_starts)
			end = max(all_ends)
			direction = all_directions[0]
		return(all_coords, start, end, direction, is_multi_part)
	except Exception as e:
		raise RuntimeError(traceback.format_exc())

def runPyHmmerForRiboProts(best_tg_gbk_file, tg_query_prots_file, ribo_norm_dir, logObject, threads=1):
	"""
	Description:
	Annotate ribosomal proteins from Hug et al. 2016 (using HMMs as provided in GToTree by Lee 2019) in a reference genome.
	********************************************************************************************************************
	Parameters:
	- best_tg_gbk_file: The full genome GenBank file for the reference/select genome.
	- tg_query_prots_file: The FASTA file to write ribosomal proteins identified to (note will append to file).
	- ribo_norm_dir: Output workspace to write intermediate files to.
	- logObject: A logging object.
	- threads: The number of threads to use [Default is 1].
	********************************************************************************************************************
	"""
	try:
		pfam_db_file = None
		try:
			zol_data_directory = str(os.getenv("ZOL_DATA_PATH")).strip()
			db_locations = None
			if zol_data_directory != 'None':
				try:
					zol_data_directory = os.path.abspath(zol_data_directory) + '/'
					db_locations = zol_data_directory + 'database_location_paths.txt'
				except:
					pass
	
			if db_locations == None or not os.path.isfile(db_locations):
				sys.stderr.write('Warning: databases do not appear to be setup or setup properly - so unable to annotate!\n')
	
			with open(db_locations) as odl:
				for line in odl:
					line = line.strip()
					if len(line.split('\t')) != 4: continue
					name, _, db_file, _ = line.split('\t')
					if name == 'pfam': 
						pfam_db_file = db_file
	
			assert(rp_db_file != None and os.path.isfile(rp_db_file))
		except:
			msg = 'Issues validating that the ribosomal proteins file is available. Downloading from Zenodo to output directory.'
			logObject.warning(msg)
			sys.stderr.write(msg + '\n')
			
			curr_path = os.path.abspath(os.getcwd()) + '/'
			download_path = '/'.join((db_locations.split('/')[:-1])) + '/'
			download_links = ['https://zenodo.org/record/7860735/files/Universal_Hug_et_al.hmm?download=1']

			# Download
			os.chdir(download_path)
			rp_db_file = download_path + 'Universal_Hug_et_al.hmm'
			try:
				for dl in download_links:
					axel_download_dbs_cmd = ['axel', '-a', '-n', str(threads), dl]
					os.system(' '.join(axel_download_dbs_cmd))
					assert(os.path.isfile(rp_db_file))
					os.system(' '.join(['hmmpress', rp_db_file]))
			except Exception as e:
				sys.stderr.write('Error occurred during downloading!\n')
				logObject.error('Error occurred during downloading with axel.\n')
				sys.exit(1)
			os.chdir(curr_path)

		hmm_lengths = {}
		z = 0
		try:
			with pyhmmer.plan7.HMMFile(rp_db_file) as hmm_file:
				for hmm in hmm_file:
					hmm_lengths[hmm.name] = len(hmm.consensus)
					z += 1
		except:
			raise RuntimeError("Problem getting HMM consensus lengths!")

		best_tg_faa_file = ribo_norm_dir + 'Reference_Genome_All_Proteins.faa'
		best_tg_faa_handle = open(best_tg_faa_file, 'w')
		try:
			with open(best_tg_gbk_file) as obtgf:
				for rec in SeqIO.parse(obtgf, 'genbank'):
					for feat in rec.features:
						if feat.type == 'CDS':
							lt = feat.qualifiers.get('locus_tag')[0]
							prot_seq = feat.qualifiers.get('translation')[0]
							best_tg_faa_handle.write('>' + lt + '\n' + prot_seq + '\n')
		except:
			raise RuntimeError("Problem processing full genome GenBank file.")
		best_tg_faa_handle.close()
		
		alphabet = pyhmmer.easel.Alphabet.amino()
		sequences = []
		with pyhmmer.easel.SequenceFile(best_tg_faa_file, digital=True, alphabet=alphabet) as seq_file:
			sequences = list(seq_file)

		reference_ribo_prots = set([])
		with pyhmmer.plan7.HMMFile(rp_db_file) as hmm_file:
			for hits in pyhmmer.hmmsearch(hmm_file, sequences, bit_cutoffs='trusted', Z=int(z), cpus=threads):
				for hit in hits:
					# solution for calcualting coverage taken from pcamargo's answer in a pyhmmer ticket on Github: https://github.com/althonos/pyhmmer/issues/27
					n_aligned_positions = len(hit.best_domain.alignment.hmm_sequence) - hit.best_domain.alignment.hmm_sequence.count(".")
					hmm_coverage = (n_aligned_positions / hmm_lengths[hit.best_domain.alignment.hmm_name])
					if hmm_coverage >= 0.25:
						reference_ribo_prots.add(hit.name.decode())
			
		tg_query_prots_handle = open(tg_query_prots_file, 'a+')
		with open(best_tg_faa_file) as obtff:
			for rec in SeqIO.parse(obtff, 'fasta'):
				if rec.id in reference_ribo_prots:
					tg_query_prots_handle.write('>Ribosomal_Protein|' + rec.id + '\n' + str(rec.seq) + '\n')
		tg_query_prots_handle.close()

	except:
		raise RuntimeError('Problem with running pyhmmer for finding ribosomal proteins in reference genome!')

def runPyHmmerForVOGforSalt(inputs):
	"""
	Description:
	Annotate a single sample's predicted proteome file for VOG HMMs.
	********************************************************************************************************************
	Parameters:
	- inputs: A list of length 5:
		- db_file: HMM database.
		- z: Size of database, used for E-value computation.
		- protein_faa: sample's proteome FASTA file.
		- annotation_result_file: Path to output file where to write annotation information.
		- threads: number of threads to use for search.
	********************************************************************************************************************
	"""
	db_file, z, protein_faa, annotation_result_file, threads = inputs
	try:
		hmm_lengths = {}
		try:
			with pyhmmer.plan7.HMMFile(db_file) as hmm_file:
				for hmm in hmm_file:
					hmm_lengths[hmm.name] = len(hmm.consensus)
		except:
			raise RuntimeError("Problem getting HMM consensus lengths!")

		alphabet = pyhmmer.easel.Alphabet.amino()
		sequences = []
		with pyhmmer.easel.SequenceFile(protein_faa, digital=True, alphabet=alphabet) as seq_file:
			sequences = list(seq_file)

		outf = open(annotation_result_file, 'w')
		with pyhmmer.plan7.HMMFile(db_file) as hmm_file:
			for hits in pyhmmer.hmmsearch(hmm_file, sequences, Z=int(z), cpus=threads):
				for hit in hits:
					# solution for calcualting coverage taken from pcamargo's answer in a pyhmmer ticket on Github: https://github.com/althonos/pyhmmer/issues/27
					n_aligned_positions = len(hit.best_domain.alignment.hmm_sequence) - hit.best_domain.alignment.hmm_sequence.count(".")
					hmm_coverage = (n_aligned_positions / hmm_lengths[hit.best_domain.alignment.hmm_name])
					outf.write('\t'.join([hits.query_name.decode(), 'NA', hit.name.decode(), 'NA', str(hit.evalue), str(hit.score), str(hmm_coverage)]) + '\n')
		outf.close()
	except:
		raise RuntimeError('Problem running pyhmmer for annotating MGEs in target genomes!')

			
def annotateMGEs(inputs):
	"""
	Description:
	Annotate MGEs for a single sample's predicted proteome file.
	********************************************************************************************************************
	Parameters:
	- inputs: A list of length 6:
		- sample: sample/genome name.
		- faa_file: path to proteome for sample/genome.
		- vog_annot_file: path to the pyhmmer results for VOG annotations (will be written to - should not exist).
		- mobsuite_annot_file: path to the DIAMOND blastp results for MOB-suite annotations (will be written to - should not exist).
		- is_annot_file: path to the DIAMOND blastp results for ISfinder annotations (will be written to - should not exist).
		- logObject: a logging object.
	********************************************************************************************************************
	"""
	sample = 'NA'
	try:
		sample, faa_file, vog_annot_file, mobsuite_annot_file, is_annot_file, logObject = inputs

		zol_data_directory = str(os.getenv("ZOL_DATA_PATH")).strip()
		db_locations = None
		conda_setup_success = None
		if zol_data_directory != 'None':
			try:
				zol_data_directory = os.path.abspath(zol_data_directory) + '/'
				db_locations = zol_data_directory + 'database_location_paths.txt'
			except:
				pass	

		if db_locations == None or not os.path.isfile(db_locations):
			sys.stderr.write('Warning: databases do not appear to be setup or setup properly - so unable to annotate!\n')

		with open(db_locations) as odl:
			for line in odl:
				line = line.strip()
				if len(line.split('\t')) != 4: continue
				name, _, db_file, z = line.split('\t')
				if name == 'vog':
					vog_pyhmmer_input = [db_file, z, faa_file, vog_annot_file, 1]
					runPyHmmerForVOGforSalt(vog_pyhmmer_input)
				elif name == 'mobsuite' or name == 'isfinder':
					annot_result_file = is_annot_file
					if name == 'mobsuite':
						annot_result_file = mobsuite_annot_file
					search_cmd = ['diamond', 'blastp', '--ignore-warnings', '-fast', '--outfmt', '6', 'qseqid', 'sseqid',
					   'pident', 'evalue', 'qcovhsp', '-p', '1', '-d', db_file, '-q', faa_file, '-o', annot_result_file]
					runCmdViaSubprocess(search_cmd, logObject, check_files=[annot_result_file])
	except Exception as e:
		sys.stderr.write('Issues with MGE annotation commands for sample %s.\n' % sample)
		logObject.error('Issues with MGE annotation commands for sample %s.' % sample)
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.exit(1)

def processMGEAnnotations(inputs):
	"""
	Description:
	Function to process MGE annotation results from DIAMOND blastp and pyhmmer searching for a single sample as well as 
	corresponding genomic GenBank file to determine summarized information for 
	********************************************************************************************************************
	Parameters:
	- inputs: A list of length 6:
		- sample: sample/genome name.
		- gbk_file: path to GenBank for sample/genome.
		- summary_file: The output file for the sample with information on MGE locations.
		- vog_annot_file: path to the pyhmmer results for VOG annotations (should already exist).
		- mobsuite_annot_file: path to the DIAMOND blastp results for MOB-suite annotations (should already exist).
		- is_annot_file: path to the DIAMOND blastp results for ISfinder annotations (should already exist).
		- logObject: a logging object.
	********************************************************************************************************************
	"""
	try:
		# CURRENTLY HARDCODED! 
		max_dmnd_annotation_evalue = 1e-5
		max_hmm_annotation_evalue = 1e-5
		min_percent_identity = 40.0 # only used for diamond
		min_coverage = 70.0
		sample, gbk_file, summary_file, vog_annot_file, mobsuite_annot_file, is_annot_file, logObject = inputs

		try:
			assert(os.path.isfile(vog_annot_file) and os.path.isfile(mobsuite_annot_file) and os.path.isfile(is_annot_file))
		except:
			sys.stderr.write('Issue validating the existence of one or more of the annotation files for sample: %s\n' % sample)
			logObject.error('Issue validating the existence of one or more of the annotation files for sample: %s' % sample)
			# TODO: Consider making this an end-program error
			return 
		
		vog_hits = set([])
		with open(vog_annot_file) as ovaf:
			for line in ovaf:
				line = line.rstrip('\n')
				if line.startswith('#'): continue
				ls = line.split('\t')
				_, _, query, _, evalue, score, coverage = ls
				evalue = decimal.Decimal(evalue)
				coverage = float(coverage)*100.0
				if evalue <= max_hmm_annotation_evalue and coverage >= min_coverage: 
					vog_hits.add(query)
		
		mobsuite_hits = set([])
		with open(mobsuite_annot_file) as omaf:
			for line in omaf:
				line = line.strip()
				query, _, pident, evalue, qcovhsp = line.split('\t')
				pident = float(pident)
				evalue = decimal.Decimal(evalue)
				qcovhsp = float(pident)
				if qcovhsp >= min_coverage and evalue <= max_dmnd_annotation_evalue and pident >= min_percent_identity:
					mobsuite_hits.add(query) 

		isfinder_hits = set([])
		with open(is_annot_file) as oiaf:
			for line in oiaf:
				line = line.strip()
				query, _, pident, evalue, qcovhsp = line.split('\t')
				pident = float(pident)
				evalue = decimal.Decimal(evalue)
				qcovhsp = float(pident)
				if qcovhsp >= min_coverage and evalue <= max_dmnd_annotation_evalue and pident >= min_percent_identity:
					isfinder_hits.add(query) 

		summary_handle = open(summary_file, 'w')
		summary_handle.write('\t'.join(['sample', 'scaffold', 'total_cds', 'mobsuite_cds_hits', 'mobsuite_cds_prop', 'vog_cds_hits', 'vog_cds_prop', 'is_cds_hits', 'is_cds_prop', 'is_cds_boundaries']) + '\n')
		with open(gbk_file) as ogbf:
			for rec in SeqIO.parse(ogbf, 'genbank'):
				scaffold_id = rec.id
				scaffold_is_boundary_coords = set([])				
				total_cds = 0
				mhits = 0
				vhits = 0
				ihits = 0
				for feat in rec.features:
					if feat.type != 'CDS': continue
					lt = feat.qualifiers.get('locus_tag')[0]
					total_cds += 1
					if lt in isfinder_hits:
						ihits += 1
						start = min([int(x.strip('>').strip('<')) for x in
									 str(feat.location)[1:].split(']')[0].split(':')]) + 1
						end = max([int(x.strip('>').strip('<')) for x in
								   str(feat.location)[1:].split(']')[0].split(':')])
						scaffold_is_boundary_coords.add(start)
						scaffold_is_boundary_coords.add(end)
					if lt in mobsuite_hits: mhits += 1
					if lt in vog_hits: vhits += 1
				if total_cds > 0:
					summary_handle.write('\t'.join( [str(x) for x in [sample, scaffold_id, total_cds, mhits, round(mhits/total_cds, 3), vhits, round(vhits/total_cds, 3), ihits, round(ihits/total_cds, 3)]] + [', '.join([str(y) for y in scaffold_is_boundary_coords]) ]) + '\n')
		summary_handle.close()

	except Exception as e:
		sys.stderr.write('Issues with MGE annotation processing for sample %s.\n' % sample)
		logObject.error('Issues with MGE annotation processing for sample %s.' % sample)
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.exit(1)

def processDiamondForGCtoRiboRatio(diamond_results_file, query_faa, fai_ind_gc_file, gc_to_ribo_aai_stat_file, logObject):
	"""
	Description:
	Function to process results from DIAMOND blastping ribosomal + focal gene-cluster proteins from reference genome 
	to the remainder of target genomes to determine the Beta-RD statistic.
	********************************************************************************************************************
	Parameters:
	- diamond_results_file: DIAMOND blastp results searching for query gene cluster and ribosomal proteins from reference
	                        genome against all target genomes.  
	- query_faa: the query FASTA file containing all ribosomal + gene-cluster protein seqeunces.
	- fai_ind_gc_file: fai individual gene cluster instances tsv.
	- gc_to_ribo_aai_stat_file: Path to output file to write statistics per gene cluster instance.
	- logObject: a logging object.
	********************************************************************************************************************
	"""
	try:
		gc_prots = defaultdict(set)
		genome_gcis = defaultdict(set)
		lt_to_gci = {}
		with open(fai_ind_gc_file) as ofigf:
			for i, line in enumerate(ofigf):
				if i == 0: continue			
				line = line.strip()
				ls = line.split('\t')
				gc_gbk = ls[1]
				gc_gbk_file = gc_gbk.split('/')[-1]
				gci = '.gbk'.join(gc_gbk_file.split('.gbk')[:-1])
				genome = '.gbk'.join(gc_gbk_file.split('.gbk')[:-1])
				if '_fai-gene-cluster' in genome:
					genome = genome.split('_fai-gene-cluster')[0]
				genome_gcis[genome].add(gci)
				with open(gc_gbk) as ogg:
					for rec in SeqIO.parse(ogg, 'genbank'):
						for feat in rec.features:
							if feat.type == 'CDS':
								lt = feat.qualifiers.get('locus_tag')[0]
								gc_prots[genome].add(lt)
								lt_to_gci[genome + '|' + lt] = gci 

		gci_query_top_hits = defaultdict(lambda: defaultdict(lambda: [[], 0.0, []]))
		with open(diamond_results_file) as odrf:
			for line in odrf:
				line = line.strip()
				ls = line.split('\t')
				query = ls[0]
				query_class = query.split('|')[0]
				hit = ls[1]
				genome, lt = hit.split('|')
				if query_class == 'Gene_Cluster_Protein' and not lt in gc_prots[genome]: continue
				identity = float(ls[2])
				bitscore = float(ls[4])
				qcovhsp = float(ls[7])
				if qcovhsp <= 25.0: continue
				if query_class == 'Ribosomal_Protein':
					gcis = genome_gcis[genome]
				else:
					gcis = set([lt_to_gci[hit]])
				for gci in gcis:
					if bitscore > gci_query_top_hits[gci][query][1]:
						gci_query_top_hits[gci][query][0] = [hit]
						gci_query_top_hits[gci][query][1] = bitscore
						gci_query_top_hits[gci][query][2] = [identity]
					elif bitscore == gci_query_top_hits[gci][query][1]:
						gci_query_top_hits[gci][query][0].append(hit)
						gci_query_top_hits[gci][query][2].append(identity)

		all_ribo_prots = set([])
		all_gc_prots = set([])
		with open(query_faa) as oqf:
			for rec in SeqIO.parse(oqf, 'fasta'):
				if rec.id.startswith('Gene_Cluster_Protein|'):
					all_gc_prots.add(rec.id)
				else:
					all_ribo_prots.add(rec.id)

		gci_top_hits = defaultdict(list)
		for gci in gci_query_top_hits:
			for query in gci_query_top_hits[gci]:
				max_identity = 0.0
				max_hit = None
				for i, hit in enumerate(gci_query_top_hits[gci][query][0]):
					if gci_query_top_hits[gci][query][2][i] >= max_identity:
						max_identity = gci_query_top_hits[gci][query][2][i] 
						max_hit = hit
				gci_top_hits[gci].append([query, bitscore, max_hit, max_identity])

		gc_aais = []
		rp_aais = []
		data = []
		for gci in gci_top_hits:
			ribo_identities = []
			gc_identities = []
			accounted_lts = set([])
			hits = gci_top_hits[gci]
			for query_hit_data in sorted(hits, key=itemgetter(1), reverse=True):
				query = query_hit_data[0]
				max_ident = query_hit_data[3]
				if not query_hit_data[2] in accounted_lts:
					if query.split('|')[0] == 'Ribosomal_Protein':
						ribo_identities.append(max_ident)
					else:
						gc_identities.append(max_ident)
				accounted_lts.add(query_hit_data[2])
			if len(ribo_identities) > 0 and len(gc_identities) > 0:
				ribo_aai = round(statistics.mean(ribo_identities),3)
				gc_aai = round(statistics.mean(gc_identities),3)
				ribo_prot_prop = round(len(ribo_identities)/len(all_ribo_prots),3)
				gc_prot_prop = round(len(gc_identities)/len(all_gc_prots),3)
				gc_aais.append(gc_aai)
				rp_aais.append(ribo_aai)
				data.append([gci, ribo_aai, gc_aai, ribo_prot_prop, gc_prot_prop])
		
		slope, intercept, r, p, se = stats.linregress(rp_aais, gc_aais)

		gc_to_ribo_aai_stat_handle = open(gc_to_ribo_aai_stat_file, 'w')
		gc_to_ribo_aai_stat_handle.write('\t'.join(['gc_instance', 'distance_to_linear_regression', 'ribo_aai', 'gc_aai', 'prop_ribo_prots', 'prop_gc_prots']) + '\n')
		for rd in data:
			gci, ribo_aai, gc_aai, ribo_prot_prop, gc_prot_prop = rd
			expected_gc_aai = ribo_aai*slope + intercept
			diff_obs_to_exp = gc_aai - expected_gc_aai
			row_data = [gci, diff_obs_to_exp, ribo_aai, gc_aai, ribo_prot_prop, gc_prot_prop]
			gc_to_ribo_aai_stat_handle.write('\t'.join([str(x) for x in row_data]) + '\n')
		gc_to_ribo_aai_stat_handle.close()
	except Exception as e:
		sys.stderr.write('Issues with processing DIAMOND results for aligning ribosomal and gene cluster proteins from the reference to the remainder of the target genomes.')
		logObject.error('Issues with processing DIAMOND results for aligning ribosomal and gene cluster proteins from the reference to the remainder of the target genomes.')
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.exit(1)


def makeGCvsRiboProtAAIScatterplot(rscript_file, gc_to_ribo_aai_stat_file, gc_ribo_aai_plot_pdf_file, logObject):
	"""
	Description:
	Function to plot the gene cluster to ribosomal protein AAI relationship as a scatterplot.
	********************************************************************************************************************
	Parameters:
	- gc_to_ribo_aai_stat_file: The gene cluster and ribosomal proteins AAI relationship file.
	- gc_ribo_aai_plot_pdf_file: The path to the resulting PDF with the scatterplot showing the GC to ribo proteins AAI
								 relationship.
	- logObject: a logging object.
	********************************************************************************************************************
	"""
	try:
		generateSaltGCvsRiboAAIPlot(rscript_file, gc_to_ribo_aai_stat_file, gc_ribo_aai_plot_pdf_file, logObject)
		plot_cmd = ['Rscript', rscript_path]
		try:
			subprocess.call(' '.join(plot_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(gc_ribo_aai_plot_pdf_file))
			logObject.info('Successfully ran: %s' % ' '.join(plot_cmd))
		except Exception as e:
			logObject.error('Had an issue running R based plotting - potentially because of R setup issues in conda: %s' % ' '.join(plot_cmd))
			sys.stderr.write('Had an issue running R based plotting - potentially because of R setup issues in conda: %s\n' % ' '.join(plot_cmd))
			logObject.error(e)
			sys.stderr.write(traceback.format_exc())
			sys.exit(1)
	except Exception as e:
		sys.stderr.write('Issues with creating a scatterplot of gene cluster vs. ribosomal protein AAIs.\n')
		logObject.error('Issues with creating a scatterplot of gene cluster vs. ribosomal protein AAIs.')
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.exit(1)

def consolidateSaltySpreadsheet(fai_ind_gc_file, genome_gbks, codoff_result_dir, mge_annot_result_file, gc_to_ribo_aai_stat_file, result_tsv_file, result_xlsx_file, logObject):
	"""
	Description:
	Consolidate results from all three analyses to inform on lateral transfer into a table and create an auto-colored
	spreadsheet XLSX file.
	********************************************************************************************************************
	Parameters:
	- fai_ind_gc_file: fai individual gene cluster instances tsv.
	- genome_gbks: A dictionary mapping a genome id/name to a full-genome GenBank file.
	- codoff_result_dir: results directory from running codoff for each gene cluster against their respective background 
				         genomes.
	- mge_annot_result_file: overview file with MGE annotations for scaffolds from target genomes.
	- gc_to_ribo_aai_stat_file: Path to output file to write statistics per gene cluster instance.
	- result_tsv_file: The path to the output file to write the consolidated table to in TSV format.
	- result_xlsx_file: The path to the output file to write the consolidated table to in XLSX format.
	- logObject: a logging object.
	********************************************************************************************************************
	"""
	try:

		header = ['gene cluster (GC) instance', 'GC gbk path', 'genome', 'scaffold', 'scaffold length (bp)', 'scaffold CDS count', 'GC CDS count', 
			      'codoff empirical p-value', 'GC AAI observed - expectation',  'GC AAI between genome and reference genome', 
				  'ribosomal protein AAI between genome and reference genome', 'distance to IS element', 'scaffold CDS proportion IS elements', 
				  'scaffold CDS proportion VOGs', 'scaffold CDS proportion plasmid-associated']

		gc_codoff_pvals = defaultdict(lambda: 'NA')
		for f in os.listdir(codoff_result_dir):
			gc = '.txt'.join(f.split('.txt')[:-1])
			gc_codoff_file = codoff_result_dir + f
			with open(gc_codoff_file) as ogcf:
				for line in ogcf:
					line = line.strip()
					ls = line.split('\t')
					if ls[0] == 'Empirical P-value': 
						gc_codoff_pvals[gc] = ls[1]

		genome_scaffold_annot_info = defaultdict(lambda: defaultdict(lambda: ['NA']*8))
		if mge_annot_result_file != None and os.path.isfile(mge_annot_result_file):
			with open(mge_annot_result_file) as omarf:
				for i, line in enumerate(omarf):
					if i == 0: continue 
					line = line.strip('\n')
					sample, scaffold, total_cds, mobsuite_cds_hits, mobsuite_cds_prop, vog_cds_hits, vog_cds_prop, is_cds_hits, is_cds_prop, is_cds_boundaries = line.split('\t')
					genome_scaffold_annot_info[sample][scaffold] = [total_cds, mobsuite_cds_hits, mobsuite_cds_prop, vog_cds_hits, vog_cds_prop, is_cds_hits, is_cds_prop, is_cds_boundaries] 

		gc_aai_stats = defaultdict(lambda: ['NA']*5)
		with open(gc_to_ribo_aai_stat_file) as ogtrasf:
			for i, line in enumerate(ogtrasf):
				if i == 0: continue
				line = line.strip('\n')
				gc, ratio_statistic, ribo_aai, gc_aai, prop_ribo_prots, prop_gc_prots = line.split('\t')
				gc_aai_stats[gc] = [ratio_statistic, ribo_aai, gc_aai, prop_ribo_prots, prop_gc_prots]

		num_rows = 0
		tsv_outf = open(result_tsv_file, 'w')
		tsv_outf.write('\t'.join(header) + '\n')
		with open(fai_ind_gc_file) as ofigf:
			for i, line in enumerate(ofigf):
				if i == 0: continue		
				line = line.strip()
				ls = line.split('\t')	
				gc_gbk = ls[1]
				gc_gbk_file = gc_gbk.split('/')[-1]
				gc = '.gbk'.join(gc_gbk_file.split('.gbk')[:-1])
				genome = '.gbk'.join(gc_gbk_file.split('.gbk')[:-1])
				if '_fai-gene-cluster' in genome:
					genome = genome.split('_fai-gene-cluster')[0]
	
				scaffold = 'NA'
				gc_cds_count = 0

				cds_coords = {}
				scaffold_lens = {}
				with open(genome_gbks[genome]) as ogg:
					for rec in SeqIO.parse(ogg, 'genbank'):
						scaffold_lens[rec.id] = len(str(rec.seq))
						for feat in rec.features:
							if feat.type == 'CDS':
								lt = feat.qualifiers.get('locus_tag')[0]
								loc_str = str(feat.location)
								all_coords, start, end, direction, is_multi_part = parseFeatureCoord(loc_str)
								cds_coords[lt] = [start, end]
				
				min_gc_coord = 1e100
				max_gc_coord = 0
				with open(gc_gbk) as ogg:
					for rec in SeqIO.parse(ogg, 'genbank'):
						scaffold = rec.id
						for feat in rec.features:
							if feat.type == "CDS":
								gc_cds_count += 1
								lt = feat.qualifiers.get('locus_tag')[0]
								cds_start = cds_coords[lt][0]
								cds_end = cds_coords[lt][1]
								if cds_start < min_gc_coord:
									min_gc_coord = cds_start
								if cds_end > max_gc_coord:
									max_gc_coord = cds_end

				scaffold_cds_count, mobsuite_cds_hits, mobsuite_cds_prop, vog_cds_hits, vog_cds_prop, is_cds_hits, is_cds_prop, is_cds_boundaries = genome_scaffold_annot_info[genome][scaffold]
				ratio_statistic, ribo_aai, gc_aai, prop_ribo_prots, prop_gc_prots = gc_aai_stats[gc]
				gc_codoff_pval = gc_codoff_pvals[gc]
			
				min_is_distance = 1e100
				for is_coord in is_cds_boundaries.split(', '):
					if is_numeric(is_coord):
						is_coord = int(is_coord)
						dist = None 
						if is_coord >= min_gc_coord and is_coord <= max_gc_coord:
							dist = 0
						else:
							dist = min([abs(is_coord-min_gc_coord), abs(is_coord-max_gc_coord)])
						if dist < min_is_distance:
							min_is_distance = dist
				if min_is_distance == 1e100:
					min_is_distance = 'NA'

				scaffold_length = scaffold_lens[scaffold]

				row = [gc, gc_gbk, genome, scaffold, scaffold_length, scaffold_cds_count, gc_cds_count, gc_codoff_pval, ratio_statistic, gc_aai, ribo_aai, min_is_distance, 
		               is_cds_prop, vog_cds_prop, mobsuite_cds_prop]
				tsv_outf.write('\t'.join([str(x) for x in row]) + '\n')
				num_rows += 1
		tsv_outf.close()
			
		# Generate Excel spreadsheet
		writer = pd.ExcelWriter(result_xlsx_file, engine='xlsxwriter')
		workbook = writer.book
		dd_sheet = workbook.add_worksheet('Data Dictionary')
		dd_sheet.write(0, 0, 'Data Dictionary describing columns of "SALT Results" spreadsheet can be found on zol\'s Wiki page at:')
		dd_sheet.write(1, 0, 'https://github.com/Kalan-Lab/zol/wiki/5.4-horizontal-or-lateral-transfer-assessment-of-gene-clusters-using-salt')

		na_format = workbook.add_format({'font_color': '#a6a6a6', 'bg_color': '#FFFFFF', 'italic': True})
		header_format = workbook.add_format({'bold': True, 'text_wrap': True, 'valign': 'top', 'fg_color': '#D7E4BC', 'border': 1})

		numeric_columns = {'total scaffold CDS count', 'GC CDS count', 'codoff empirical p-value', 'GC AAI observed - expectation', 'distance to IS element', 
         			       'scaffold CDS proportion IS elements', 'scaffold CDS proportion VOGs', 'scaffold CDS proportion plasmid-associated', 
						   'GC AAI between genome and reference genome', 'ribosomal protein AAI between genome and reference genome'}

		results_df = loadTableInPandaDataFrame(result_tsv_file, numeric_columns)
		results_df.to_excel(writer, sheet_name='SALT Results', index=False, na_rep="NA")
		worksheet =  writer.sheets['SALT Results']
		worksheet.conditional_format('A2:BA' + str(num_rows+1), {'type': 'cell', 'criteria': '==', 'value': '"NA"', 'format': na_format})
		worksheet.conditional_format('A1:BA1', {'type': 'cell', 'criteria': '!=', 'value': 'NA', 'format': header_format})

		# codoff p-value
		worksheet.conditional_format('H2:H' + str(num_rows+1), {'type': '3_color_scale', 'min_color': "#f07878", 'mid_color': '#f7ca81', 'max_color': "#f7de99", "min_value": 0.0, "mid_value": 0.05, "max_value": 1.0, 'min_type': 'num', 'mid_type': 'num', 'max_type': 'num'})

		# diff
		worksheet.conditional_format('I2:I' + str(num_rows+1), {'type': '3_color_scale', 'min_color': "#917967", 'mid_color': '#FFFFFF', 'max_color': "#ba8dc9", "min_value": -100.0, "mid_value": 0.0, "max_value": 100.0, 'min_type': 'num', 'mid_type': 'num', 'max_type': 'num'})

		# AAI columns 
		worksheet.conditional_format('J2:J' + str(num_rows+1), {'type': '2_color_scale', 'min_color': "#dfedf5", 'max_color': "#8fa9b8", "min_value": 0.0, "max_value": 100.0, 'min_type': 'num', 'max_type': 'num'})
		worksheet.conditional_format('K2:K' + str(num_rows+1), {'type': '2_color_scale', 'min_color': "#dfedf5", 'max_color': "#8fa9b8", "min_value": 0.0, "max_value": 100.0, 'min_type': 'num', 'max_type': 'num'})

		# dist to IS/transposon
		worksheet.conditional_format('L2:L' + str(num_rows+1), {'type': '2_color_scale', 'min_color': "#ed87ad", 'max_color': "#f5c6d8", "min_value": 0.0, "max_value": 10000.0, 'min_type': 'num', 'max_type': 'num'})

		# IS/VOG/plasmid
		worksheet.conditional_format('M2:M' + str(num_rows+1), {'type': '2_color_scale', 'min_color': "#98ad93", 'max_color': "#78b56b", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})
		worksheet.conditional_format('N2:N' + str(num_rows+1), {'type': '2_color_scale', 'min_color': "#98ad93", 'max_color': "#78b56b", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})
		worksheet.conditional_format('O2:O' + str(num_rows+1), {'type': '2_color_scale', 'min_color': "#98ad93", 'max_color': "#78b56b", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})

		workbook.close()

	except Exception as e:
		sys.stderr.write('Issues with creating the final salt XLSX spreadsheet.\n')
		logObject.error('Issues with creating the final salt XLSX spreadsheet.')
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.exit(1)

# R functions
def clusterHeatmapR(medlen_data_file, heatmap_data_file, pdf_file, height, width, rscript_path, logObject):
	try:

		rph = open(rscript_path, 'w')
		rph.write('library(ggplot2)\n')
		rph.write('library(cowplot)\n\n')
	
		rph.write('medlen.data_file <- "' + medlen_data_file + '"\n')
		rph.write('heatmap.data_file <- "' + heatmap_data_file + '"\n')
		rph.write('pdf_file <- "' + pdf_file + '"\n')
		rph.write('height <- as.numeric(' + str(height) + ')\n')
		rph.write('width <- as.numeric(' + str(width) + ')\n\n')

		rph.write('medlen.data <- read.table(medlen.data_file, header=T, sep="\\t")\n')
		rph.write('heatmap.data <- read.table(heatmap.data_file, header=T, sep="\\t")\n\n')

		rph.write('gg_ml <- ggplot(medlen.data, aes(x = reorder(og, og_order), y = med_length)) + theme_classic() + xlab("") + \n')
		rph.write('ylab("Median Length\n(kbp)") + geom_bar(stat="identity", fill="black") + \n')
		rph.write('theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())\n\n')

		rph.write('gg_hm <- ggplot(heatmap.data, aes(x = reorder(og, og_order), y = genbank, fill=as.factor(og_presence), label=copy_count)) + \n')
		rph.write('theme_classic() + xlab("ortholog group IDs in Consensus Order") + ylab("") + geom_tile(color="white", show.legend=F) + geom_text() + \n')
		rph.write('theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + \n')
		rph.write('scale_fill_manual(values=c("#FFFFFF", "#889cbd"))\n\n')

		rph.write('pdf(pdf_file, height=height, width=width)\n')
		rph.write('print(plot_grid(gg_ml, gg_hm, ncol=1, axis="l", align="v", rel_heights=c(1, 4)))\n')
		rph.write('dev.off()\n')
		rph.close()
	except Exception as e:
		sys.stderr.write('Issues with creating Rscript clusterHeatmap.R.\n')
		logObject.error('Issues with creating Rscript clusterHeatmap.R.')
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.exit(1)		

def phyloHeatmapR(phylo_tree_file, heatmap_data_file, pdf_file, height, width, rscript_path, logObject):
	try:
		rph = open(rscript_path, 'w')
		rph.write('library(ggtree)\n')
		rph.write('library(ggplot2)\n')
		rph.write('library(ape)\n')
		rph.write('library(dplyr)\n')
		rph.write('library(aplot)\n\n')

		rph.write('phylo.tree_file <- "' + phylo_tree_file + '"\n')
		rph.write('heatmap.data_file <- "' + heatmap_data_file + '"\n')
		rph.write('pdf_file <- "' + pdf_file + '"\n')
		rph.write('pdf.height <- as.numeric(' + str(height) +')\n')
		rph.write('pdf.width <- as.numeric(' + str(width) + ')\n\n')

		rph.write('phylo.tree <- read.tree(phylo.tree_file)\n')
		rph.write('heatmap.data <- read.table(heatmap.data_file, header=T, sep="\\t")\n\n')

		rph.write('pdf(pdf_file, height=pdf.height, width=pdf.width)\n')
		rph.write('gg_tr <- ggtree(phylo.tree)\n')
		rph.write('gg_hm <- ggplot(heatmap.data, aes(x=query_prot_id, y=label, fill=bitscore)) + \n')
		rph.write('theme_classic() + scale_fill_gradient(low="grey", high="black", na.value="white") + \n')
		rph.write('xlab("Query Proteins/Homolog-groups") + ylab("") + geom_tile(color="white", show.legend=F) + \n')
		rph.write('theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))\n')
		rph.write('gg_hm %>% insert_left(gg_tr, width=0.4)\n')
		rph.write('dev.off()\n')
		rph.close()
	except Exception as e:
		sys.stderr.write('Issues with creating Rscript phyloHeatmap.R.\n')
		logObject.error('Issues with creating Rscript phyloHeatmap.R.')
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.exit(1)		

def plotSegmentsR(que_info_file, seq_data_file, pdf_file, pdf_file2, height, width, rscript_path, logObject):
	try:
		rph = open(rscript_path, 'w')

		rph.write('library(ggplot2)\n')
		rph.write('library(gridExtra)\n\n')

		rph.write('que.info.file <- "' + que_info_file + '"\n')
		rph.write('seg.data.file <- "' + seq_data_file + '"\n')
		rph.write('pdf_file <- "' + pdf_file + '"\n')
		rph.write('pdf_file2 <- "' + pdf_file2 + '"\n')
		rph.write('height <- as.numeric(' + str(height) + ')\n')
		rph.write('width <- as.numeric(' + str(width) + ')\n\n')

		rph.write('naming.data <- read.table(que.info.file, header=T, sep="\\t")\n')
		rph.write('pdf(pdf_file2, height=30, width=10)\n')
		rph.write('grid.table(naming.data)\n')
		rph.write('dev.off()\n\n')

		rph.write('segment.data <- read.table(seg.data.file, header=T, sep="\\t")\n')
		rph.write('colors <- c("#000000", "#FFFFFF")\n')
		rph.write('names(colors) <- c("True", "False")\n\n')

		rph.write('samples <- unique(segment.data$sample)\n')
		rph.write('pdf(pdf_file, height=height, width=width)\n')
		rph.write('for (s in samples) {\n')
		rph.write('sample.segment.data <- segment.data[segment.data$sample==s,]\n')
		rph.write('print(unique(sample.segment.data$segment_title))\n')
		rph.write('g<-ggplot(sample.segment.data, aes(x=reorder(gene, gene_order), y=sql_ratio, fill=identity, color=key)) + \n')
		rph.write('geom_bar(stat="identity", position="dodge", size=1.5) + geom_hline(yintercept=1.0, color="blue", linetype=2) + \n')
		rph.write('ylab("CDS to Query Length Ratio") + facet_grid(.~segment_title, space="free", scales="free",  labeller = label_wrap_gen(width = 50, multi_line = TRUE)) + \n')
		rph.write('xlab("CDS Classifications") + theme_bw() + theme(legend.position="bottom", axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) + \n')
		rph.write('scale_color_manual(values=colors) + scale_fill_gradient(limits=c(0.0,100.0), breaks=c(0.0, 50.0, 100.0), low="grey", high="red")\n')
		rph.write('print(g)\n')
		rph.write('}\n')
		rph.write('dev.off()\n\n')

		rph.close()
	except Exception as e:
		sys.stderr.write('Issues with creating Rscript plotSegments.R.\n')
		logObject.error('Issues with creating Rscript plotSegments.R.')
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.exit(1)		

def plotTinyAAIR(info_file, pdf_file, rscript_path, logObject):
	try:
		rph = open(rscript_path, 'w')

		rph.write('library(ggplot2)\n')
		rph.write('library(gridExtra)\n')

		rph.write('info.file <- "' + info_file + '"\n') 
		rph.write('pdf_file <- "' + pdf_file + '"\n\n')
		
		rph.write('info.dat <- read.table(info.file, header=T, sep="\\t")\n\n')

		rph.write('pdf(pdf_file, height=10, width=10)\n')
		rph.write('ggplot(info.dat, aes(x=AAI, y=Prop_Genes_Found, color=Mean_Syntenic_Correlation)) + \n')
		rph.write('geom_point(alpha=0.7) + theme_bw() + scale_color_gradient(low="#e6ffbd", high="#1754b0") + \n')
		rph.write('guides(color=guide_legend("Syntenic\nCorrelation\nto\nQuery")) + \n')
		rph.write('xlab("Average Amino-Acid Identity") + ylab("Proportion of Query Proteins with Match")\n')
		rph.write('dev.off()\n')

		rph.close()
	except:
		sys.stderr.write('Issues with creating Rscript plotTinyAAI.R.\n')
		logObject.error('Issues with creating Rscript plotTinyAAI.R.')
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.exit(1)		

def generateSyntenicVisualR(input_file, pdf_file, height, width, rscript_path):
	try:
		rph = open(rscript_path, 'w')
	
		rph.write('library(ggplot2)\n')
		rph.write('library(gggenes)\n')

		rph.write('input.file <- "' + input_file + '"\n')
		rph.write('pdf.file <- "' + pdf_file + '"\n')
		rph.write('height <- as.numeric(' + str(height) + ')\n')
		rph.write('width <- as.numeric(' + str(width) + ')\n\n')
	
		rph.write('input.data <- read.table(file=input.file, sep="\\t", header=T)\n\n')

		rph.write('pdf(pdf.file, height=height, width=width)\n')
		rph.write('ggplot(input.data, aes(xmin=Start, xmax = End, y = "", forward = Direction, label=SC)) + \n')
		rph.write('geom_gene_arrow(aes(fill=Metric)) + theme_classic() + \n')
		rph.write('scale_fill_gradient2(low="#e05c74", mid="#f2f2f2", high="#2087b3", breaks =c(-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0),\n')
		rph.write('labels=c("", "-2", "", "0", "", "2", ""), limits=c(-3,3), na.value="grey50", guide="colourbar", aesthetics="fill") + \n')
		rph.write('geom_gene_label(align="centre", min.size=5) + theme(legend.position="bottom")\n')
		rph.write('dev.off()\n')

		rph.close()
	except:
		sys.stderr.write('Issues with creating/running Rscript %s.\n' % rscript_path)
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)		

def generateSaltGCvsRiboAAIPlot(rscript_path, input_data_file, pdf_file, logObject):
	try:
		rph = open(rscript_path, 'w')

		rph.write('library(ggplot2)\n\n')

		rph.write('input.data_file <- "' + input_data_file + '"\n')
		rph.write('pdf_file <- "' + pdf_file  + '"\n\n')
			
		rph.write('dat <- read.table(input.data_file, header=T, sep="\\t")\n')
		rph.write('pdf(pdf_file, height=5, width=5)\n')
		rph.write('ggplot(dat, aes(x=ribo_aai, y=gc_aai)) + geom_point(alpha=0.7) +\n')
		rph.write('geom_smooth(method="lm", formula= y~x, linetype=2, color="red") +\n')
		rph.write('geom_abline(slope=1, yintercept=0, linetype=2, color="grey") + theme_bw() +\n')
		rph.write('xlab("Ribosomal Protein AAI") + ylab("Gene Cluster AAI")\n')
		rph.write('dev.off()\n')
	except Exception as e:
		sys.stderr.write('Issues with creating Rscript generateSaltGCVsRiboAAIPlot.R.\n')
		logObject.error('Issues with creating Rscript generateSaltGCVsRiboAAIPlot.R.\n')
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.exit(1)

def generateNjTree(rscript_path, input_dist_file, output_tree_file, logObject):
	try:
		rph.write('library(ape)\n\n')

		rph.write('input.dist.file <- "' + input_dist_file + '"\n')
		rph.write('output.nwk.file <- "' + output_tree_file + '"\n\n')

		rph.write('dat <- read.table(input.dist.file, header=T, sep="\\t", row.names=1)\n')
		rph.write('d <- as.dist(as.matrix(dat))\n')
		rph.write('njt <- nj(d)\n')
		rph.write('write.tree(njt, file=output.nwk.file)\n')

	except Exception as e:
		sys.stderr.write('Issues with creating Rscript generateNjTree.R.\n')
		logObject.error('Issues with creating Rscript generateNjTree.R.')
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.exit(1)
