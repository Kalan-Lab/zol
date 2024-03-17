import os
import sys
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
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
from zol import fai
import statistics

version = pkg_resources.require("zol")[0].version

valid_alleles = set(['A', 'C', 'G', 'T'])

# code for setup and finding location of programs based on conda vs. bioconda installation
zol_exec_directory = str(os.getenv("ZOL_EXEC_PATH")).strip()
conda_setup_success = None
nj_tree_prog = None
if zol_exec_directory != 'None':
	try:
		zol_exec_directory = os.path.abspath(zol_exec_directory) + '/'
		nj_tree_prog = zol_exec_directory + 'njTree.R'
		conda_setup_success = True
	except:
		conda_setup_success = False
if zol_exec_directory == 'None' or conda_setup_success == False:
	nj_tree_prog = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/' + 'zol/njTree.R'

if nj_tree_prog == None or not os.path.isfile(nj_tree_prog):
	sys.stderr.write('Issues in setup of the zol-suite (in util.py) - please describe your installation process and post an issue on GitHub!\n')
	sys.exit(1)

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
		with open(full_genbank_file) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				if not rec.id == scaffold: continue
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
						start = min([int(x.strip('>').strip('<')) for x in
									 str(feature.location)[1:].split(']')[0].split(':')]) + 1
						end = max(
							[int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
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

def setupReadyDirectory(directories):
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
				os.system('rm -rf %s' % d)
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
def parseGenbankForCDSProteinsAndDNA(gbk, logObject, allow_edge_cds=True):
	"""
	Description:
	This function parses GenBank for CDS protein and nucleotide sequences.
	********************************************************************************************************************
	Parameters:
	- gbk: Path to the GenBank file.
	- logObject: A logging object.
	- allow_edge_cds: Whether to regard CDS features near scaffold edges.
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
					if feature.type != 'CDS': continue
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
								additional_proteomes_directory, additional_genbanks_directory, logObject, cpus=1,
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
	- cpus: The number of CPUs to use.
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

		p = multiprocessing.Pool(cpus)
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
								cpus=1, locus_tag_length=3, gene_calling_method="pyrodigal", meta_mode=False,
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
	- cpus: The number of CPUs to use.
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

		p = multiprocessing.Pool(cpus)
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
							 logObject, cpus=1, locus_tag_length=3, avoid_locus_tags=set([]),
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
	- cpus: The number of CPUs to use.
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

		p = multiprocessing.Pool(cpus)
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
def parseSampleGenomes(genome_listing_file, format_assess_dir, format_predictions_file, logObject, cpus=1):
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
	- cpus: The number of CPUs to use.
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

		p = multiprocessing.Pool(cpus)
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

def parseGbk(gbk, prefix, logObject, use_either_lt_or_pi=False):
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

def createNJTree(additional_genbanks_directory, species_tree, workspace_dir, logObject, cpus=1):
	"""
	Description:
	Function to create a species tree using ANI estimates + a neighbor joining approach.
	********************************************************************************************************************
	Parameters:
	- additional_genbanks_directory: Directory with full genomes in GenBank format.
	- workspace_dir: Workspace directory.
	- logObject: A logging object.
	- cpus: The number of CPUs to use.
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
		p = multiprocessing.Pool(cpus)
		p.map(convertGenomeGenBankToFasta, conversion_inputs)
		p.close()

		# run skani triangle
		skani_result_file = workspace_dir + 'Skani_Triangle_Edge_Output.txt'
		skani_triangle_cmd = ['skani', 'triangle', '-E', '-l', all_genomes_listing_file,'-t', str(cpus), '-o', skani_result_file]
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

		unrooted_tree_file = workspace_dir + 'Unrooted_Species_Tree.nwk'
		nj_tree_cmd = ['Rscript', nj_tree_prog, dist_matrix_file, unrooted_tree_file]
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
				os.system('rm -f %s' % df)
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
							   prop_key_prots_needed=0.0, cpus=1):
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
	- cpus: The number of CPUs to use.
	********************************************************************************************************************
	"""
	try:
		genome_wide_tsv_result_file = workspace_dir + 'total_gcs.tsv'
		target_genome_dmnd_db = target_genomes_db + 'Target_Genomes_DB.dmnd'

		# run diamond (use's fai function: runDiamondBlastp)
		diamond_results_file = workspace_dir + 'DIAMOND_Results.txt'
		fai.runDiamondBlastp(target_genome_dmnd_db, query_protein_fasta, workspace_dir, logObject,
					 diamond_sensitivity=blastp_mode, evalue_cutoff=evalue_cutoff, compute_query_coverage=True, 
					 cpus=cpus)

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

def determineFaiParamRecommendataions(genbanks, ortho_matrix_file, hg_prot_dir, outdir, logObject, cpus=1):
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
		
		p = multiprocessing.Pool(cpus)
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
						loc_str = feature.location
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
		prf_handle.write('Parameter recommendations - CPUs set to 4 by default\n')
		prf_handle.write('please provide the path to the prepTG database yourself!\n')
		prf_handle.write('=============================================================\n')
		prf_handle.write('Lenient / Sensitive Recommendations for Exploratory Analysis:\n')
		fai_cmd = ['fai', '--cpus', '4', '--output_dir', 'fai_Search_Results/', '--draft_mode', 
			       '--evalue_cutoff', str(max(maximum_of_maximum_evalues, 1e-10)), '--min_prop', str(max(prop_ref_cds_nc-0.25, 0.1)), 
				   '--syntenic_correlation_threshold', '0.0', '--max_genes_disconnect', str(max_distance_between_ncs+3)]
		prf_handle.write(' '.join(fai_cmd) + '\n')
		prf_handle.write('-------------------------------------------------------------\n')
		prf_handle.write('Strict / Specific Recommendations:\n')
		fai_cmd = ['fai', '--cpus', '4', '--output_dir', 'fai_Search_Results/', '--draft_mode', '--filter_paralogs',
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