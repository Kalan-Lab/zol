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
from scipy import stats
import numpy as np
import gzip
import pathlib
import copy
import itertools
import multiprocessing
import pickle
import resource

valid_alleles = set(['A', 'C', 'G', 'T'])
curr_dir = os.path.abspath(pathlib.Path(__file__).parent.resolve()) + '/'
main_dir = '/'.join(curr_dir.split('/')[:-2]) + '/'

def memory_limit(mem):
	max_virtual_memory = mem*1000000000
	soft, hard = resource.getrlimit(resource.RLIMIT_AS)
	resource.setrlimit(resource.RLIMIT_AS, (max_virtual_memory, hard))
	print(resource.getrlimit(resource.RLIMIT_AS))
def cleanUpSampleName(original_name):
	return original_name.replace('#', '').replace('*', '_').replace(':', '_').replace(';', '_').replace(' ',
																										'_').replace(
		':', '_').replace('|', '_').replace('"', '_').replace("'", '_').replace("=", "_").replace('-', '_').replace('(',
																													'').replace(
		')', '').replace('/', '').replace('\\', '').replace('[', '').replace(']', '').replace(',', '')

def readInAnnotationFilesForExpandedSampleSet(expansion_listing_file, full_dir, logObject=None):
	"""
	Function to read in GenBank paths from expansion listing file and load into dictionary with keys corresponding to sample IDs.
	:param expansion_listing_file: tab-delimited file with three columns: (1) sample ID (2) Genbank path (3) predicted proteome path.
	:param logObject: python logging object handler.
	:return sample_annotation_data: dictionary of dictionaries with primary keys as sample names and secondary keys as either "genbank" or "predicted_proteome", with final values being paths to corresponding files.
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
		raise RuntimeError(traceback.format_exc())

def createGenbank(full_genbank_file, new_genbank_file, scaffold, start_coord, end_coord):
	"""
	Function to prune full genome-sized GenBank for only features in BGC of interest.
	:param full_genbank_file: Prokka generated GenBank file for full genome.
	:param new_genbank_file: Path to BGC specific Genbank to be created
	:param scaffold: Scaffold identifier.
	:param start_coord: Start coordinate.
	:param end_coord: End coordinate.
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
		raise RuntimeError(traceback.format_exc())

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

def is_integer(x):
	try:
		x = int(x)
		return True
	except:
		return False

def is_genbank(gbk, check_for_cds=False):
	"""
	Function to check in Genbank file is correctly formatted.
	"""
	try:
		recs = 0
		cds_flag = False
		assert (gbk.endswith('.gbk') or gbk.endswith('.gbff') or gbk.endswith('.gbk.gz') or gbk.endswith('.gbff.gz'))
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

def checkValidGenBank(genbank_file, quality_assessment=False, draft_assessment=False, use_either_lt_or_pi=False):
	try:
		number_of_cds = 0
		lt_has_comma = False
		lt_is_none = False
		seqs = ''
		recs = 0
		edgy_cds = False
		with open(genbank_file) as ogbk:
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

def parseGenbankForCDSProteinsAndDNA(gbk_path, logObject, allow_edge_cds=True):
	try:
		proteins = {}
		nucleotides = {}
		upstream_regions = {}
		with open(gbk_path) as ogbk:
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

					#try:
					#	final_upstream_region = feature.qualifiers.get('orf_upstream')[0]
					#except:
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

def default_to_regular(d):
	"""
	Function to convert defaultdict of defaultdict to dict of dict
	Taken from Martijn Pieters response in StackOverflow:
	https://stackoverflow.com/questions/26496831/how-to-convert-defaultdict-of-defaultdicts-of-defaultdicts-to-dict-of-dicts-o
	"""
	if isinstance(d, defaultdict):
		d = {k: default_to_regular(v) for k, v in d.items()}
	return d

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

def convertGenbankToCDSProtsFasta(genbank, protein, logObject, use_either_lt_or_pi=False):
	try:
		print(genbank)
		prot_handle = open(protein, 'w')
		with open(genbank) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type != 'CDS': continue
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
		sys.exit(1)
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
		if x == '< 3 segregating sites!':
			return(x)
		else:
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

def checkCoreHomologGroupsExist(ortho_matrix_file):
	"""
	Function to check that a core ortholog/homolog group exists across homologous gene-clusters.
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
								additional_proteomes_directory, addtitional_genbanks_directory, logObject, cpus=1,
								locus_tag_length=3):
	try:
		alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		possible_locustags = sorted(list([''.join(list(x)) for x in list(itertools.product(alphabet, repeat=locus_tag_length))]))

		miniprot_cmds = []
		for i, sample in enumerate(sample_genomes):
			sample_assembly = sample_genomes[sample]
			sample_locus_tag = ''.join(list(possible_locustags[i]))

			sample_mp_db = additional_miniprot_outdir + sample + '.mpi'
			sample_mp_gff = additional_miniprot_outdir + sample + '.gff'
			sample_mp_gbk = addtitional_genbanks_directory + sample + '.gbk'
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
				sample_mp_gbk = addtitional_genbanks_directory + sample + '.gbk'
				sample_mp_faa = additional_proteomes_directory + sample + '.faa'
				assert (os.path.isfile(sample_mp_gbk) and os.path.isfile(sample_mp_faa))
			except:
				raise RuntimeError("Unable to validate successful genbank/predicted-proteome creation for sample %s" % sample)
	except Exception as e:
		logObject.error("Problem with creating commands for running miniprot or convertMiniprotGffToGbkAndProt.py. Exiting now ...")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

def processGenomesUsingProdigal(sample_genomes, prodigal_outdir, prodigal_proteomes, prodigal_genbanks, logObject, cpus=1,
				   locus_tag_length=3, use_prodigal=False, meta_mode=False, avoid_locus_tags=set([])):
	"""
	Void function to run Prodigal based gene-calling and annotations.
	:param sample_genomes: dictionary with keys as sample names and values as genomic assembly paths.
	:param prodigal_outdir: full path to directory where Prokka results will be written.
	:param prodigal_proteomes: full path to directory where Prokka generated predicted-proteome FASTA files will be moved after prodigal has run.
	:param prodigal_genbanks: full path to directory where Prokka generated Genbank (featuring predicted CDS) files will be moved after prodigal has run.
	:param taxa: name of the taxonomic clade of interest.
	:param logObject: python logging object handler.
	:param cpus: number of cpus to use in multiprocessing Prokka cmds.
	:param locus_tag_length: length of locus tags to generate using unique character combinations.
	Note length of locus tag must be 3 beause this is substituting for base lsaBGC analysis!!
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

			prodigal_cmd = ['runProdigalAndMakeProperGenbank.py', '-i', sample_assembly, '-s', sample,
							'-l', sample_locus_tag, '-o', prodigal_outdir]
			if use_prodigal:
				prodigal_cmd += ['-p']
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
				raise RuntimeError(
					"Unable to validate successful genbank/predicted-proteome creation for sample %s" % sample)
	except Exception as e:
		logObject.error(
			"Problem with creating commands for running prodigal via script runProdigalAndMakeProperGenbank.py. Exiting now ...")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())

def processGenomesAsGenbanks(sample_genomes, proteomes_directory, genbanks_directory, gene_name_mapping_outdir,
							 logObject, cpus=1, locus_tag_length=3, avoid_locus_tags=set([]),
							 rename_locus_tags=False):
	"""
	Extracts CDS/proteins from existing Genbank files and recreates
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
				raise RuntimeError(
					"Unable to validate successful genbank/predicted-proteome creation for sample %s" % sample)
	except Exception as e:
		logObject.error("Problem with processing existing Genbanks to (re)create genbanks/proteomes. Exiting now ...")
		logObject.error(traceback.format_exc())
		raise RuntimeError(traceback.format_exc())
	return sample_genomes

def determineGenomeFormat(inputs):
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
		sys.exit(1)

def parseSampleGenomes(genome_listing_file, format_assess_dir, format_predictions_file, logObject, cpus=1):
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
		raise RuntimeError(traceback.format_exc())

def filterRecordsNearScaffoldEdge(genbank_file, filt_genbank_file, logObject, quality_assessment=False):
	try:
		number_of_cds = 0
		seqs = ""
		recs = 0
		recs_with_edgy_cds = set([])
		with open(genbank_file) as ogbk:
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
			with open(genbank_file) as ogbk:
				for rec_it, rec in enumerate(SeqIO.parse(ogbk, 'genbank')):
					if rec_it in recs_with_edgy_cds: continue
					SeqIO.write(rec, out_handle, 'genbank')
			out_handle.close()

	except Exception as e:
		sys.stderr.write('Issue parsing GenBank %s and CDS locus tag renaming.\n' % genbank_file)
		logObject.error('Issue parsing GenBank %s and CDS locus tag renaming.' % genbank_file)
		sys.stderr.write(str(e) + '\n')
		raise RuntimeError(traceback.format_exc())
		sys.exit(1)

def renameCDSLocusTag(genbank_file, lt, rn_genbank_file, logObject, quality_assessment=False, draft_assessment=False):
	try:
		number_of_cds = 0
		seqs = ""
		recs = 0
		recs_with_edgy_cds = set([])
		with open(genbank_file) as ogbk:
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
			with open(genbank_file) as ogbk:
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
		sys.stderr.write('Issue parsing GenBank %s and CDS locus tag renaming.\n' % genbank_file)
		logObject.error('Issue parsing GenBank %s and CDS locus tag renaming.' % genbank_file)
		sys.stderr.write(str(e) + '\n')
		raise RuntimeError(traceback.format_exc())
		sys.exit(1)


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
					location = {'scaffold': rec.id, 'start': start, 'end': end, 'direction': dir}
					gc_gene_locations[prefix + '|' + lt] = location
		return gc_gene_locations
	except Exception as e:
		sys.stderr.write('Issue parsing GenBank %s\n' % gbk_path)
		logObject.error('Issue parsing GenBank %s' % gbk_path)
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def determinePossibleLTs():
	alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
	possible_locustags = sorted(list([''.join(list(lt)) for lt in itertools.product(alphabet, repeat=4)]))
	return possible_locustags

def gatherAnnotationFromDictForHomoloGroup(hg, db, dict):
	try:
		assert(db in dict)
		annot_set_filt = set([x for x in dict[db][hg][0] if x.strip() != ''])
		assert(len(annot_set_filt) > 0)
		return('; '.join(annot_set_filt) + ' (' + str(max(dict[db][hg][1])) + ')')
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

def chunks(lst, n):
	"""
    Yield successive n-sized chunks from lst.
    Solution taken from: https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    """
	for i in range(0, len(lst), n):
		yield lst[i:i + n]

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


	sample_data = [gene_location, dict(scaffold_genes), boundary_genes, dict(gene_id_to_order), dict(gene_order_to_id)]

	with open(pkl_result_file, 'wb') as pickle_file:
		pickle.dump(sample_data, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
