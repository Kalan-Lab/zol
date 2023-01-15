import os
import sys
from Bio import SeqIO
import subprocess
from collections import defaultdict
import multiprocessing
from zol import util
import itertools
import math
from operator import itemgetter
import statistics
import json
import random
import xlsxwriter
import pandas as pd
import traceback

zol_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
stag_prog = zol_main_directory + 'STAG_1.0.0/stag'
def partitionSequencesByHomologGroups(ortho_matrix_file, prot_dir, nucl_dir, hg_prot_dir, hg_nucl_dir, logObject):
	try:
		g_to_hg = {}
		samples = []
		with open(ortho_matrix_file) as oomf:
			for i, line in enumerate(oomf):
				line = line.strip()
				ls = line.split('\t')
				if i == 0:
					samples = ls[1:]
					continue
				hg = ls[0]
				for j, gs in enumerate(ls[1:]):
					for g in gs.split(', '):
						g = g.strip()
						g_to_hg[g] = hg
		for pf in os.listdir(prot_dir):
			pfile = prot_dir + pf
			with open(pfile) as opf:
				for rec in SeqIO.parse(opf, 'fasta'):
					hg = g_to_hg[rec.id]
					hpf_handle = open(hg_prot_dir + hg + '.faa', 'a+')
					hpf_handle.write('>' + rec.description + '\n' + str(rec.seq) + '\n')
					hpf_handle.close()
		for nf in os.listdir(nucl_dir):
			nfile = nucl_dir + nf
			with open(nfile) as onf:
				for rec in SeqIO.parse(onf, 'fasta'):
					hg = g_to_hg[rec.id]
					hnf_handle = open(hg_nucl_dir + hg + '.fna', 'a+')
					hnf_handle.write('>' + rec.description + '\n' + str(rec.seq) + '\n')
					hnf_handle.close()
	except Exception as e:
		sys.stderr.write('Issues with partitioning sequences to homolog groups.\n')
		logObject.error('Issues with partitioning sequences to homolog groups.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)


def createProteinAlignments(prot_dir, prot_algn_dir, logObject, use_super5=False, cpus=1):
	try:
		for pf in os.listdir(prot_dir):
			prefix = '.faa'.join(pf.split('.faa')[:-1])
			prot_file = prot_dir + pf
			prot_algn_file = prot_algn_dir + prefix + '.msa.faa'
			align_cmd = ['muscle', '-align', prot_file, '-output', prot_algn_file, '-amino', '-threads', str(cpus)]
			if use_super5:
				align_cmd = ['muscle', '-super5', prot_file, '-output', prot_algn_file, '-amino', '-threads', str(cpus)]
			try:
				subprocess.call(' '.join(align_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
								executable='/bin/bash')
				assert(os.path.isfile(prot_algn_file))
				logObject.info('Successfully ran: %s' % ' '.join(align_cmd))
			except Exception as e:
				logObject.error('Had an issue running MUSCLE: %s' % ' '.join(align_cmd))
				sys.stderr.write('Had an issue running MUSCLE: %s\n' % ' '.join(align_cmd))
				logObject.error(e)
				sys.exit(1)
	except Exception as e:
		sys.stderr.write('Issues with creating protein alignments.\n')
		logObject.error('Issues with creating protein alignments.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def createCodonAlignments(prot_algn_dir, nucl_dir, codo_algn_dir, logObject, cpus=1):
	try:
		pal2nal_cmds = []
		for paf in os.listdir(prot_algn_dir):
			prefix = '.msa.faa'.join(paf.split('.msa.faa')[:-1])
			prot_algn_file = prot_algn_dir + paf
			nucl_file = nucl_dir + prefix + '.fna'
			codo_algn_file = codo_algn_dir + prefix + '.msa.fna'
			pal2nal_cmds.append(['pal2nal.pl', prot_algn_file, nucl_file, '-output', 'fasta', '>', codo_algn_file, logObject])
		p = multiprocessing.Pool(cpus)
		p.map(util.multiProcess, pal2nal_cmds)
		p.close()
	except Exception as e:
		sys.stderr.write('Issues with creating codon alignments.\n')
		logObject.error('Issues with creating codon alignments.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def trimAlignments(prot_algn_dir, codo_algn_dir, prot_algn_trim_dir, codo_algn_trim_dir, logObject, cpus=1):
	try:
		trim_cmds = []
		for paf in os.listdir(prot_algn_dir):
			prefix = '.msa.faa'.join(paf.split('.msa.faa')[:-1])
			prot_algn_file = prot_algn_dir + paf
			prot_algn_trim_file = prot_algn_trim_dir + paf
			codo_algn_file = codo_algn_dir + prefix + '.msa.fna'
			codo_algn_trim_file = codo_algn_trim_dir + prefix + '.msa.fna'
			trim_cmds.append(['trimal', '-in', prot_algn_file, '-out', prot_algn_trim_file, '-keepseqs', '-gt', '0.9', logObject])
			trim_cmds.append(['trimal', '-in', codo_algn_file, '-out', codo_algn_trim_file, '-keepseqs', '-gt', '0.9', logObject])
		p = multiprocessing.Pool(cpus)
		p.map(util.multiProcess, trim_cmds)
		p.close()
	except Exception as e:
		sys.stderr.write('Issues with trimming protein/codon alignments.\n')
		logObject.error('Issues with trimming protein/codon alignments.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def createGeneTrees(codo_algn_trim_dir, tree_dir, logObject, cpus=1):
	try:
		fasttree_cmds = []
		for catf in os.listdir(codo_algn_trim_dir):
			if not catf.endswith('.msa.fna'): continue
			prefix = '.msa.faa'.join(catf.split('.msa.fna')[:-1])
			codo_algn_trim_file = codo_algn_trim_dir + catf
			tree_file = tree_dir + prefix + '.tre'
			fasttree_cmds.append(['fasttree', '-nt', codo_algn_trim_file, '>', tree_file, logObject])
		p = multiprocessing.Pool(cpus)
		p.map(util.multiProcess, fasttree_cmds)
		p.close()
	except Exception as e:
		sys.stderr.write('Issues with creating gene-trees.\n')
		logObject.error('Issues with creating gene-trees.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def createProfileHMMsAndConsensusSeqs(prot_algn_dir, phmm_dir, cons_dir, logObject, cpus=1):
	try:
		hmmbuild_cmds = []
		hmmemit_cmds = []
		for paf in os.listdir(prot_algn_dir):
			prefix = '.msa.faa'.join(paf.split('.msa.faa')[:-1])
			prot_algn_file = prot_algn_dir + paf
			prot_hmm_file = phmm_dir + prefix + '.hmm'
			prot_cons_file = cons_dir + prefix + '.cons.faa'
			hmmbuild_cmds.append(['hmmbuild', '--amino', '--cpu', '2', '-n', prefix, prot_hmm_file, prot_algn_file, logObject])
			hmmemit_cmds.append(['hmmemit', '-c', '-o', prot_cons_file, prot_hmm_file, logObject])
		p = multiprocessing.Pool(cpus)
		p.map(util.multiProcess, hmmbuild_cmds)
		p.map(util.multiProcess, hmmemit_cmds)
		p.close()
	except Exception as e:
		sys.stderr.write('Issues with creating profile HMMs and consensus sequences.\n')
		logObject.error('Issues with creating profile HMMs and consensus sequences.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def annotateConsensusSequences(protein_faa, annotation_dir, logObject, cpus=1, min_annotation_evalue=1e-10):
	db_locations = zol_main_directory + 'db/database_location_paths.txt'
	try:
		individual_cpus = 1
		pool_size = cpus
		if cpus > 9:
			individual_cpus = math.floor(cpus/9)
			pool_size = 9
		assert(os.path.isfile(db_locations))
		search_cmds = []
		name_to_info_file = {}
		hmm_based_annotations = set([])
		with open(db_locations) as odls:
			for line in odls:
				line = line.strip()
				name, annot_info_file, db_file = line.split('\t')
				name_to_info_file[name] = annot_info_file
				annotation_result_file = annotation_dir + name + '.txt'
				if db_file.endswith('.hmm'):
					hmm_based_annotations.add(name)
					search_cmd = ['hmmscan', '--cpu', str(individual_cpus), '--tblout', annotation_result_file,
								  db_file, protein_faa, logObject]
				elif db_file.endswith('.dmnd'):
					search_cmd = ['diamond', 'blastp', '-p', str(individual_cpus), '-d', db_file, '-q', protein_faa,
								  '-o', annotation_result_file, logObject]
				search_cmds.append(search_cmd)

		p = multiprocessing.Pool(pool_size)
		p.map(util.multiProcess, search_cmds)
		p.close()

		annotations = defaultdict(lambda: defaultdict(lambda: ['NA', 'NA'])) # db -> query -> [hit descriptions, evalue]
		for rf in os.listdir(annotation_dir):
			db_name = rf.split('.txt')[0]
			annot_info_file = name_to_info_file[db_name]

			id_to_description = {}
			with open(annot_info_file) as oaif:
				for line in oaif:
					line = line.strip()
					ls = line.split('\t')
					id_to_description[ls[0]] = ls[1]

			best_hits_by_evalue = defaultdict(lambda: [[], min_annotation_evalue])
			if db_name in hmm_based_annotations:
				# parse HMM based results from HMMER3
				with open(annotation_dir + rf) as oarf:
					for line in oarf:
						line = line.rstrip('\n')
						if line.startswith('#'): continue
						ls = line.split()
						query = ls[2]
						hit = ls[0]
						evalue = float(ls[4])
						if db_name != 'pfam':
							if evalue < best_hits_by_evalue[query][1]:
								best_hits_by_evalue[query] = [[hit], evalue]
							elif evalue == best_hits_by_evalue[query][1]:
								best_hits_by_evalue[query][0].append(hit)
						else:
							if evalue < min_annotation_evalue:
								best_hits_by_evalue[query][0].append(hit)
			else:
				# parse DIAMOND BLASTp based results
				with open(annotation_dir + rf) as oarf:
					for line in oarf:
						line = line.strip()
						ls = line.split('\t')
						query = ls[0]
						hit = ls[1]
						evalue = float(ls[-2])
						if evalue < best_hits_by_evalue[query][1]:
							best_hits_by_evalue[query] = [[hit], evalue]
						elif evalue == best_hits_by_evalue[query][1]:
							best_hits_by_evalue[query][0].append(hit)
			with open(protein_faa) as opf:
				for rec in SeqIO.parse(opf, 'fasta'):
					if rec.id in best_hits_by_evalue:
						annotations[db_name][rec.id] = [[id_to_description[x] for x in best_hits_by_evalue[rec.id][0]], best_hits_by_evalue[rec.id][1]]
		return(default_to_regular(annotations))
	except Exception as e:
		sys.stderr.write('Issues with creating profile HMMs and consensus sequences.\n')
		logObject.error('Issues with creating profile HMMs and consensus sequences.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def default_to_regular(d):
	"""
	Function taken from
	https://stackoverflow.com/questions/26496831/how-to-convert-defaultdict-of-defaultdicts-of-defaultdicts-to-dict-of-dicts-o
	"""
	try:
		if isinstance(d, defaultdict):
			d = {k: default_to_regular(v) for k, v in d.items()}
		return d
	except:
		sys.exit(1)

def determineConsensusOrderOfHGs(genbanks, ortho_matrix_file, logObject):
	try:
		gc_gene_to_hg = {}
		core_hgs = set([])
		with open(ortho_matrix_file) as omf:
			for i, line in enumerate(omf):
				if i == 0: continue
				line = line.rstrip('\n')
				ls = line.split('\t')
				hg = ls[0]
				sample_count = 0
				for lts in ls[1:]:
					for lt in lts.split(', '):
						gc_gene_to_hg[lt] = hg
					if lts.strip() != '':
						sample_count += 1
				if sample_count/float(len(ls[1:])) == 1.0:
					core_hgs.add(hg)

		gc_gene_counts = defaultdict(int)
		gc_genes = defaultdict(set)
		gc_gene_locations = {}
		for gbk in genbanks:
			prefix = '.'.join(gbk.split('/')[-1].split('.')[:-1])
			gene_locations = util.parseGbk(gbk, prefix, logObject)
			for g in gene_locations:
				gc_gene_locations[g] = gene_locations[g]
			gc_gene_counts[gbk] = len(gc_gene_locations)
			gc_genes[gbk] = set(gc_gene_locations.keys())

		ref_bgc = None
		for i, item in enumerate(sorted(gc_gene_counts.items(), key=itemgetter(1), reverse=True)):
			ref_bgc = item[0]
			break

		gcs_ordered = [ref_bgc] + sorted(list(set(gc_genes.keys()).difference(set([ref_bgc]))))
		ref_hg_directions = {}
		hg_pair_scpus = defaultdict(int)
		hg_preceding_scpus = defaultdict(lambda: defaultdict(int))
		hg_following_scpus = defaultdict(lambda: defaultdict(int))
		all_hgs = set(['start', 'end'])
		direction_forward_support = defaultdict(int)
		direction_reverse_support = defaultdict(int)
		for i, gc in enumerate(gcs_ordered):
			curr_gc_genes = gc_genes[gc]
			hg_directions = {}
			hg_lengths = defaultdict(list)
			hg_starts = {}
			for g in sorted(curr_gc_genes):
				ginfo = gc_gene_locations[g]
				gstart = ginfo['start']
				gend = ginfo['end']
				if g in gc_gene_to_hg:
					hg = gc_gene_to_hg[g]
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
					hg_after = hgs[j + 1]
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
			if hps[0][0] in core_hgs or hps[0][1] in core_hgs:
				anchor_edge = hps[0]
				break
		try:
			assert (anchor_edge != None)
		except:
			sys.stderr.write("Unexpected error, no anchor edge found, could be because no protocore homolog group exists, which shouldn't be the case!\n")
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
					# print(hg + '\t' + str(best_score) + '\t'+ relative_pos + '\t' + str(neighboriest_hg) + '\t' + str(neighboriest_hg_index))

					if relative_pos == 'before':
						ordered_hgs_list.insert(neighboriest_hg_index, hg)
					elif relative_pos == 'after':
						ordered_hgs_list.insert(neighboriest_hg_index + 1, hg)
					accounted_hgs.add(hg)
					not_accounted_hgs = all_hgs.difference(accounted_hgs)
					progress_made = True
					break

			if not progress_made:
				break
		# these shouldn't really exist but just append them to the end if they do
		unaccountable_hgs = all_hgs.difference(accounted_hgs)
		ordered_hgs_list += list(sorted(unaccountable_hgs))

		hg_order_scores = {}
		i = 1
		for hg in ordered_hgs_list:
			if not hg in set(['start', 'end']):
				consensus_direction = '+'
				if direction_forward_support[hg] >= direction_reverse_support[hg]: consensus_direction = '-'
				hg_order_scores[hg] = [i, consensus_direction]
				i += 1
		return hg_order_scores
	except Exception as e:
		sys.stderr.write('Issues in attempting to calculate order score for each homolog group.\n')
		logObject.error("Issues in attempting to calculate order score for each homolog group.")
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def determineHGStats(orthogroup_matrix_file, hg_nucl_dir, logObject):
	try:
		hg_single_copy_status = {}
		hg_prop_samples = {}
		hg_lts = defaultdict(set)
		with open(orthogroup_matrix_file) as omf:
			for i, line in enumerate(omf):
				line = line.rstrip('\n')
				ls = line.split('\t')
				if i == 0: continue
				hg = ls[0]
				is_single_copy = True
				sample_count = 0
				for lts in ls[1:]:
					if ',' in lts:
						is_single_copy = False
					if lts.strip() != '':
						sample_count += 1
					for lt in lts.split(', '):
						if lt.strip() == '': continue
						hg_lts[hg].add(lt)
				hg_single_copy_status[hg] = is_single_copy
				hg_prop_samples[hg] = sample_count/float(len(ls[1:]))

		hg_median_lengths = {}
		for f in os.listdir(hg_nucl_dir):
			hg = f.split('.fna')[0]
			lengths = []
			with open(hg_nucl_dir + f) as ohpf:
				for rec in SeqIO.parse(ohpf, 'fasta'):
					lengths.append(len(str(rec.seq)))
			hg_median_lengths[hg] = statistics.mean(lengths)
		return([hg_single_copy_status, hg_prop_samples, hg_median_lengths, dict(hg_lts)])
	except Exception as e:
		logObject.error('Issues with determining basic stats for homolog groups.')
		sys.stderr.write('Issues with determining basic stats for homolog groups.\n')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def individualHyphyRun(inputs):
	hg, hg_codo_algn_file, hg_full_codo_tree_file, gard_output, best_gard_output, fubar_outdir, skip_gard, gard_mode, logObject = inputs
	try:
		input_gbks_with_hg = set([])
		with open(hg_codo_algn_file) as ohcaf:
			for rec in SeqIO.parse(ohcaf, 'fasta'):
				input_gbks_with_hg.add(rec.id.split("|")[0])

		assert(len(input_gbks_with_hg) > 0)
		if len(input_gbks_with_hg) < 4:
			return

		if skip_gard:
			fubar_cmd = ['hyphy', 'CPU=1', 'fubar', '--alignment', hg_codo_algn_file, '--tree', hg_full_codo_tree_file]
			try:
				subprocess.call(' '.join(fubar_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
								executable='/bin/bash')
				assert(os.path.isfile(hg_codo_algn_file + '.FUBAR.json'))
				os.system('mv %s %s' % (hg_codo_algn_file + '.FUBAR.json', fubar_outdir))
				logObject.info('Successfully ran: %s' % ' '.join(fubar_cmd))
			except Exception as e:
				logObject.error('Had an issue running FUBAR: %s' % ' '.join(fubar_cmd))
				sys.stderr.write('Had an issue running FUBAR: %s\n' % ' '.join(fubar_cmd))
				logObject.error(e)
				sys.exit(1)
		else:
			gard_cmd = ['hyphy', 'CPU=1', 'gard', '--mode', gard_mode, '--alignment', hg_codo_algn_file,
							  '--output', gard_output, '--output-lf', best_gard_output]

			try:
				subprocess.call(' '.join(gard_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
								executable='/bin/bash')
				assert(os.path.isfile(best_gard_output))
				logObject.info('Successfully ran: %s' % ' '.join(gard_cmd))
			except Exception as e:
				logObject.error('Had an issue running FUBAR: %s' % ' '.join(gard_cmd))
				sys.stderr.write('Had an issue running FUBAR: %s\n' % ' '.join(gard_cmd))
				logObject.error(e)
				sys.exit(1)

			fubar_cmd = ['hyphy', 'CPU=1', 'fubar', '--alignment', best_gard_output]
			print(' '.join(fubar_cmd))
			try:
				subprocess.call(' '.join(fubar_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
								executable='/bin/bash')
				assert(os.path.isfile(best_gard_output + '.FUBAR.json'))
				os.system('mv %s %s' % (best_gard_output + '.FUBAR.json', fubar_outdir))
				logObject.info('Successfully ran: %s' % ' '.join(fubar_cmd))
			except Exception as e:
				logObject.error('Had an issue running FUBAR: %s' % ' '.join(fubar_cmd))
				sys.stderr.write('Had an issue running FUBAR: %s\n' % ' '.join(fubar_cmd))
				logObject.error(e)
				sys.exit(1)

	except Exception as e:
		sys.stderr.write('Issues with running HYPHY based analyses for homolog group %s\n' % hg)
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)
		sys.stderr.write(traceback.format_exc())


def runHyphyAnalyses(codo_algn_dir, tree_dir, gard_results_dir, fubar_results_dir, logObject, cpus=1, skip_gard=False, gard_mode="Faster"):
	try:
		hyphy_inputs = []
		for caf in os.listdir(codo_algn_dir):
			if not caf.endswith('.msa.fna'): continue
			hg = caf.split('.msa.fna')[0]
			hg_codo_algn_file = codo_algn_dir + caf
			hg_full_codo_tree_file = tree_dir + hg + '.tre'
			gard_output = gard_results_dir + hg + '.json'
			best_gard_output = gard_results_dir + hg + '.best'
			fubar_output = fubar_results_dir + hg + '.json'

			hyphy_inputs.append([hg, hg_codo_algn_file, hg_full_codo_tree_file, gard_output, best_gard_output, fubar_results_dir,
							skip_gard, gard_mode, logObject])

		p = multiprocessing.Pool(cpus)
		p.map(individualHyphyRun, hyphy_inputs)
		p.close()

		gard_partitions = {}
		for f in os.listdir(gard_results_dir):
			if f.endswith('.json'):
				hg = f.split('.json')[0]
				gard_json_result = gard_results_dir + f
				with open(gard_json_result) as ogjr:
					gard_results = json.load(ogjr)
				number_of_partitions = len(gard_results['trees'])
				gard_partitions[hg] = number_of_partitions

		fubar_sel_props = {}
		fubar_sel_sites = {}
		for f in os.listdir(fubar_results_dir):
			if f.endswith('.json'):
				hg = f.split('.msa.fna.FUBAR.json')[0]
				fubar_json_result = fubar_results_dir + f
				with open(fubar_json_result) as ofjr:
					fubar_results = json.load(ofjr)
				pos_selected_sites = 0
				neg_selected_sites = 0
				for partition in fubar_results['MLE']['content']:
					for site_mle_info in fubar_results['MLE']['content'][partition]:
						alpha, beta, diff, prob_agb, prob_alb, bayesfactor, _, _ = site_mle_info
						if prob_agb >= 0.9:
							neg_selected_sites += 1
						if prob_alb >= 0.9:
							pos_selected_sites += 1
				tot_selected_sites = pos_selected_sites + neg_selected_sites
				prop_selected_sites_positive = 'NA'
				if tot_selected_sites >= 1:
					prop_selected_sites_positive = float(pos_selected_sites)/float(neg_selected_sites+pos_selected_sites)
				fubar_sel_props[hg] = prop_selected_sites_positive
				fubar_sel_sites[hg] = tot_selected_sites
		return([gard_partitions, fubar_sel_props, fubar_sel_sites])
	except Exception as e:
		sys.stderr.write('Issues with running GARD or FUBAR analyses.\n')
		logObject.error('Issues with running GARD or FUBAR analyses.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def runGeneTreeCongruenceAnalysis(genbanks, tree_dir, gtc_results_dir, logObject, subsample_max=1000):
	try:
		# Run STAG by David Emms and Steve Kelley
		stag_map_file = gtc_results_dir + 'Species_Mapping.txt'
		smf_handle = open(stag_map_file, 'w')
		prefices = []
		for gbk in genbanks:
			prefix = '.'.join(gbk.split('/')[-1].split('.')[:-1])
			prefices.append(prefix)
			smf_handle.write(prefix + '|*\t' + prefix + '\n')
		smf_handle.close()

		stag_cmd = [stag_prog, stag_map_file, tree_dir]
		result_gc_tree = None
		try:
			subprocess.call(' '.join(stag_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			par_dir_with_species_tree = '/'.join(os.path.abspath(tree_dir).split('/')[:-1]) + '/'
			for sd in os.listdir(par_dir_with_species_tree):
				if sd.startswith("STAG"):
					result_gc_tree = par_dir_with_species_tree + sd + '/SpeciesTree.tre'
			assert(os.path.isfile(result_gc_tree))
			logObject.info('Successfully ran: %s' % ' '.join(stag_cmd))
		except Exception as e:
			logObject.error('Had an issue running STAG: %s' % ' '.join(stag_cmd))
			sys.stderr.write('Had an issue running STAG: %s\n' % ' '.join(stag_cmd))
			logObject.error(e)
			sys.exit(1)

		all_pairwise_combinations = itertools.combinations(sorted(prefices), 2)
		all_pairwise_combinations = sorted(list(all_pairwise_combinations))
		num_total_pairwise = len(all_pairwise_combinations)
		subsample_size = min(subsample_max, num_total_pairwise)
		selected_pairs = random.sample(all_pairwise_combinations, k=subsample_size)

		gc_pw_info = util.calculateSelectDistances(result_gc_tree, selected_pairs)
		hg_to_gc_congruence = {}
		for hgt in os.listdir(tree_dir):
			hg = hgt.split('.tre')[0]
			outf = gtc_results_dir + hg + '.txt'
			util.computeCongruence(hg, tree_dir + hgt, gc_pw_info, selected_pairs, outf, logObject)

		for f in os.listdir(gtc_results_dir):
			if f == 'Species_Mapping.txt': continue
			hg = f.split('.txt')[0]
			with open(gtc_results_dir + f) as ogf:
				for line in ogf:
					line = line.strip()
					hg_to_gc_congruence[hg] = line
		return([result_gc_tree, hg_to_gc_congruence])
	except Exception as e:
		sys.stderr.write('Issues with running gene tree to gene-cluster tree congruence analysis.\n')
		logObject.error('Issues with running gene tree to gene-cluster tree congruence analysis.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def runTajimasDAnalysis(codo_algn_trim_dir, logObject):
	try:
		hg_tajimas_d = {}
		for catf in os.listdir(codo_algn_trim_dir):
			if not catf.endswith('.msa.fna'): continue
			hg = catf.split('.msa.fna')[0]
			codo_algn_trimmed_file = codo_algn_trim_dir + catf
			codo_sequences = []
			with open(codo_algn_trimmed_file) as ocatf:
				for rec in SeqIO.parse(ocatf, 'fasta'):
					codo_sequences.append(str(rec.seq))
			# alignment must be > 100 bp and have 4 or more sequences
			if len(codo_sequences) >= 4 and len(codo_sequences[0]) >= 100:
				taj_d = calculateTajimasD(codo_sequences)
			else:
				taj_d = 'NA'
			hg_tajimas_d[hg] = taj_d
		return(hg_tajimas_d)
	except Exception as e:
		sys.stderr.write('Issues with calculating Tajima\'s D for homolog groups.\n')
		logObject.error('Issues with calculating Tajima\s D for homolog groups.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

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

def consolidateReport(hg_stats, annotations, evo_stats, final_report_xlsx, final_report_tsv, logObject):
	"""
	dict_keys(['pfam', 'vfdb', 'paperblast', 'pgap', 'vog', 'isfinder', 'card', 'mibig'])
	dict_keys(['hg_single_copy_status', 'hg_prop_samples', 'hg_median_lengths', 'hg_order_scores'])
	dict_keys(['tajimas_d', 'gard_partitions', 'fubar_sel_props', 'fubar_sel_sites', 'gene_tree_congruence'])
	"""
	try:
		header = ['Homolog Group (HG) ID', 'HG is Single Copy?', 'Proportion of GenBanks with HG',
				  'HG Median Length (bp)', 'HG Consensus Order', 'HG Consensus Direction', 'Tajima\'s D',
				  'GARD Partitions Based on Recombination Breakpoints',
				  'Number of Sites Identified as Under Positive or Negative Selection by FUBAR',
				  'Proportion of Sites Under Selection which are Positive', 'HG Tree to Gene Cluster Tree Congruence',
				  'KO Annotation (E-value)', 'PGAP Annotation (E-value)', 'PaperBLAST Annotation (E-value)',
				  'CARD Annotation (E-value)', 'IS Finder (E-value)', 'MI-BiG Annotation (E-value)',
				  'VOG Annotation (E-value)',  'VFDB Annotation (E-value)', 'Pfam Domains', 'CDS Locus Tags']

		frt_handle = open(final_report_tsv, 'w')
		frt_handle.write('\t'.join(header) + '\n')
		# traverse HG in consensus order
		for hg_tup in sorted(hg_stats['hg_order_scores'].items(), key=lambda e: e[1][0]):
			hg = hg_tup[0]
			hg_scs = hg_stats['hg_single_copy_status'][hg]
			hg_cons = hg_stats['hg_prop_samples'][hg]
			hg_mlen = hg_stats['hg_median_lengths'][hg]
			hg_lts = '; '.join(hg_stats['hg_locus_tags'][hg])
			hg_ordr = hg_tup[1][0]
			hg_dire = hg_tup[1][1]
			hg_tajd = util.gatherValueFromDictForHomologGroup(hg, evo_stats['tajimas_d'])
			hg_gpar = util.gatherValueFromDictForHomologGroup(hg, evo_stats['gard_partitions'])
			hg_ssit = util.gatherValueFromDictForHomologGroup(hg, evo_stats['fubar_sel_sites'])
			hg_spro = util.gatherValueFromDictForHomologGroup(hg, evo_stats['fubar_sel_props'])
			hg_gtcc = util.gatherValueFromDictForHomologGroup(hg, evo_stats['gene_tree_congruence'])
			ko_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, annotations['ko'])
			pgap_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, annotations['pgap'])
			pb_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, annotations['paperblast'])
			card_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, annotations['card'])
			isf_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, annotations['isfinder'])
			mibig_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, annotations['mibig'])
			vog_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, annotations['vog'])
			vfdb_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, annotations['vfdb'])
			pfam_annots = 'NA'
			if hg in annotations['pfam']:
				pfam_annots = '; '.join(annotations['pfam'][hg][0])

			row = [hg, hg_scs, hg_cons, hg_mlen, hg_ordr, hg_dire, hg_tajd, hg_gpar, hg_ssit, hg_spro, hg_gtcc,
				   ko_annot, pgap_annot, pb_annot, card_annot, isf_annot, mibig_annot, vog_annot, vfdb_annot,
				   pfam_annots, hg_lts]
			row = [str(x) for x in row]
			frt_handle.write('\t'.join(row) + '\n')
		frt_handle.close()

		# Generate Excel spreadsheet
		writer = pd.ExcelWriter(final_report_xlsx, engine='xlsxwriter')
		workbook = writer.book
		dd_sheet = workbook.add_worksheet('Data Dictionary')
		dd_sheet.write(0, 0, 'WARNING: some evolutionary statistics are experimental - evaluate with caution!!!')
		dd_sheet.write(1, 0, 'Data Dictionary describing columns of "Overview" spreadsheets can be found on zol\'s Wiki at:')
		dd_sheet.write(2, 0, 'https://github.com/Kalan-Lab/zol/wiki/XXXXX')

		numeric_columns = set(['Proportion of GenBanks with HG', 'HG Median Length (bp)', 'HG Consensus Order',
							   'Tajima\'s D', 'GARD Partitions Based on Recombination Breakpoints',
							    'Number of Sites Identified as Under Positive or Negative Selection by FUBAR',
		'Proportion of Sites Under Selection which are Positive", "HG Tree to Gene Cluster Tree Congruence'])
		results_df = util.loadTableInPandaDataFrame(final_report_tsv, numeric_columns)
		results_df.to_excel(writer, sheet_name='ZoL Results', index=False, na_rep="NA")
		workbook.close()

	except Exception as e:
		sys.stderr.write('Issues creating consolidated results files.\n')
		logObject.error('Issues creating consolidated results files.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def runTreemmer(gene_cluster_tree, max_for_visualization, logObject):
	representative_genbanks = set([])
	try:
		retain_listing_file = gene_cluster_tree + '_trimmed_list_X_' + str(max_for_visualization)
		treemmer_cmd = ['Treemmer_v0.3.py', '-X=' + str(max_for_visualization), gene_cluster_tree]
		try:
			subprocess.call(' '.join(treemmer_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(retain_listing_file))
			logObject.info('Successfully ran: %s' % ' '.join(treemmer_cmd))
		except Exception as e:
			logObject.error('Had an issue running Treemmer: %s' % ' '.join(treemmer_cmd))
			sys.stderr.write('Had an issue running Treemmer: %s\n' % ' '.join(treemmer_cmd))
			logObject.error(e)
			sys.exit(1)
		with open(retain_listing_file) as orlf:
			for line in orlf:
				line = line.strip()
				representative_genbanks.add(line)
	except Exception as e:
		sys.stderr.write('Issues running Treemmer to reduce the set of input GenBanks to a representative set for easier visualization with clinker.\n')
		logObject.error('Issues running Treemmer to reduce the set of input GenBanks to a representative set for easier visualization with clinker.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)
	return(representative_genbanks)

def runClinker(genbanks, hg_lts, representative_genbanks, fin_dir, work_dir, logObject):
	try:
		renamed_lts_dir = work_dir + 'GenBanks_LTs_Renamed/'
		input_gbks = []
		for gbk in genbanks:
			gbk_name = '.'.join(gbk.split('/')[-1].split('.')[:-1])
			if not gbk_name in representative_genbanks: continue
			updated_gbk = renamed_lts_dir + gbk
			ug_handle = open(updated_gbk, 'w')
			with open(gbk) as ogbk:
				for rec in SeqIO.parse(ogbk, 'fasta'):
					for feature in rec.features:
						if not feature.type == 'CDS': continue
						lt = feature.qualifiers.get('locus_tag')[0]
						feature.qualifiers.get('locus_tag')[0] = gbk_name + '|' + lt
					SeqIO.write(rec, ug_handle, 'genbank')
			ug_handle.close()
			input_gbks.append(updated_gbk)

		hg_mapping_file = work_dir + 'Locus_Tag_to_Homolog_Group_Mapping.csv'
		hmf_handle = open(hg_mapping_file, 'w')
		for hg in hg_lts:
			for lt in hg_lts[hg]:
				if lt.split('|')[0] in representative_genbanks:
					hmf_handle.write(lt + ',' + hg + '\n')
		hmf_handle.close()

		result_html = fin_dir + 'clinker_Visual_of_Representative_Gene_Clusters.html'
		clinker_cmd = ['clinker', '-gf', hg_mapping_file, '-p', result_html] + input_gbks
		try:
			subprocess.call(' '.join(clinker_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(result_html))
			logObject.info('Successfully ran: %s' % ' '.join(clinker_cmd))
		except Exception as e:
			logObject.error('Had an issue running clinker: %s' % ' '.join(clinker_cmd))
			sys.stderr.write('Had an issue running clinker: %s\n' % ' '.join(clinker_cmd))
			logObject.error(e)
			sys.exit(1)
	except Exception as e:
		sys.stderr.write('Issues running clinker or creating inputs for clinker.\n')
		logObject.error('Issues running clinker or creating inputs for clinker.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)
