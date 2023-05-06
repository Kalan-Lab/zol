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
import decimal
import pickle
import shutil
from scipy import stats

zol_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
plot_prog = zol_main_directory + 'zol/clusterHeatmap.R'

def reinflateOrthoGroups(ortho_matrix_file, prot_dir, rog_dir, logObject, cpus=1):
	try:
		reps = set([])
		rep_to_hg = {}
		with open(ortho_matrix_file) as oomf:
			for i, line in enumerate(oomf):
				if i == 0: continue
				line = line.strip('\n')
				ls = line.split('\t')
				hg = ls[0]
				for pids in ls[1:]:
					for pid in pids.split(','):
						pid = pid.strip()
						reps.add(pid)
						rep_to_hg[pid] = hg

		comp_prot_file = rog_dir + 'All_Proteins.faa'
		os.system('find %s -type f -name "*.faa" -exec cat {} + >> %s' % (prot_dir, comp_prot_file))

		cdhit_nr_prefix = rog_dir + 'CD-HIT_Results'
		cdhit_cluster_file = cdhit_nr_prefix + '.clstr'
		cdhit_cmd = ['cd-hit', '-i', comp_prot_file, '-o', cdhit_nr_prefix, '-c', '0.98', '-aL', '0.95', '-aS', '0.95',
					 '-d', '0', '-T', str(cpus), '-M', '20000']

		try:
			subprocess.call(' '.join(cdhit_cmd), shell=True, stdout=subprocess.DEVNULL,
							stderr=subprocess.DEVNULL, executable='/bin/bash')
			assert(os.path.isfile(cdhit_nr_prefix))
			assert(os.path.isfile(cdhit_cluster_file))
		except Exception as e:
			logObject.error("Issue with running: %s" % ' '.join(cdhit_cmd))
			logObject.error(e)
			raise RuntimeError(e)

		cluster_rep = None
		clust_proteins = defaultdict(set)
		protein_to_clust = {}
		tmp = []
		with open(cdhit_cluster_file) as occf:
			for line in occf:
				line = line.strip()
				if line.startswith('>'):
					if len(tmp) > 0 and cluster_rep != None:
						for lt in tmp:
							protein_to_clust[lt] = cluster_rep
							clust_proteins[cluster_rep].add(lt)
					tmp = []
					cluster_rep = None
					continue
				ls = line.split()
				lt = ls[2][1:-3]
				if line.endswith(' *'):
					cluster_rep = lt
				tmp.append(lt)
		if len(tmp) > 0 and cluster_rep != None:
			for lt in tmp:
				protein_to_clust[lt] = cluster_rep
				clust_proteins[cluster_rep].add(lt)

		all_samples = set([])
		for f in os.listdir(prot_dir):
			all_samples.add('.faa'.join(f.split('.faa')[:-1]))

		inflated_og_matrix_file = rog_dir + 'Orthogroups.tsv'
		iomf_handle = open(inflated_og_matrix_file, 'w')
		accounted = set([])
		with open(ortho_matrix_file) as oomf:
			for i, line in enumerate(oomf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					iomf_handle.write('Sample\t' + '\t'.join(sorted(list(all_samples))) + '\n')
					continue
				hg = ls[0]
				cluster_obs = set([])
				all_pids_by_sample = defaultdict(set)
				for pids in ls[1:]:
					for pid in pids.split(','):
						pid = pid.strip()
						if pid == '': continue
						cluster_obs.add(protein_to_clust[pid])
						all_pids_by_sample[pid.split('|')[0]].add(pid)
				for clust in cluster_obs:
					for pid in clust_proteins[clust]:
						if pid in reps and rep_to_hg[pid] != hg:
							sys.stderr.write('Warning: The protein %s is a representative of homolog group %s, but can potentially belong to multiple. Skipping its incorporation for homolog group %s.\n' % (pid, rep_to_hg[pid], hg))
							logObject.warning('The protein %s is a representative of a homolog group %s, but can potentially belong to multiple. Skipping its incorporation for homolog group %s.' % (pid, rep_to_hg[pid], hg))
						elif pid in accounted:
							sys.stderr.write('Warning: The protein %s has already been clustered into a homolog group, but can potentially belong to multiple. Skipping its incorporation for homolog group %s.\n' % (pid, hg))
							logObject.warning('The protein %s has already been clustered into a homolog group, but can potentially belong to multiple. Skipping its incorporation for homolog group %s.' % (pid, hg))
						else:
							all_pids_by_sample[pid.split('|')[0]].add(pid)
						accounted.add(pid)
				row = [hg]
				for sample in sorted(list(all_samples)):
					samp_pids_for_hg = ''
					if sample in all_pids_by_sample:
						samp_pids_for_hg = ', '.join(sorted(all_pids_by_sample[sample]))
					row.append(samp_pids_for_hg)
				iomf_handle.write('\t'.join(row) + '\n')
		iomf_handle.close()

	except Exception as e:
		sys.stderr.write('Issues with reinflation of orthogroups to full gene-cluster set.\n')
		logObject.error('Issues with reinflation of orthogroups to full gene-cluster set.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def dereplicateUsingSkani(genbanks, focal_genbanks, derep_dir, kept_dir, logObject, skani_identiy_threshold=99.0, skani_coverage_threshold=95.0, mcl_inflation=None, cpus=1):
	derep_genbanks = set([])
	try:
		full_nucl_seq_dir = derep_dir + 'FASTAs/'
		util.setupReadyDirectory([full_nucl_seq_dir])
		fasta_listing_file = derep_dir + 'Gene_Clusters_FASTA_Listing.txt'
		flf_handle = open(fasta_listing_file, 'w')
		longest_seq = defaultdict(int)
		tot_seq = defaultdict(int)
		for gbk in genbanks:
			gbk_prefix = None
			if gbk.endswith('.gbk') or gbk.endswith('.gbff') or gbk.endswith('.genbank'):
				gbk_prefix = '.'.join(gbk.split('/')[-1].split('.')[:-1])
			assert(gbk_prefix != None)
			gbk_fasta_file = full_nucl_seq_dir + gbk_prefix + '.fasta'
			flf_handle.write(gbk_fasta_file + '\n')
			gbk_fasta_handle = open(gbk_fasta_file, 'w')
			with open(gbk) as ogbk:
				for rec in SeqIO.parse(ogbk, 'genbank'):
					gbk_fasta_handle.write('>' + rec.id + '\n' + str(rec.seq) + '\n')
					if longest_seq[gbk_prefix] < len(str(rec.seq)):
						longest_seq[gbk_prefix] = len(str(rec.seq))
					tot_seq[gbk_prefix] += len(str(rec.seq))
			gbk_fasta_handle.close()
		flf_handle.close()

		skani_sketch_db = derep_dir + 'skani_sketch/'
		skani_sketch_cmd = ['skani', 'sketch', '-t', str(cpus), '-l', fasta_listing_file, '-o', skani_sketch_db]
		try:
			subprocess.call(' '.join(skani_sketch_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isdir(skani_sketch_db))
			logObject.info('Successfully ran skani: %s' % ' '.join(skani_sketch_cmd))
		except Exception as e:
			logObject.error('Had an issue running skani: %s' % ' '.join(skani_sketch_cmd))
			sys.stderr.write('Had an issue running skani: %s\n' % ' '.join(skani_sketch_cmd))
			logObject.error(e)
			sys.exit(1)

		skani_result_file = derep_dir + 'skani_results.tsv'
		skani_dist_cmd = ['skani', 'dist', '-t', str(cpus), '-q', skani_sketch_db + '*', '-r', skani_sketch_db + '*',
						  '--min-af', str(skani_coverage_threshold), '-o', skani_result_file]

		try:
			subprocess.call(' '.join(skani_dist_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(skani_result_file))
			logObject.info('Successfully ran skani: %s' % ' '.join(skani_dist_cmd))
		except Exception as e:
			logObject.error('Had an issue running skani: %s' % ' '.join(skani_dist_cmd))
			sys.stderr.write('Had an issue running skani: %s\n' % ' '.join(skani_dist_cmd))
			logObject.error(e)
			sys.exit(1)

		similar_pairs_file = derep_dir + 'Similar_Pairs.txt'
		similar_pairs_handle = open(similar_pairs_file, 'w')
		all_gcs = set([])
		paired_gcs = set([])
		visited = set([])
		with open(skani_result_file) as osrf:
			for i, line in enumerate(osrf):
				if i == 0: continue
				line = line.strip()
				f1, f2, ani, _, _, _, _ = line.split('\t')
				s1 = '.fasta'.join(f1.split('/')[-1].split('.fasta')[:-1])
				s2 = '.fasta'.join(f2.split('/')[-1].split('.fasta')[:-1])
				if float(ani) >= skani_identiy_threshold:
					pair_tup = sorted([s1, s2])
					if not tuple(pair_tup) in visited:
						if mcl_inflation == None:
							similar_pairs_handle.write(pair_tup[0] + '\t' + pair_tup[1] + '\n')
						else:
							similar_pairs_handle.write(pair_tup[0] + '\t' + pair_tup[1] + '\t' + str(ani) + '\n')
						if s1 != s2:
							paired_gcs.add(s1)
							paired_gcs.add(s2)
						else:
							all_gcs.add(s1)
					visited.add(tuple(pair_tup))
		similar_pairs_handle.close()

		focal = None
		if focal_genbanks != None and os.path.isfile(focal_genbanks):
			focal = set([])
			with open(focal_genbanks) as ofgf:
				for line in ofgf:
					gbk = line.strip()
					gbk_prefix = None
					if gbk.endswith('.gbk') or gbk.endswith('.gbff') or gbk.endswith('.genbank'):
						gbk_prefix = '.'.join(gbk.split('/')[-1].split('.')[:-1])
					assert(gbk_prefix != None)
					focal.add(gbk_prefix)

		clusters_file = derep_dir + 'Cluster_Families.txt'
		if mcl_inflation == None:
			clust_cmd = ['slclust', '<', similar_pairs_file, '>', clusters_file]
		else:
			clust_cmd = ['mcl', similar_pairs_file, '--abc', '-I', str(mcl_inflation), '-o', clusters_file,
						 '-te', str(cpus)]

		try:
			subprocess.call(' '.join(clust_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
			assert (os.path.isfile(clusters_file))
			logObject.info('Successfully ran: %s' % ' '.join(clust_cmd))
		except Exception as e:
			logObject.error('Had an issue running slclust or mcl: %s' % ' '.join(clust_cmd))
			sys.stderr.write('Had an issue running slclust or mcl: %s\n' % ' '.join(clust_cmd))
			logObject.error(e)
			sys.exit(1)

		representatives = set([])
		rep_genbank_members = defaultdict(set)
		with open(clusters_file) as ocf:
			for line in ocf:
				line = line.strip()
				gcc = line.split()
				cluster_gc_stats = []
				gcs = set([])
				for gc in gcc:
					cluster_gc_stats.append([gc, longest_seq[gc], tot_seq[gc]])
					gcs.add(gc)
				focal_not_relevant = False
				if focal == None or len(gcs.intersection(focal)) == 0:
					focal_not_relevant = True
				rep = [x[0] for x in sorted(cluster_gc_stats, key=itemgetter(1,2), reverse=True) if focal_not_relevant or x[0] in focal][0]
				rep_genbank_members[rep] = gcs
				representatives.add(rep)

		for gc in all_gcs.difference(paired_gcs):
			rep_genbank_members[gc] = set([gc])
			representatives.add(gc)

		for gbk in genbanks:
			gbk_prefix = None
			if gbk.endswith('.gbk') or gbk.endswith('.gbff') or gbk.endswith('.genbank'):
				gbk_prefix = '.'.join(gbk.split('/')[-1].split('.')[:-1])
			assert (gbk_prefix != None)
			if not gbk_prefix in representatives: continue
			shutil.copy(gbk, kept_dir)
			derep_genbanks.add(kept_dir + gbk.split('/')[-1])

		num_gbk = len(derep_genbanks)
		if num_gbk == 0:
			sys.stderr.write('Issues with dereplicating GenBanks! All GenBanks deemed as redundant ...\n')
			logObject.error('Issues with dereplicating GenBanks! All GenBanks deemed as redundant ...')
		else:
			sys.stdout.write('Found %d GenBanks retained after dereplication.\n' % num_gbk)
			logObject.info('Found %d GenBanks retained after dereplication.' % num_gbk)

	except Exception as e:
		sys.stderr.write('Issues with run skani based dereplication of input GenBanks.\n')
		logObject.error('Issues with run skani based dereplication of input GenBanks.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return([derep_genbanks, rep_genbank_members])

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
					if not rec.id in g_to_hg: continue
					hg = g_to_hg[rec.id]
					hpf_handle = open(hg_prot_dir + hg + '.faa', 'a+')
					hpf_handle.write('>' + rec.description + '\n' + str(rec.seq) + '\n')
					hpf_handle.close()
		for nf in os.listdir(nucl_dir):
			nfile = nucl_dir + nf
			with open(nfile) as onf:
				for rec in SeqIO.parse(onf, 'fasta'):
					if not rec.id in g_to_hg: continue
					hg = g_to_hg[rec.id]
					hnf_handle = open(hg_nucl_dir + hg + '.fna', 'a+')
					hnf_handle.write('>' + rec.description + '\n' + str(rec.seq) + '\n')
					hnf_handle.close()
	except Exception as e:
		sys.stderr.write('Issues with partitioning sequences to homolog groups.\n')
		logObject.error('Issues with partitioning sequences to homolog groups.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)


def partitionAndCreateUpstreamNuclAlignments(ortho_matrix_file, nucl_upstr_dir, hg_upst_dir, upst_algn_dir, logObject, cpus=1, use_super5=False):
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

		for uf in os.listdir(nucl_upstr_dir):
			ufile = nucl_upstr_dir + uf
			with open(ufile) as ouf:
				for rec in SeqIO.parse(ouf, 'fasta'):
					if not rec.id in g_to_hg: continue
					hg = g_to_hg[rec.id]
					hpf_handle = open(hg_upst_dir + hg + '.fna', 'a+')
					hpf_handle.write('>' + rec.description + '\n' + str(rec.seq) + '\n')
					hpf_handle.close()

		for pf in os.listdir(hg_upst_dir):
			prefix = '.fna'.join(pf.split('.fna')[:-1])
			upst_file = hg_upst_dir + pf
			if os.path.getsize(upst_file) == 0: continue
			min_seq_len = 10000
			with open(upst_file) as ouf:
				for rec in SeqIO.parse(ouf, 'fasta'):
					if len(str(rec.seq)) < min_seq_len:
						min_seq_len = len(str(rec.seq))
			if min_seq_len < 10: continue
			upst_algn_file = upst_algn_dir + prefix + '.msa.fna'
			align_cmd = ['muscle', '-align', upst_file, '-output', upst_algn_file, '-nt', '-threads', str(cpus)]
			if use_super5:
				align_cmd = ['muscle', '-super5', upst_file, '-output', upst_algn_file, '-nt', '-threads', str(cpus)]
			try:
				subprocess.call(' '.join(align_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
								executable='/bin/bash')
				assert(os.path.isfile(upst_algn_file))
				logObject.info('Successfully ran: %s' % ' '.join(align_cmd))
			except Exception as e:
				logObject.error('Had an issue running MUSCLE: %s' % ' '.join(align_cmd))
				sys.stderr.write('Had an issue running MUSCLE: %s\n' % ' '.join(align_cmd))
				logObject.error(e)
				sys.exit(1)

	except Exception as e:
		sys.stderr.write('Issues with partitioning/aligning upstream sequences.\n')
		logObject.error('Issues with partitioning/aligning upstream sequences.')
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
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def refineGeneCalling(custom_database, hg_prot_dir, hg_nucl_dir, refine_workpace_dir, logObject, use_super5=True, cpus=1):
	hg_prot_refined_dir = refine_workpace_dir + 'HG_Protein_Refined_Sequences/'
	hg_nucl_refined_dir = refine_workpace_dir + 'HG_Nucleotide_Refined_Sequences/'
	try:
		tmp_hg_prot_seq_dir = refine_workpace_dir + 'HG_Prot_Seqs_with_References/'
		tmp_hg_prot_aln_dir = refine_workpace_dir + 'HG_Prot_Alns_with_References/'
		util.setupReadyDirectory([hg_prot_refined_dir, hg_nucl_refined_dir, tmp_hg_prot_seq_dir, tmp_hg_prot_aln_dir])

		concat_proteins_faa_file = refine_workpace_dir + 'All_Proteins.faa'
		cpff_handle = open(concat_proteins_faa_file, 'w')
		for f in os.listdir(hg_prot_dir):
			with open(hg_prot_dir + f) as ohf:
				for rec in SeqIO.parse(ohf, 'fasta'):
					cpff_handle.write('>' + f.split('.faa')[0] + '|' + rec.id + '\n' + str(rec.seq) + '\n')
		cpff_handle.close()

		dmnd_db = refine_workpace_dir + 'All_Proteins.dmnd'
		blastp_file = refine_workpace_dir + 'Blast_Results.txt'
		makedb_cmd = ['diamond', 'makedb', '--ignore-warnings', '--in', concat_proteins_faa_file, '-d', dmnd_db]
		search_cmd = ['diamond', 'blastp', '--ignore-warnings', '-p', str(cpus), '-d', dmnd_db, '-q', custom_database,
					  '-o', blastp_file]
		try:
			subprocess.call(' '.join(makedb_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(dmnd_db))
			logObject.info('Successfully ran: %s' % ' '.join(makedb_cmd))
		except Exception as e:
			logObject.error('Had an issue running DIAMOND makedb: %s' % ' '.join(makedb_cmd))
			sys.stderr.write('Had an issue running DIAMOND makedb: %s\n' % ' '.join(makedb_cmd))
			logObject.error(e)
			sys.exit(1)

		try:
			subprocess.call(' '.join(search_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(blastp_file))
			logObject.info('Successfully ran: %s' % ' '.join(search_cmd))
		except Exception as e:
			logObject.error('Had an issue running DIAMOND blastp: %s' % ' '.join(search_cmd))
			sys.stderr.write('Had an issue running DIAMOND blastp: %s\n' % ' '.join(search_cmd))
			logObject.error(e)
			sys.exit(1)

		best_query_hit = defaultdict(lambda: [set([]), 0.0])
		with open(blastp_file) as obf:
			for line in obf:
				line = line.strip()
				ls = line.split('\t')
				bitscore = float(ls[11])
				query = ls[0]
				hit_og = ls[1].split('|')[0]
				if best_query_hit[query][1] < bitscore:
					best_query_hit[query] = [set([hit_og]), bitscore]
				elif best_query_hit[query][1] == bitscore:
					best_query_hit[query][0].add(hit_og)

		cd_seqs = {}
		with open(custom_database) as ocd:
			for rec in SeqIO.parse(ocd, 'fasta'):
				cd_seqs[rec.id] = str(rec.seq)

		hg_matched = set([])
		for qb in sorted(best_query_hit.items(), key=lambda x: x[1][1], reverse=True):
			query = qb[0]
			assert(len(best_query_hit[query][0])==1)
			hg = list(best_query_hit[query][0])[0]
			if hg in hg_matched: continue
			hg_matched.add(hg)
			hg_prot_faa = hg_prot_dir + hg + '.faa'
			hg_nucl_fna = hg_nucl_dir + hg + '.fna'
			hg_seq_with_que = tmp_hg_prot_seq_dir + hg + '.faa'
			hg_aln_with_que = tmp_hg_prot_aln_dir + hg + '.msa.faa'

			hswq_handle = open(hg_seq_with_que, 'w')
			hswq_handle.write('>' + query + '\n' + cd_seqs[query] + '\n')
			with open(hg_prot_faa) as ohpf:
				for rec in SeqIO.parse(ohpf, 'fasta'):
					hswq_handle.write('>' + rec.id + '\n' + str(rec.seq) + '\n')
			hswq_handle.close()

			align_cmd = ['muscle', '-align', hg_seq_with_que, '-output', hg_aln_with_que, '-amino', '-threads', str(cpus)]
			if use_super5:
				align_cmd = ['muscle', '-super5', hg_seq_with_que, '-output', hg_aln_with_que, '-amino', '-threads', str(cpus)]
			try:
				subprocess.call(' '.join(align_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
								executable='/bin/bash')
				assert(os.path.isfile(hg_aln_with_que))
				logObject.info('Successfully ran: %s' % ' '.join(align_cmd))
			except Exception as e:
				logObject.error('Had an issue running MUSCLE: %s' % ' '.join(align_cmd))
				sys.stderr.write('Had an issue running MUSCLE: %s\n' % ' '.join(align_cmd))
				logObject.error(e)
				sys.exit(1)

			msa_avoid_pos = set([])
			with open(hg_aln_with_que) as ohawq:
				for rec in SeqIO.parse(ohawq, 'fasta'):
					if rec.id == query:
						for pos, bp in enumerate(str(rec.seq)):
							if bp == '-':
								msa_avoid_pos.add(pos)

			nucl_refine_seq_pos = defaultdict(set)
			hpr_handle = open(hg_prot_refined_dir + hg + '.faa', 'w')
			with open(hg_aln_with_que) as ohawq:
				for rec in SeqIO.parse(ohawq, 'fasta'):
					if rec.id != query:
						refined_prot_seq = ''
						ref_pos = 0
						for msa_pos, bp in enumerate(str(rec.seq)):
							if bp != '-' and not msa_pos in msa_avoid_pos:
								refined_prot_seq += bp
								nucl_refine_seq_pos[rec.id].add(ref_pos*3 + 0)
								nucl_refine_seq_pos[rec.id].add(ref_pos*3 + 1)
								nucl_refine_seq_pos[rec.id].add(ref_pos*3 + 2)
							if bp != '-':
								ref_pos += 1
						hpr_handle.write('>' + rec.id + '\n' + refined_prot_seq + '\n')
			hpr_handle.close()

			hnr_handle = open(hg_nucl_refined_dir + hg + '.fna', 'w')
			with open(hg_nucl_fna) as ohnf:
				for rec in SeqIO.parse(ohnf, 'fasta'):
					refined_nucl_seq = ''
					for pos, bp in enumerate(str(rec.seq)):
						if pos in nucl_refine_seq_pos[rec.id]:
							refined_nucl_seq += bp
					hnr_handle.write('>' + rec.id + '\n' + refined_nucl_seq + '\n')
			hnr_handle.close()

	except Exception as e:
		sys.stderr.write('Issues with refining gene-calling/exons using reference proteins.\n')
		logObject.error('Issues with refining gene-calling/exons using reference proteins.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

	return([hg_prot_refined_dir, hg_nucl_refined_dir])

def annotateCustomDatabase(protein_faa, custom_protein_db_faa, annotation_dir, logObject, cpus=1, max_annotation_evalue=1e-5):
	custom_annotations = {}
	try:
		custom_annot_dir = annotation_dir + 'Custom_Annotation/'
		util.setupReadyDirectory([custom_annot_dir])
		dmnd_db = custom_annot_dir + 'Custom.dmnd'
		blastp_file = custom_annot_dir + 'Custom.txt'
		makedb_cmd = ['diamond', 'makedb', '--ignore-warnings', '--in', custom_protein_db_faa, '-d', dmnd_db]
		search_cmd = ['diamond', 'blastp', '--ignore-warnings', '-p', str(cpus), '-d', dmnd_db, '-q', protein_faa,
					  '-o', blastp_file]

		try:
			subprocess.call(' '.join(makedb_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(dmnd_db))
			logObject.info('Successfully ran: %s' % ' '.join(makedb_cmd))
		except Exception as e:
			logObject.error('Had an issue running DIAMOND makedb: %s' % ' '.join(makedb_cmd))
			sys.stderr.write('Had an issue running DIAMOND makedb: %s\n' % ' '.join(makedb_cmd))
			logObject.error(e)
			sys.exit(1)

		try:
			subprocess.call(' '.join(search_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(blastp_file))
			logObject.info('Successfully ran: %s' % ' '.join(search_cmd))
		except Exception as e:
			logObject.error('Had an issue running DIAMOND blastp: %s' % ' '.join(search_cmd))
			sys.stderr.write('Had an issue running DIAMOND blastp: %s\n' % ' '.join(search_cmd))
			logObject.error(e)
			sys.exit(1)

		id_to_description = {}
		with open(custom_protein_db_faa) as ocpdf:
			for rec in SeqIO.parse(ocpdf, 'fasta'):
				id_to_description[rec.id] = rec.description

		best_hits_by_bitscore = defaultdict(lambda: [[], [], 0.0])
		with open(blastp_file) as obf:
			for line in obf:
				line = line.strip()
				ls = line.split('\t')
				que, hit = ls[:2]
				eval = float(ls[10])
				if eval > max_annotation_evalue: continue
				bitscore = float(ls[11])
				if bitscore > best_hits_by_bitscore[que][2]:
					best_hits_by_bitscore[que] = [[hit], [eval], bitscore]
				elif bitscore == best_hits_by_bitscore[que][2]:
					best_hits_by_bitscore[que][0].append(hit)
					best_hits_by_bitscore[que][1].append(eval)

		for que in best_hits_by_bitscore:
			custom_annotations[que] = [[id_to_description[x] for x in best_hits_by_bitscore[que][0]], best_hits_by_bitscore[que][1]]

	except Exception as e:
		sys.stderr.write('Issues with annotating using custom database.\n')
		logObject.error('Issues with annotating using custom database.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return(custom_annotations)

def annotateConsensusSequences(protein_faa, annotation_dir, logObject, cpus=1, max_annotation_evalue=1e-5):
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
					search_cmd = ['diamond', 'blastp', '--ignore-warnings', '-p', str(individual_cpus), '-d', db_file,
								  '-q', protein_faa, '-o', annotation_result_file, logObject]
				search_cmds.append(search_cmd)

		p = multiprocessing.Pool(pool_size)
		p.map(util.multiProcess, search_cmds)
		p.close()

		annotations = defaultdict(lambda: defaultdict(lambda: ['NA', 'NA'])) # db -> query -> [hit descriptions, evalue]
		for rf in os.listdir(annotation_dir):
			if not rf.endswith('.txt'): continue
			db_name = rf.split('.txt')[0]
			annot_info_file = name_to_info_file[db_name]

			id_to_description = defaultdict(lambda: 'NA')
			with open(annot_info_file) as oaif:
				for line in oaif:
					line = line.strip()
					ls = line.split('\t')
					id_to_description[ls[0]] = ls[1]

			# or by_score if HMMscan - lets avoid a second variable
			best_hits_by_bitscore = defaultdict(lambda: [[], [], 0.0])
			if db_name in hmm_based_annotations:
				# parse HMM based results from HMMER3
				with open(annotation_dir + rf) as oarf:
					for line in oarf:
						line = line.rstrip('\n')
						if line.startswith('#'): continue
						ls = line.split()
						query = ls[2]
						hit = ls[0]
						evalue = decimal.Decimal(ls[4])
						score = float(ls[5])
						if evalue > max_annotation_evalue: continue
						if db_name != 'pfam':
							if score > best_hits_by_bitscore[query][2]:
								best_hits_by_bitscore[query] = [[hit], [evalue], score]
							elif score == best_hits_by_bitscore[query][2]:
								best_hits_by_bitscore[query][0].append(hit)
								best_hits_by_bitscore[query][1].append(evalue)
						else:
							if evalue < max_annotation_evalue:
								best_hits_by_bitscore[query][0].append(hit)
			else:
				# parse DIAMOND BLASTp based results
				with open(annotation_dir + rf) as oarf:
					for line in oarf:
						line = line.strip()
						ls = line.split('\t')
						query = ls[0]
						hit = ls[1]
						bitscore = float(ls[11])
						evalue = decimal.Decimal(ls[10])
						if evalue > max_annotation_evalue: continue
						if bitscore > best_hits_by_bitscore[query][2]:
							best_hits_by_bitscore[query] = [[hit], [evalue], bitscore]
						elif bitscore == best_hits_by_bitscore[query][2]:
							best_hits_by_bitscore[query][0].append(hit)
							best_hits_by_bitscore[query][1].append(evalue)

			with open(protein_faa) as opf:
				for rec in SeqIO.parse(opf, 'fasta'):
					if rec.id in best_hits_by_bitscore:
						annotations[db_name][rec.id] = [[id_to_description[x] for x in best_hits_by_bitscore[rec.id][0]], best_hits_by_bitscore[rec.id][1]]
		return(default_to_regular(annotations))
	except Exception as e:
		sys.stderr.write('Issues with annotating consensus sequences.\n')
		logObject.error('Issues with annotating consensus sequences.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
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
		single_copy_core_hgs = set([])
		with open(ortho_matrix_file) as omf:
			for i, line in enumerate(omf):
				if i == 0: continue
				line = line.rstrip('\n')
				ls = line.split('\t')
				hg = ls[0]
				sample_count = 0
				sc_sample_count = 0
				for lts in ls[1:]:
					for lt in lts.split(', '):
						gc_gene_to_hg[lt] = hg
					if lts.strip() != '':
						sample_count += 1
					if lts.strip() != '' and not ', ' in lts:
						sc_sample_count += 1
				if sample_count/float(len(ls[1:])) == 1.0:
					core_hgs.add(hg)
					if sc_sample_count/float(len(ls[1:])) == 1.0:
						single_copy_core_hgs.add(hg)

		gc_gene_counts = defaultdict(int)
		gc_genes = defaultdict(set)
		scaff_relative_gene_locations = {}
		for gbk in genbanks:
			prefix = '.'.join(gbk.split('/')[-1].split('.')[:-1])
			gene_locations = util.parseGbk(gbk, prefix, logObject)
			scaffolds = set([gene_locations[x]['scaffold'] for x in gene_locations])
			for scaff in scaffolds:
				gc_gene_locations = {}
				for g in gene_locations:
					if gene_locations[g]['scaffold'] == scaff:
						gc_gene_locations[g] = gene_locations[g]
						scaff_relative_gene_locations[g] = gene_locations[g]
				gc_gene_counts[gbk + "|" + scaff] = len(gc_gene_locations)
				gc_genes[gbk + "|" + scaff] = set(gc_gene_locations.keys())

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
				ginfo = scaff_relative_gene_locations[g]
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
			if hps[0][0] in core_hgs and hps[0][1] in core_hgs:
				anchor_edge = hps[0]
				break
		try:
			assert (anchor_edge != None)
		except:
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

def determineHGStats(orthogroup_matrix_file, hg_nucl_dir, logObject, representative_associated_members=None, impute_broad_conservation=False):
	try:
		hg_single_copy_status = {}
		hg_prop_samples = {}
		hg_lts = defaultdict(set)
		samples = []
		with open(orthogroup_matrix_file) as omf:
			for i, line in enumerate(omf):
				line = line.rstrip('\n')
				ls = line.split('\t')
				if i == 0:
					samples = ls[1:]
					continue
				hg = ls[0]
				is_single_copy = True
				sample_count = 0
				weighted_count = 0
				total_weighted_count = 0
				for j, lts in enumerate(ls[1:]):
					samp = samples[j]
					if representative_associated_members != None and impute_broad_conservation:
						total_weighted_count += len(representative_associated_members[samp])
					if ',' in lts:
						is_single_copy = False
					if lts.strip() != '':
						sample_count += 1
						if representative_associated_members != None and impute_broad_conservation:
							weighted_count += len(representative_associated_members[samp])
					for lt in lts.split(', '):
						if lt.strip() == '': continue
						hg_lts[hg].add(lt)
				hg_single_copy_status[hg] = is_single_copy
				if representative_associated_members != None and impute_broad_conservation:
					hg_prop_samples[hg] = weighted_count/float(total_weighted_count)
				else:
					hg_prop_samples[hg] = sample_count/float(len(ls[1:]))
		hg_median_lengths = {}
		hg_median_gcskew = {}
		hg_median_gc = {}
		for f in os.listdir(hg_nucl_dir):
			hg = f.split('.fna')[0]
			lengths = []
			gcs = []
			gc_skews = []
			with open(hg_nucl_dir + f) as ohpf:
				for rec in SeqIO.parse(ohpf, 'fasta'):
					seq = str(rec.seq)
					tot_bases = len(seq)
					g = sum([1 for x in seq if x == 'C'])
					c = sum([1 for x in seq if x == 'G'])
					gc_sum = g+c
					gc = gc_sum/tot_bases
					gcs.append(gc)
					gc_skew = (g-c)/gc_sum
					gc_skews.append(gc_skew)
					lengths.append(len(str(rec.seq)))
			hg_median_lengths[hg] = statistics.median(lengths)
			hg_median_gcskew[hg] = statistics.median(gc_skews)
			hg_median_gc[hg] = statistics.median(gcs)
		return([hg_single_copy_status, hg_prop_samples, hg_median_lengths, hg_median_gcskew, hg_median_gc, dict(hg_lts)])
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

		unique_seqs = set([])
		align_len = 0
		with open(hg_codo_algn_file) as ohcaf:
			for rec in SeqIO.parse(ohcaf, 'fasta'):
				unique_seqs.add(str(rec.seq))
				align_len = len(str(rec.seq))

		if len(unique_seqs) == 1:
			return

		if align_len <= 200:
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
				logObject.error('Had an issue running GARD: %s' % ' '.join(gard_cmd))
				sys.stderr.write('Had an issue running GARD: %s\n' % ' '.join(gard_cmd))
				logObject.error(e)
				sys.exit(1)

			fubar_cmd = ['hyphy', 'CPU=1', 'fubar', '--alignment', best_gard_output]
			#print(' '.join(fubar_cmd))
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
		fubar_deba = {}
		for f in os.listdir(fubar_results_dir):
			if f.endswith('.json'):
				hg = f.split('.msa.fna.FUBAR.json')[0]
				if f.endswith('.best.FUBAR.json'):
					hg = f.split('.best.FUBAR.json')[0]
				fubar_json_result = fubar_results_dir + f
				try:
					with open(fubar_json_result) as ofjr:
						fubar_results = json.load(ofjr)
					pos_selected_sites = 0
					neg_selected_sites = 0
					sum_deba = 0
					tot_sites= 0
					for partition in fubar_results['MLE']['content']:
						for site_mle_info in fubar_results['MLE']['content'][partition]:
							tot_sites += 1
							alpha, beta, diff, prob_agb, prob_alb, bayesfactor, _, _ = site_mle_info
							sum_deba += (beta - alpha)
							if prob_agb >= 0.9:
								neg_selected_sites += 1
							if prob_alb >= 0.9:
								pos_selected_sites += 1
					tot_selected_sites = pos_selected_sites + neg_selected_sites
					prop_selected_sites_positive = 'NA'
					if tot_selected_sites >= 1:
						prop_selected_sites_positive = float(pos_selected_sites)/float(neg_selected_sites+pos_selected_sites)
					fubar_sel_props[hg] = prop_selected_sites_positive
					fubar_sel_sites[hg] = tot_selected_sites # /float(tot_sites) TODO: make this proportion - more useful!!
					avg_deba = "NA"
					if tot_sites > 0:
						avg_deba = sum_deba/float(tot_sites)
					fubar_deba[hg] = avg_deba
					# TODO: process "grid" field in FUBAR results to get most probable dN/dS ratio
				except:
					fubar_sel_props[hg] = 'NA'
					fubar_sel_sites[hg] = 'NA'
					fubar_deba[hg] = 'NA'

		return([gard_partitions, fubar_sel_props, fubar_sel_sites, fubar_deba])
	except Exception as e:
		sys.stderr.write('Issues with running GARD or FUBAR analyses.\n')
		logObject.error('Issues with running GARD or FUBAR analyses.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def determineSeqSimProteinAlignment(inputs):
	use_only_core = True # hardcoded true at the moment
	hg, protein_alignment_file, outf, logObject = inputs
	protein_sequences = {}
	with open(protein_alignment_file) as ocaf:
		for rec in SeqIO.parse(ocaf, 'fasta'):
			protein_sequences[rec.id] = str(rec.seq).upper()

	pair_seq_matching = defaultdict(lambda: defaultdict(lambda: 0.0))
	for i, g1 in enumerate(sorted(protein_sequences)):
		s1 = g1.split('|')[0]
		g1s = protein_sequences[g1]
		for j, g2 in enumerate(sorted(protein_sequences)):
			if i >= j: continue
			s2 = g2.split('|')[0]
			if s1 == s2: continue
			g2s = protein_sequences[g2]
			tot_comp_pos = 0
			match_pos = 0
			for pos, g1a in enumerate(g1s):
				g2a = g2s[pos]
				if g1a != '-' or g2a != '-':
					if not use_only_core or (use_only_core and g1a != '-' and g2a != '-'):
						tot_comp_pos += 1
						if g1a == g2a:
							match_pos += 1
			general_matching_percentage = 0.0
			if tot_comp_pos > 0:
				general_matching_percentage = float(match_pos) / float(tot_comp_pos)
			if pair_seq_matching[s1][s2] < general_matching_percentage and pair_seq_matching[s2][s1] < general_matching_percentage:
				pair_seq_matching[s1][s2] = general_matching_percentage
				pair_seq_matching[s2][s1] = general_matching_percentage

	pair_seq_matching_normal = default_to_regular(pair_seq_matching)
	with open(outf, 'wb') as pickle_file:
		pickle.dump(pair_seq_matching_normal, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)

def computeBetaRDgc(prot_algn_dir, evo_dir, logObject, cpus=1):
	"""
	Note, Beta-RD gene-cluster statistic here is being computed in a different manner than what we did in lsaBGC.
	"""
	brd_results_dir = evo_dir + 'BetaRDgc_Calculations/'
	util.setupReadyDirectory([brd_results_dir])
	hg_med_beta_rd = {}
	hg_max_beta_rd = {}
	try:
		inputs = []
		for f in os.listdir(prot_algn_dir):
			hg = f.split('.msa.faa')[0]
			outf = brd_results_dir + hg + '.sims.pkl'
			inputs.append([hg, prot_algn_dir + f, outf, logObject])

		p = multiprocessing.Pool(cpus)
		p.map(determineSeqSimProteinAlignment, inputs)
		p.close()

		hg_sims_dict = {}
		gc_wide_sims_dict = defaultdict(lambda: defaultdict(list))
		for f in os.listdir(brd_results_dir):
			hg = f.split('.sims.pkl')[0]
			sims_dict = None
			with open(brd_results_dir + f, 'rb') as handle:
				sims_dict = pickle.load(handle)
			if len(sims_dict) < 2: continue
			hg_sims_dict[hg] = sims_dict
			for i, s1 in enumerate(sorted(sims_dict)):
				for j, s2 in enumerate(sorted(sims_dict)):
					if s1 != s2 and s2 in sims_dict[s1]:
						gc_wide_sims_dict[s1][s2].append(sims_dict[s1][s2])
					else:
						gc_wide_sims_dict[s1][s2].append(0.0)

		for hg in hg_sims_dict:
			Brdgc = []
			for i, s1 in enumerate(sorted(hg_sims_dict[hg])):
				for j, s2 in enumerate(sorted(hg_sims_dict[hg])):
					if i >= j: continue
					Brdgc.append(hg_sims_dict[hg][s1][s2] / float(statistics.median(gc_wide_sims_dict[s1][s2])))
			if len(Brdgc) > 0:
				hg_med_beta_rd[hg] = statistics.median(Brdgc)
				hg_max_beta_rd[hg] = max(Brdgc)
	except Exception as e:
		sys.stderr.write('Issues with calculating Beta-RD gene-cluster for homolog groups.\n')
		logObject.error('Issues with calculating Beta-RD gene-cluster for homolog groups.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return([hg_med_beta_rd, hg_max_beta_rd])

def runEntropyAnalysis(codo_algn_trim_dir, upst_algn_dir, evo_dir, logObject, cpus=1):
	try:
		entropy_res_dir = evo_dir + 'Entropy_Calculations/'
		util.setupReadyDirectory([entropy_res_dir])
		inputs = []
		for f in os.listdir(codo_algn_trim_dir):
			hg = f.split('.msa.fna')[0]
			caf = codo_algn_trim_dir + f
			outf = entropy_res_dir + hg + '_codon.txt'
			inputs.append([hg, caf, outf, logObject])

		for f in os.listdir(upst_algn_dir):
			hg = f.split('.msa.fna')[0]
			uaf = upst_algn_dir + f
			outf = entropy_res_dir + hg + '_upstream.txt'
			inputs.append([hg, uaf, outf, logObject])

		p = multiprocessing.Pool(cpus)
		p.map(calculateMSAEntropy, inputs)
		p.close()

		hg_entropy = {}
		hg_upst_entropy = {}
		for f in os.listdir(entropy_res_dir):
			with open(entropy_res_dir + f) as oef:
				for line in oef:
					line = line.strip()
					hg, ep = line.split('\t')
					if f.endswith('_upstream.txt'):
						hg_upst_entropy[hg] = ep
					elif f.endswith('_codon.txt'):
						hg_entropy[hg] = ep
		return([hg_entropy, hg_upst_entropy])
	except Exception as e:
		sys.stderr.write('Issues with calculating entropy for homolog groups or their upstream regions.\n')
		logObject.error('Issues with calculating entropy for homolog groups or their upstream regions.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def calculateMSAEntropy(inputs):
	hg, nucl_algn_fasta, outf, logObject = inputs
	try:
		seqs = []
		with open(nucl_algn_fasta) as onaf:
			for rec in SeqIO.parse(onaf, 'fasta'):
				seqs.append(list(str(rec.seq)))
		accounted_sites = 0
		all_entropy = 0.0
		for tup in zip(*seqs):
			als = list(tup)
			missing_prop = sum([1 for al in als if not al in set(['A', 'C', 'G', 'T'])])/float(len(als))
			if missing_prop >= 0.1: continue
			filt_als = [al for al in als if al in set(['A', 'C', 'G', 'T'])]
			a_freq = sum([1 for al in filt_als if al == 'A'])/float(len(filt_als))
			c_freq = sum([1 for al in filt_als if al == 'C'])/float(len(filt_als))
			g_freq = sum([1 for al in filt_als if al == 'G'])/float(len(filt_als))
			t_freq = sum([1 for al in filt_als if al == 'T'])/float(len(filt_als))
			site_entropy = stats.entropy([a_freq, c_freq, g_freq, t_freq],base=4)
			all_entropy += site_entropy
			accounted_sites += 1
		avg_entropy = "NA"
		if accounted_sites > 0:
			avg_entropy = all_entropy/accounted_sites
		outf_handle = open(outf, 'w')
		outf_handle.write(hg + '\t' + str(avg_entropy) + '\n')
		outf_handle.close()
	except Exception as e:
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def calculateAmbiguity(codo_algn_dir, codo_algn_trim_dir, logObject):
	full_amb_prop = {}
	trim_amb_prop = {}
	try:
		for caf in os.listdir(codo_algn_dir):
			if not caf.endswith('.msa.fna'): continue
			hg = caf.split('.msa.fna')[0]
			codo_algn_file = codo_algn_dir + caf
			codo_sequences = []
			with open(codo_algn_file) as ocaf:
				for rec in SeqIO.parse(ocaf, 'fasta'):
					codo_sequences.append(list(str(rec.seq)))
			tot = 0
			amb = 0
			for al in zip(*codo_sequences):
				tot += 1
				all = list(al)
				amb_site_prop = sum([1 for x in all if not x in set(['A', 'C', 'G', 'T'])])/float(len(all))
				if amb_site_prop >= 0.1:
					amb += 1
			amb_prop = float(amb)/float(tot)
			full_amb_prop[hg] = amb_prop

		for catf in os.listdir(codo_algn_trim_dir):
			if not catf.endswith('.msa.fna'): continue
			hg = catf.split('.msa.fna')[0]
			codo_algn_trimmed_file = codo_algn_trim_dir + catf
			codo_sequences = []
			with open(codo_algn_trimmed_file) as ocatf:
				for rec in SeqIO.parse(ocatf, 'fasta'):
					codo_sequences.append(list(str(rec.seq)))
			tot = 0
			amb = 0
			for al in zip(*codo_sequences):
				tot += 1
				all = list(al)
				amb_site_prop = sum([1 for x in all if not x in set(['A', 'C', 'G', 'T'])])/float(len(all))
				if amb_site_prop >= 0.1:
					amb += 1
			amb_prop = float(amb)/float(tot)
			trim_amb_prop[hg] = amb_prop

	except Exception as e:
		sys.stderr.write('Issues with calculating ambiguity for full or trimmed codon alignments.\n')
		logObject.error('Issues with calculating ambiguity for full or trimmed codon alignments.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)
	return([full_amb_prop, trim_amb_prop])

def runTajimasDAnalysisPerHG(inputs):
	hg, trim_codon_align, outf, logObject = inputs
	try:
		outf_handle = open(outf, 'w')
		codo_sequences = []
		with open(trim_codon_align) as ocatf:
			for rec in SeqIO.parse(ocatf, 'fasta'):
				codo_sequences.append(str(rec.seq))
		# at least 4 sequences and 60 bp in filtered alignment
		if len(codo_sequences) >= 4 and len(codo_sequences[0]) > 60:
			taj_d, seg_sites = calculateTajimasD(codo_sequences)
			seg_sites_prop = seg_sites/len(codo_sequences[0])
		else:
			taj_d = 'NA'
			seg_sites_prop = 'NA'
		outf_handle.write('\t'.join([hg, str(taj_d), str(seg_sites_prop)]) + '\n')
	except Exception as e:
		sys.stderr.write('Issues with calculating Tajima\'s D for homolog group %s.\n' % hg)
		logObject.error('Issues with calculating Tajima\'s D for homolog group %s.' % hg)
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def runTajimasDAnalysis(codo_algn_trim_dir, evo_dir, logObject, cpus=1):
	try:
		tajd_resdir = evo_dir + 'TajimasD_and_SegSites_Calculations/'
		util.setupReadyDirectory([tajd_resdir])

		inputs = []
		for catf in os.listdir(codo_algn_trim_dir):
			if not catf.endswith('.msa.fna'): continue
			hg = catf.split('.msa.fna')[0]
			trim_codon_align = codo_algn_trim_dir + catf
			outf = tajd_resdir + hg + '.txt'
			inputs.append([hg, trim_codon_align, outf, logObject])

		p = multiprocessing.Pool(cpus)
		p.map(runTajimasDAnalysisPerHG, inputs)
		p.close()

		hg_tajimas_d = {}
		hg_seg_sites_prop = {}
		for f in os.listdir(tajd_resdir):
			with open(tajd_resdir + f) as otf:
				for line in otf:
					line = line.strip()
					hg, tajd, ssp = line.split('\t')
					hg_tajimas_d[hg] = tajd
					hg_seg_sites_prop[hg] = ssp
		return([hg_tajimas_d, hg_seg_sites_prop])
	except Exception as e:
		sys.stderr.write('Issues with calculating Tajima\'s D for homolog groups.\n')
		logObject.error('Issues with calculating Tajima\'s D for homolog groups.')
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
		return ([D, S])
	else:
		return (["< 3 segregating sites!", S])


def compareFocalAndComparatorGeneClusters(focal_genbank_ids, comparator_genbank_ids, codo_algn_trim_dir, upst_algn_dir, logObject, representative_associated_members=None, impute_broad_conservation=False):
	comp_stats = {}
	try:
		total_foc_broad = set([])
		total_com_broad = set([])
		if impute_broad_conservation and representative_associated_members != None:
			for gc in focal_genbank_ids:
				total_foc_broad.add(gc)
				for orthogc in representative_associated_members[gc]:
					total_foc_broad.add(orthogc)
			for gc in comparator_genbank_ids:
				total_com_broad.add(gc)
				for orthogc in representative_associated_members[gc]:
					total_com_broad.add(orthogc)
		for f in os.listdir(codo_algn_trim_dir):
			hg = f.split('.msa.fna')[0]
			codo_algn_trim_file = codo_algn_trim_dir + f
			focal_samps_with_hg = set([])
			focal_samps_with_hg_broad = set([])
			compa_samps_with_hg = set([])
			compa_samps_with_hg_broad = set([])
			focal_seqs = []
			compa_seqs = []
			with open(codo_algn_trim_file) as opatf:
				for rec in SeqIO.parse(opatf, 'fasta'):
					sample = rec.id.split('|')[0]
					seq = str(rec.seq)
					if sample in focal_genbank_ids:
						focal_samps_with_hg.add(sample)
						focal_samps_with_hg_broad.add(sample)
						if impute_broad_conservation and representative_associated_members != None:
							for orthogc in representative_associated_members[sample]:
								focal_samps_with_hg_broad.add(orthogc)
						focal_seqs.append(seq)
					elif sample in comparator_genbank_ids:
						compa_samps_with_hg.add(sample)
						compa_samps_with_hg_broad.add(sample)
						if impute_broad_conservation and representative_associated_members != None:
							for orthogc in representative_associated_members[sample]:
								compa_samps_with_hg_broad.add(orthogc)
						compa_seqs.append(seq)

			diff_between = 0
			pw_between = 0
			diff_foc_within = 0
			pw_foc_within = 0
			for i, s1 in enumerate(focal_seqs):
				for j, s2 in enumerate(focal_seqs):
					if i >= j: continue
					diff_foc_within += sum(1 for a, b in zip(s1, s2) if a != b and a != '-' and b != '-')
					pw_foc_within += 1

			for i, s1 in enumerate(focal_seqs):
				for j, s2 in enumerate(compa_seqs):
					diff_between += sum(1 for a, b in zip(s1, s2) if a != b and a != '-' and b != '-')
					pw_between += 1

			# Fst estimated according to Hudson, Slatkin and Maddison 1989
			# which is closely related to Nei and Chesser 1983.
			# While the derivations were specific to diploid organisms,
			# the concept of the estimation can more simply be applied
			# to haploid and that is what is assumed here.
			pi_between, pi_within, fst = ['NA']*3
			if pw_between > 0 and pw_foc_within > 0:
				pi_between = diff_between/float(pw_between)
				pi_within = (diff_foc_within)/float(pw_foc_within)
				if pi_between > 0:
					fst = 1.0 - (float(pi_within)/float(pi_between))

			#pi_foc = diff_foc_within/float(pw_foc_within)
			#pi_com = diff_com_within/float(pw_com_within)

			if impute_broad_conservation and representative_associated_members != None:
				prop_foc_with = len(focal_samps_with_hg_broad)/float(len(total_foc_broad))
				prop_com_with = len(compa_samps_with_hg_broad)/float(len(total_com_broad))
			else:
				prop_foc_with = len(focal_samps_with_hg)/float(len(focal_genbank_ids))
				prop_com_with = len(compa_samps_with_hg)/float(len(comparator_genbank_ids))

			upst_algn_file = upst_algn_dir + hg + '.msa.fna'
			upst_fst = 'NA'
			if os.path.isfile(upst_algn_file) and os.path.getsize(upst_algn_file) > 0:
				focal_samps_with_hg = set([])
				compa_samps_with_hg = set([])
				focal_seqs = []
				compa_seqs = []
				with open(upst_algn_file) as oaf:
					for rec in SeqIO.parse(oaf, 'fasta'):
						sample = rec.id.split('|')[0]
						seq = str(rec.seq)
						if sample in focal_genbank_ids:
							focal_samps_with_hg.add(sample)
							focal_seqs.append(seq)
						elif sample in comparator_genbank_ids:
							compa_samps_with_hg.add(sample)
							compa_seqs.append(seq)

				diff_between = 0
				pw_between = 0
				diff_foc_within = 0
				diff_com_within = 0
				pw_foc_within = 0
				pw_com_within = 0
				for i, s1 in enumerate(focal_seqs):
					for j, s2 in enumerate(focal_seqs):
						if i >= j: continue
						diff_foc_within += sum(1 for a, b in zip(s1, s2) if a != b and a != '-' and b != '-')
						pw_foc_within += 1

				for i, s1 in enumerate(compa_seqs):
					for j, s2 in enumerate(compa_seqs):
						if i >= j: continue
						diff_com_within += sum(1 for a, b in zip(s1, s2) if a != b and a != '-' and b != '-')
						pw_com_within += 1

				for i, s1 in enumerate(focal_seqs):
					for j, s2 in enumerate(compa_seqs):
						diff_between += sum(1 for a, b in zip(s1, s2) if a != b and a != '-' and b != '-')
						pw_between += 1

				#print([diff_between, pw_between, pw_foc_within, pw_com_within, diff_foc_within, diff_com_within])
				if pw_between > 0 and pw_foc_within > 0:
					pi_between = diff_between/float(pw_between)
					pi_within = (diff_foc_within)/float(pw_foc_within)
					if pi_between > 0:
						upst_fst = 1.0 - (float(pi_within)/float(pi_between))

			comp_stats[hg] = {'prop_foc_with': prop_foc_with, 'prop_com_with': prop_com_with, 'fst': fst, 'fst_upst': upst_fst}

	except Exception as e:
		sys.stderr.write('Issues with performing comparative analyses between user-defined Gene Cluster groups.\n')
		logObject.error('Issues with performing comparative analyses between user-defined Gene Cluster groups.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return (comp_stats)

def consolidateReport(consensus_prot_seqs_faa, comp_stats, hg_stats, annotations, evo_stats, final_report_xlsx, final_report_tsv, logObject, run_hyphy=False, ces=False):
	"""
	dict_keys(['pfam', 'vfdb', 'paperblast', 'pgap', 'vog', 'isfinder', 'card', 'mibig'])
	dict_keys(['hg_single_copy_status', 'hg_prop_samples', 'hg_median_lengths', 'hg_order_scores'])
	dict_keys(['tajimas_d', 'gard_partitions', 'fubar_sel_props', 'fubar_sel_sites', 'gene_tree_congruence'])
	"""
	try:
		# Note to self, eventually quit being lazy and conditionally display all columns (e.g. FUBAR columns) when requested by user
		# to avoid having columns with NA values.
		header = ['Homolog Group (HG) ID', 'HG is Single Copy?', 'Proportion of Total Gene Clusters with HG',
				  'HG Median Length (bp)', 'HG Consensus Order', 'HG Consensus Direction']
		if comp_stats != None:
			header += ['Proportion of Focal Gene Clusters with HG', 'Proportion of Comparator Gene Clusters with HG',
					   'Fixation Index', 'Upstream Region Fixation Index']
		header += ['Tajima\'s D', 'Proportion of Filtered Codon Alignment is Segregating Sites', 'Entropy',
				   'Upstream Region Entropy', 'Median Beta-RD-gc', 'Max Beta-RD-gc',
				   'Proportion of sites which are highly ambiguous in codon alignment',
				   'Proportion of sites which are highly ambiguous in trimmed codon alignment', 'Median GC',
				   'Median GC Skew']
		if run_hyphy:
			header += ['GARD Partitions Based on Recombination Breakpoints',
			           'Number of Sites Identified as Under Positive or Negative Selection by FUBAR',
				       'Average delta(Beta, Alpha) by FUBAR across sites',
				       'Proportion of Sites Under Selection which are Positive']
		header += ['Custom Annotation (E-value)', 'KO Annotation (E-value)', 'PGAP Annotation (E-value)',
				   'PaperBLAST Annotation (E-value)', 'CARD Annotation (E-value)', 'IS Finder (E-value)',
				   'MI-BiG Annotation (E-value)', 'VOG Annotation (E-value)',  'VFDB Annotation (E-value)',
				   'Pfam Domains', 'CDS Locus Tags', 'HG Consensus Sequence']

		seqs = {}
		with open(consensus_prot_seqs_faa) as ocpsf:
			for rec in SeqIO.parse(ocpsf, 'fasta'):
				seqs[rec.id] = str(rec.seq)

		frt_handle = open(final_report_tsv, 'w')
		frt_handle.write('\t'.join(header) + '\n')
		# traverse HG in consensus order
		num_rows = 1
		for hg_tup in sorted(hg_stats['hg_order_scores'].items(), key=lambda e: e[1][0]):
			hg = hg_tup[0]
			if not hg in hg_stats['hg_median_lengths']: continue
			hg_scs = hg_stats['hg_single_copy_status'][hg]
			hg_cons = hg_stats['hg_prop_samples'][hg]
			hg_mlen = hg_stats['hg_median_lengths'][hg]
			hg_lts = '; '.join(sorted(hg_stats['hg_locus_tags'][hg]))
			hg_full_amb = hg_stats['hg_full_ambiguity'][hg]
			hg_trim_amb = hg_stats['hg_trim_ambiguity'][hg]
			hg_gc = hg_stats['hg_median_gc'][hg]
			hg_gcs = hg_stats['hg_median_gcskew'][hg]
			hg_ordr = hg_tup[1][0]
			hg_dire = '"' + hg_tup[1][1] + '"'
			hg_tajd, hg_entr, hg_upst_entr, hg_segs, hg_gpar, hg_ssit, hg_deba, hg_spro, hg_med_brdgc, hg_max_brdgc, fst, fst_upst = ['NA']*12
			if hg_scs == True or ces:
				hg_tajd = util.gatherValueFromDictForHomologGroup(hg, evo_stats['tajimas_d'])
				hg_entr = util.gatherValueFromDictForHomologGroup(hg, evo_stats['entropy'])
				hg_upst_entr = util.gatherValueFromDictForHomologGroup(hg, evo_stats['entropy_upst'])
				hg_segs = util.gatherValueFromDictForHomologGroup(hg, evo_stats['segregating_sites_prop'])
				hg_gpar = util.gatherValueFromDictForHomologGroup(hg, evo_stats['gard_partitions'])
				hg_ssit = util.gatherValueFromDictForHomologGroup(hg, evo_stats['fubar_sel_sites'])
				hg_spro = util.gatherValueFromDictForHomologGroup(hg, evo_stats['fubar_sel_props'])
				hg_deba = util.gatherValueFromDictForHomologGroup(hg, evo_stats['fubar_dba'])
				hg_med_brdgc = util.gatherValueFromDictForHomologGroup(hg, evo_stats['median_beta_rd_gc'])
				hg_max_brdgc = util.gatherValueFromDictForHomologGroup(hg, evo_stats['max_beta_rd_gc'])
			cust_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, 'custom', annotations)
			ko_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, 'ko', annotations)
			pgap_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, 'pgap', annotations)
			pb_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, 'paperblast', annotations)
			card_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, 'card', annotations)
			isf_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, 'isfinder', annotations)
			mibig_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, 'mibig', annotations)
			vog_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, 'vog', annotations)
			vfdb_annot = util.gatherAnnotationFromDictForHomoloGroup(hg, 'vfdb', annotations)
			pfam_annots = 'NA'
			if 'pfam' in annotations and hg in annotations['pfam']:
				pfam_annots = '; '.join(annotations['pfam'][hg][0])
			con_seq = seqs[hg]
			row = [hg, hg_scs, hg_cons, hg_mlen, hg_ordr, hg_dire]
			if comp_stats != None:
				fp = comp_stats[hg]['prop_foc_with']
				cp = comp_stats[hg]['prop_com_with']
				if hg_scs == True or ces:
					fst = comp_stats[hg]['fst']
					fst_upst = comp_stats[hg]['fst_upst']
				row += [fp, cp, fst, fst_upst]
			row += [hg_tajd, hg_segs, hg_entr, hg_upst_entr, hg_med_brdgc, hg_max_brdgc, hg_full_amb, hg_trim_amb,
					hg_gc, hg_gcs]
			if run_hyphy:
				row += [hg_gpar, hg_ssit, hg_deba, hg_spro]
			row += [cust_annot, ko_annot, pgap_annot, pb_annot, card_annot, isf_annot, mibig_annot, vog_annot,
					vfdb_annot, pfam_annots, hg_lts, con_seq]
			row = [str(x) for x in row]
			frt_handle.write('\t'.join(row) + '\n')
			num_rows += 1
		frt_handle.close()

		# Generate Excel spreadsheet
		writer = pd.ExcelWriter(final_report_xlsx, engine='xlsxwriter')
		workbook = writer.book
		dd_sheet = workbook.add_worksheet('Data Dictionary')
		dd_sheet.write(0, 0, 'WARNING: some evolutionary statistics are experimental - evaluate with caution!!!')
		dd_sheet.write(1, 0, 'Data Dictionary describing columns of "Overview" spreadsheets can be found on zol\'s Wiki at:')
		dd_sheet.write(2, 0, 'https://github.com/Kalan-Lab/zol/wiki/3.-more-info-on-zol#explanation-of-report')

		numeric_columns = {'Proportion of Total Gene Clusters with HG', 'Proportion of Focal Gene Clusters with HG',
						   'Proportion of Comparator Gene Clusters with HG', 'Fixation Index',
						   'Upstream Region Fixation Index', 'HG Median Length (bp)', 'HG Consensus Order',
						   'Tajima\'s D', 'Entropy', 'Upstream Region Entropy',
						   'GARD Partitions Based on Recombination Breakpoints',
						   'Number of Sites Identified as Under Positive or Negative Selection by FUBAR',
						   'Proportion of Sites Under Selection which are Positive', 'Median Beta-RD-gc',
						   'Max Beta-RD-gc', 'Proportion of Filtered Codon Alignment is Segregating Sites',
						   'Proportion of sites which are highly ambiguous in codon alignment',
						   'Proportion of sites which are highly ambiguous in trimmed codon alignment',
						   'Average delta(Beta, Alpha) by FUBAR across sites', 'Median GC', 'Median GC Skew'}

		warn_format = workbook.add_format({'bg_color': '#bf241f', 'bold': True, 'font_color': '#FFFFFF'})
		na_format = workbook.add_format({'font_color': '#a6a6a6', 'bg_color': '#FFFFFF', 'italic': True})
		header_format = workbook.add_format({'bold': True, 'text_wrap': True, 'valign': 'top', 'fg_color': '#D7E4BC', 'border': 1})

		results_df = util.loadTableInPandaDataFrame(final_report_tsv, numeric_columns)
		results_df.to_excel(writer, sheet_name='ZoL Results', index=False, na_rep="NA")
		worksheet =  writer.sheets['ZoL Results']
		worksheet.conditional_format('B2:B' + str(num_rows), {'type': 'cell', 'criteria': '==', 'value': '"False"', 'format': warn_format})
		worksheet.conditional_format('A2:BA' + str(num_rows), {'type': 'cell', 'criteria': '==', 'value': '"NA"', 'format': na_format})
		worksheet.conditional_format('A1:BA1', {'type': 'cell', 'criteria': '!=', 'value': 'NA', 'format': header_format})

		# prop gene-clusters with hg
		worksheet.conditional_format('C2:C' + str(num_rows), {'type': '2_color_scale', 'min_color': "#f7de99", 'max_color': "#c29006", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})
		# gene-lengths
		worksheet.conditional_format('D2:D' + str(num_rows), {'type': '2_color_scale', 'min_color': "#a3dee3", 'max_color': "#1ebcc9", "min_value": 100, "max_value": 2500, 'min_type': 'num', 'max_type': 'num'})
		if comp_stats != None:
			# prop focal gene-clusters with hg
			worksheet.conditional_format('G2:G' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#f7de99", 'max_color': "#c29006",'min_type': 'num', 'max_type': 'num',
										  "min_value": 0.0, "max_value": 1.0})
			# prop comparator gene-clusters with hg
			worksheet.conditional_format('H2:H' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#f7de99", 'max_color': "#c29006",'min_type': 'num', 'max_type': 'num',
										  "min_value": 0.0, "max_value": 1.0})
			# fst
			worksheet.conditional_format('I2:I' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#aecaf5", 'max_color': "#6198ed",'min_type': 'num', 'max_type': 'num',
										  "min_value": 0.0, "max_value": 1.0})

			# upstream region fst
			worksheet.conditional_format('J2:J' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#aecaf5", 'max_color': "#6198ed",'min_type': 'num', 'max_type': 'num',
										  "min_value": 0.0, "max_value": 1.0})

			# taj-d
			worksheet.conditional_format('K2:K' + str(num_rows),
										 {'type': '3_color_scale', 'min_color': "#f7a09c", "mid_color": "#e0e0e0", 'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										  'max_color': "#87cefa", "min_value": -2.0, "mid_value": 0.0, "max_value": 2.0})
			# prop seg sites
			worksheet.conditional_format('L2:L' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#eab3f2", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#a37ba8", "min_value": 0.0, "max_value": 1.0})

			# entropy
			worksheet.conditional_format('M2:M' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#f7a8bc", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#fa6188", "min_value": 0.0, "max_value": 1.0})

			# upstream region entropy
			worksheet.conditional_format('N2:N' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#f7a8bc", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#fa6188", "min_value": 0.0, "max_value": 1.0})

			# beta-rd gc
			worksheet.conditional_format('O2:O' + str(num_rows),
										 {'type': '3_color_scale', 'min_color': "#fac087", "mid_color": "#e0e0e0", 'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										 'max_color': "#9eb888", "min_value": 0.75, "mid_value": 1.0, "max_value": 1.25})
			# max beta-rd gc
			worksheet.conditional_format('P2:P' + str(num_rows),
										 {'type': '3_color_scale', 'min_color': "#fac087", "mid_color": "#e0e0e0", 'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										  'max_color': "#9eb888", "min_value": 0.75, "mid_value": 1.0, "max_value": 1.25})

			# ambiguity full ca
			worksheet.conditional_format('Q2:Q' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#ed8c8c", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#ab1616", "min_value": 0.0, "max_value": 1.0})

			# ambiguity trim ca
			worksheet.conditional_format('R2:R' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#ed8c8c", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#ab1616", "min_value": 0.0, "max_value": 1.0})

			# GC
			worksheet.conditional_format('S2:S' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#abffb7", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#43bf55", "min_value": 0.0, "max_value": 1.0})

			# GC Skew
			worksheet.conditional_format('T2:T' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#c7afb4", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#965663", "min_value": -2.0, "max_value": 2.0})

		else:
			# taj-d
			worksheet.conditional_format('G2:G' + str(num_rows),
										 {'type': '3_color_scale', 'min_color': "#f7a09c", "mid_color": "#e0e0e0",'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										  'max_color': "#87cefa", "min_value": -2.0, "mid_value": 0.0, "max_value": 2.0})

			# prop seg sites
			worksheet.conditional_format('H2:H' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#eab3f2", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#a37ba8", "min_value": 0.0, "max_value": 1.0})

			# entropy
			worksheet.conditional_format('I2:I' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#f7a8bc", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#fa6188", "min_value": 0.0, "max_value": 1.0})

			# upstream region entropy
			worksheet.conditional_format('J2:J' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#f7a8bc", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#fa6188", "min_value": 0.0, "max_value": 1.0})

			# median beta-rd gc
			worksheet.conditional_format('K2:K' + str(num_rows),
										 {'type': '3_color_scale', 'min_color': "#fac087", "mid_color": "#e0e0e0",'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										  'max_color': "#9eb888", "min_value": 0.75, "mid_value": 1.0, "max_value": 1.25})
			# max beta-rd gc
			worksheet.conditional_format('L2:L' + str(num_rows),
										 {'type': '3_color_scale', 'min_color': "#fac087", "mid_color": "#e0e0e0",'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										  'max_color': "#9eb888", "min_value": 0.75, "mid_value": 1.0, "max_value": 1.25})

			# ambiguity full ca
			worksheet.conditional_format('M2:M' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#ed8c8c", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#ab1616", "min_value": 0.0, "max_value": 1.0})

			# ambiguity trim ca
			worksheet.conditional_format('N2:N' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#ed8c8c", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#ab1616", "min_value": 0.0, "max_value": 1.0})

			# GC
			worksheet.conditional_format('O2:O' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#abffb7", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#43bf55", "min_value": 0.0, "max_value": 1.0})

			# GC Skew
			worksheet.conditional_format('P2:P' + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#c7afb4", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#965663", "min_value": -2.0, "max_value": 2.0})

		worksheet.autofilter('A1:BA' + str(num_rows))
		worksheet.filter_column(2, 'x >= 0.1')
		workbook.close()

	except Exception as e:
		sys.stderr.write('Issues creating consolidated results files.\n')
		logObject.error('Issues creating consolidated results files.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def plotHeatmap(hg_stats, genbanks, plot_result_pdf, work_dir, logObject, height=7, width=10, full_genbank_labels=False):
	try:
		representative_genbanks = set([])
		for gbk in genbanks:
			gbk_prefix = None
			if gbk.endswith('.gbk') or gbk.endswith('.gbff') or gbk.endswith('.genbank'):
				gbk_prefix = '.'.join(gbk.split('/')[-1].split('.')[:-1])
			assert (gbk_prefix != None)
			representative_genbanks.add(gbk_prefix)

		# create input tracks
		ml_track_file = work_dir + 'HG_Median_Length_Info.txt'
		hm_track_file = work_dir + 'HG_Heatmap_Info.txt'
		ml_track_handle = open(ml_track_file, 'w')
		hm_track_handle = open(hm_track_file, 'w')
		ml_track_handle.write('og\tog_order\tmed_length\n')
		hm_track_handle.write('og\tog_order\tgenbank\tog_presence\tcopy_count\n')
		gn_lab_keys = set([])
		gn_labs = set([])
		for hg_tup in sorted(hg_stats['hg_order_scores'].items(), key=lambda e: e[1][0]):
			hg = hg_tup[0]
			if not hg in hg_stats['hg_median_lengths']: continue
			hg_mlen = hg_stats['hg_median_lengths'][hg]
			hg_lts = hg_stats['hg_locus_tags'][hg]
			hg_ordr = hg_tup[1][0]
			sample_copy_counts = defaultdict(int)
			for lt in hg_lts:
				gn = lt.split('|')[0]
				if not gn in representative_genbanks: continue
				sample_copy_counts[gn] += 1
			if sum(sample_copy_counts.values()) == 0: continue
			ml_track_handle.write(hg + '\t' + str(hg_ordr) + '\t' + str(float(hg_mlen)/1000.0) + '\n')
			for gn in representative_genbanks:
				pres = '0'
				copy_count = ''
				if sample_copy_counts[gn] > 0:
					pres = '1'
					if sample_copy_counts[gn] > 1:
						copy_count = str(sample_copy_counts[gn])
				gn_label = gn
				if not full_genbank_labels:
					gn_label = gn
					if len(gn) >= 21:
						gn_label = gn[:20]
				gn_lab_keys.add(tuple([gn, gn_label]))
				gn_labs.add(gn_label)
				hm_track_handle.write(hg + '\t' + str(hg_ordr) + '\t' + gn_label + '\t' + pres + '\t' + str(copy_count) + '\n')
		ml_track_handle.close()
		hm_track_handle.close()

		try:
			assert(len(gn_labs) == len(gn_lab_keys))
		except:
			logObject.info('Non-unique labels resulted from truncating GenBank names. Please rerun zol with the "--full_genbank_labels" argument.')
			sys.stderr.write('Non-unique labels resulted from truncating GenBank names. Please rerun zol with the "--full_genbank_labels" argument.\n')
			sys.exit(1)

		plot_cmd = ['Rscript', plot_prog, ml_track_file, hm_track_file, plot_result_pdf, str(height), str(width)]
		try:
			subprocess.call(' '.join(plot_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(plot_result_pdf))
			logObject.info('Successfully ran: %s' % ' '.join(plot_cmd))
		except Exception as e:
			logObject.error('Had an issue running R based plotting - potentially because of R setup issues in conda: %s' % ' '.join(plot_cmd))
			sys.stderr.write('Had an issue running R based plotting - potentially because of R setup issues in conda: %s\n' % ' '.join(plot_cmd))
			logObject.error(e)
			sys.exit(1)
	except Exception as e:
		sys.stderr.write('Issues creating visualizations.\n')
		logObject.error('Issues creating visualizations.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)