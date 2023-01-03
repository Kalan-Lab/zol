import os
import sys
from Bio import SeqIO
import subprocess
from collections import defaultdict
import multiprocessing
from zol import util
import itertools
import math

zol_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'

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

def createGeneTrees(prot_algn_trim_dir, tree_dir, logObject, cpus=1):
	try:
		fasttree_cmds = []
		for patf in os.listdir(prot_algn_trim_dir):
			prefix = '.msa.faa'.join(patf.split('.msa.faa')[:-1])
			prot_algn_trim_file = prot_algn_trim_dir + patf
			tree_file = tree_dir + prefix + '.tre'
			fasttree_cmds.append(['fasttree', prot_algn_trim_file, '>', tree_file, logObject])
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

def annotateConsensusSequences(protein_faa, annotation_dir, logObject, cpus=1):
	"""
	pfam	NA	/home/rauf/Projects/KalanLab/Develop_zol/coding/zol/db/Pfam-A.hmm
	ko	/home/rauf/Projects/KalanLab/Develop_zol/coding/zol/db/ko_list	/home/rauf/Projects/KalanLab/Develop_zol/coding/zol/db/profile.hmm
	pgap	/home/rauf/Projects/KalanLab/Develop_zol/coding/zol/db/hmm_PGAP.tsv	/home/rauf/Projects/KalanLab/Develop_zol/coding/zol/db/PGAP.hmm
	mibig	NA	/home/rauf/Projects/KalanLab/Develop_zol/coding/zol/db/mibig.dmnd
	vfdb	NA	/home/rauf/Projects/KalanLab/Develop_zol/coding/zol/db/vfdb.dmnd
	card	NA	/home/rauf/Projects/KalanLab/Develop_zol/coding/zol/db/card.dmnd
	vog	/home/rauf/Projects/KalanLab/Develop_zol/coding/zol/db/vog.annotations.tsv	/home/rauf/Projects/KalanLab/Develop_zol/coding/zol/db/vog.dmnd
	paperblast	NA	/home/rauf/Projects/KalanLab/Develop_zol/coding/zol/db/paperblast.dmnd
	isfinder	NA	/home/rauf/Projects/KalanLab/Develop_zol/coding/zol/db/isfinder.dmnd
	"""
	db_locations = zol_main_directory + 'db/database_location_paths.txt'
	try:
		individual_cpus = 1
		pool_size = cpus
		if cpus > 9:
			individual_cpus = math.floor(cpus/pool_size)
			pool_size = 9
		assert(os.path.isfile(db_locations))
		search_cmds = []
		annotation_descriptions = defaultdict(lambda: defaultdict(lambda: 'NA'))
		with open(db_locations) as odls:
			for line in odls:
				line = line.strip()
				name, annot_info_file, db_file = line.split('\t')

				if annot_info_file != 'NA':
					if name == 'ko':
						with open(annot_info_file) as oaif:
							for i, line in enumerate(oaif):
								if i == 0: continue
								line = line.strip()
								ls = line.split('\t')
								ko = ls[0]
								if ls[1] == '-': continue
								description = ls[-1]
								annotation_descriptions['ko'][ko] = description
					elif name == 'pgap':
						with open(annot_info_file) as opil:
							for i, line in enumerate(opil):
								if i == 0: continue
								line = line.strip()
								ls = line.split('\t')
								label = ls[2]
								description = ls[10]
								annotation_descriptions['pgap'][label] = description
					elif name == 'vog':

				annotation_result_file = annotation_dir + name + '.txt'
				if db_file.endswith('.hmm'):
					search_cmd = ['hmmsearch', '--cpu', str(individual_cpus), '--tblout', annotation_result_file,
								  db_file, protein_faa, logObject]
				elif db_file.endswith('.dmnd'):
					search_cmd = ['diamond', 'blastp', '-d', db_file, '-q', protein_faa, '-o', annotation_result_file,
								  logObject]
				search_cmds.append(search_cmd)

		p = multiprocessing.Pool(pool_size)
		p.map(util.multiProcess, search_cmds)
		p.close()

		for rf in os.listdir(annotation_dir):
			with open(annotation_dir) +

	except Exception as e:
		sys.stderr.write('Issues with creating profile HMMs and consensus sequences.\n')
		logObject.error('Issues with creating profile HMMs and consensus sequences.')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

def runFubarAnalysis(cod):
	pass
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




