import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import subprocess
from collections import defaultdict
import multiprocessing
import tqdm
from zol import util, data_dictionary
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
import pyhmmer
import itertools
import csv

def determine_ranges(i):
	"""
	Deteremine continuous ranges from a list of integers. Answer taken from the answers by
	user97370 and bossylobster on StackOverflow: 
	https://stackoverflow.com/questions/4628333/converting-a-list-of-integers-into-range-in-python
	"""
	for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
		b = list(b)
		yield b[0][1], b[-1][1]

def split_by_idx(S, list_of_indices):
        """
        Function taken from https://stackoverflow.com/questions/10851445/splitting-a-string-by-list-of-indices
        """
        left, right = 0, list_of_indices[0]
        yield S[left:right]
        left = right
        for right in list_of_indices[1:]:
                yield S[left:right]
                left = right
        yield S[left:]

def performSLFullProteinClustering(ortho_matrix_file, workspace_dir, full_prot_cluster_file, logObject):
	"""
	Description:
	Runs single-linkage clustering of full proteins based on domain/intra-domain sequence orthogroups.
	"""
	try:
		observed_fp_pairs = set([])
		observed_fp = set([])
		pair_file = workspace_dir + 'Full_Protein_Pairs.txt'				
		opf_handle = open(pair_file, 'w')
		with open(ortho_matrix_file) as oomf:
			for i, line in enumerate(oomf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0: continue
				else:
					og = ls[0]
					flat_fps = []
					for lts in ls[1:]:
						lts = lts.strip()
						if lts == '': continue
						for lt in lts.split(', '):
							fp = '|'.join(lt.split('|')[:-2])
							flat_fps.append(fp)
							observed_fp.add(fp)
					for j, fp1 in enumerate(sorted(flat_fps)):
						for k, fp2 in enumerate(sorted(flat_fps)):
							if j <= k: continue
							pair_tuple = tuple(sorted([fp1, fp2]))
							if pair_tuple in observed_fp_pairs: continue
							opf_handle.write(fp1 + '\t' + fp2 + '\n')
							observed_fp_pairs.add(pair_tuple)
		opf_handle.close()

		clusters_file = workspace_dir + 'Full_Protein_Clusters.txt'
		clust_cmd = ['slclust', '<', pair_file, '>', clusters_file]

		try:
			subprocess.call(' '.join(clust_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
			assert (os.path.isfile(clusters_file))
			logObject.info('Successfully ran: %s' % ' '.join(clust_cmd))
		except Exception as e:
			logObject.error('Had an issue running slclust: %s' % ' '.join(clust_cmd))
			sys.stderr.write('Had an issue running slclust: %s\n' % ' '.join(clust_cmd))
			logObject.error(e)
			sys.exit(1)
		
		cluster_id = 1
		res_handle = open(full_prot_cluster_file, 'w')
		clustered_set = set([])
		with open(clusters_file) as ocf:
			for line in ocf:
				line = line.strip()
				ls = line.split()
				for prot in ls:
					clustered_set.add(prot)
					res_handle.write(prot + '\t' + 'SL_' + str(cluster_id) + '\n')
				cluster_id += 1
		
		for prot in observed_fp:
			if not prot in clustered_set:
				res_handle.write(prot + '\t' + 'SL_' + str(cluster_id) + '\n')
				cluster_id += 1
		res_handle.close()

	except:
		sys.stderr.write(traceback.format_exc() + '\n')
		logObject.error(traceback.format_exc())
		sys.exit(1)

def mapChunkProteinCoordsToFeatureCoords(start_coord, end_coord, tg_seq_chunk, tg_coord_info, name, evalue, eukaryotic_gene_cluster_flag, logObject):
	"""
	Description:
	Determines coordinates of chopped protein region from feature nucleotide coordinate info - which could include 
	skipped regions corresponding to introns.
	********************************************************************************************************************
	Parameters:
	- start_coord: The start coordinate of the chunk along the protein sequence.
	- end_coord: The end coordinate of the chuk along the protein sequence.
	- tg_seq_chunk: The protein chunk - for validation
	- tg_coord_info: Coordinate info for the full feature from the GenBank file.
	- name: The unique domain identifier.
	- evalue: The E-value from pyhmmer for the domain hit. NA if not a domain hit but inter-domain region or full protein.
	- logObject: A logging object.
	********************************************************************************************************************
	"""
	
	try:
		#print(name)
		#print(evalue)
		start_coord = start_coord-1
		nucl_seq, all_coords, start, end, direction, is_multi_part = tg_coord_info
		chunk_nucl_seq = ''
		nucl_start_coord = start_coord*3
		nucl_end_coord = end_coord*3
		chunk_coords = []
		nucl_coord = 0
	
		if direction == '+':	
			for sc, ec, dc in sorted(all_coords, key=itemgetter(0), reverse=False):
				for pos in range(sc, ec+1):
					if nucl_coord >= nucl_start_coord and nucl_coord < nucl_end_coord:
						chunk_nucl_seq += nucl_seq[nucl_coord]
						chunk_coords.append(pos)
					nucl_coord += 1
		else:
			for sc, ec, dc in sorted(all_coords, key=itemgetter(0), reverse=True):
				coord_list = list(set(range(sc, ec+1)))
				for pos in sorted(coord_list, reverse=True):
					if nucl_coord >= nucl_start_coord and nucl_coord < nucl_end_coord:
						chunk_nucl_seq += nucl_seq[nucl_coord]
						chunk_coords.append(pos)
					nucl_coord += 1		
		translated_prot_seq = str(Seq(chunk_nucl_seq).translate())
		
		"""
		useful code for debugging
		if not len(chunk_nucl_seq)%3 == 0:
			print(start)
			print(start_coord+1)

			print(end)
			print(end_coord)
			print(nucl_start_coord)
			print(nucl_end_coord)

			print(direction)
			print('------')
			print(nucl_seq)
			print('-----')
			print(chunk_nucl_seq)
			print(translated_prot_seq)
			print(tg_seq_chunk)
			sys.exit(1)
		"""
		
		# account for alterante initiator codons "GTG" and "TTG"
		if not eukaryotic_gene_cluster_flag and start_coord == 0 and translated_prot_seq[0] in set(['V', 'L']):
			translated_prot_seq = 'M' + translated_prot_seq[1:]
			msg = "Warning: changing starting protein residue from V or L to M. This is expected for bacteria, but if you are seeing this with running a eukaryotic dataset - it indicates an issue!"
			logObject.warning(msg)
			sys.stderr.write(msg + '\n')

		#print(name)
		#print(direction)
		#print(translated_prot_seq)
		#print(tg_seq_chunk)
		assert(translated_prot_seq == tg_seq_chunk)

		fstrand = 1
		if direction == '-':
			fstrand = -1

		range_coords = determine_ranges(sorted(chunk_coords))
		coord_feat_locs = []
		
		for rc in range_coords:
			coord_feat_locs.append(FeatureLocation(rc[0], rc[1], strand=fstrand))
	
		summed_coord_feat_locs = sum(coord_feat_locs)
		if direction == '-':
			summed_coord_feat_locs = sum(coord_feat_locs[::-1])
		feature = SeqFeature(summed_coord_feat_locs, type='cCDS')

		feature.qualifiers['translation'] = Seq(tg_seq_chunk)
		feature.qualifiers['paf_nucl_seq'] = chunk_nucl_seq
		feature.qualifiers['locus_tag'] = name
		feature.qualifiers['ccds_pyhmmer_evalue'] = evalue
		return(feature)
	except:
		sys.stderr.write(traceback.format_exc() + '\n')
		logObject.error(traceback.format_exc())
		raise RuntimeError()

def createChoppedGenbank(inputs):
	"""
	Description:
	Create a chopped CDS GenBank file from a regular GenBank file - core function for batchCreateChoppedGenbanks().
	********************************************************************************************************************
	Parameters:
	- inputs:
		- gbk: Input GenBank file.
		- prot_file: The protbin/cgceome file to temporary create to perfrom domain annotation using pyhmmer.
		- mapping_file: The mapping file between 
		- ccds_gbk_file: The chopped up GenBank file.
		- pfam_db_file: The Pfam HMM DB file.
		- pfam_z: The Pfam record count - for accurate E-value estimation.
		- minimal_length: The minimum length in amino acids for a domain matching or intra-domain region to be kept and
	                      tagged as a chopped CDS feature.
		- logObject: A logging object.
	- threads: The number of threads to use [Default is 1].
	********************************************************************************************************************
	"""
	gbk, prot_file, mapping_file, ccds_gbk_file, pfam_db_file, pfam_z, minimal_length, eukaryotic_gene_cluster_flag, logObject = inputs
	try:
		lt_coord_info = {}
		pf_handle = open(prot_file, 'w')
		feat_record = {}
		with open(gbk) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				full_sequence = str(rec.seq)
				for feature in rec.features:
					if feature.type == 'CDS':
						lt = feature.qualifiers.get('locus_tag')[0]
						seq = feature.qualifiers.get('translation')[0]
						nucl_seq = None
						all_coords, start, end, direction = util.processLocationString(str(feature.location))

						nucl_seq = ''
						for sc, ec, dc in sorted(all_coords, key=itemgetter(0), reverse=False):
							if ec >= len(full_sequence):
								nucl_seq += full_sequence[sc - 1:]
							else:
								nucl_seq += full_sequence[sc - 1:ec]					
						if direction == '-':
							nucl_seq = str(Seq(nucl_seq).reverse_complement())
							
						feat_record[lt] = rec.id
						lt_coord_info[lt] = [nucl_seq, all_coords, start, end, direction, is_multi_part]
						pf_handle.write('>' + lt + '\n' + str(seq) + '\n')
		pf_handle.close()
				
		# align Pfam domains and remove overlap similar to BiG-SCAPE
		alphabet = pyhmmer.easel.Alphabet.amino()
		sequences = []
		with pyhmmer.easel.SequenceFile(prot_file, digital=True, alphabet=alphabet) as seq_file:
			sequences = list(seq_file)

		target_dom_hits = defaultdict(list)
		with pyhmmer.plan7.HMMFile(pfam_db_file) as hmm_file:
			for hits in pyhmmer.hmmsearch(hmm_file, sequences, bit_cutoffs="trusted", Z=int(pfam_z), cpus=1):
				for hit in hits:
					for domain in hit.domains.included:
						target_dom_hits[hit.name.decode()].append([hits.query.name.decode(), domain.alignment.target_from, domain.alignment.target_to, domain.score, domain.i_evalue])

		# chop up FASTA based on mostly non-overlapping domains, 10% leaway is given
		map_handle = open(mapping_file, 'w')
		breakpoints = defaultdict(list)
		dom_start_names = defaultdict(lambda: 'NA')
		dom_start_evalues = defaultdict(lambda: 'NA')
		for tg in target_dom_hits:
			tg_dom_name_iter = defaultdict(int)
			accounted_coords = set([])
			for dom_align_info in sorted(target_dom_hits[tg], key=itemgetter(3), reverse=True):
				dom_name, start, end, score, i_evalue = dom_align_info
				overlap_coords = accounted_coords.intersection(set(range(start, end+1)))
				if (len(overlap_coords)/float(end-start+1) >= 0.1) or (len(overlap_coords) >= minimal_length): continue
				accounted_coords = accounted_coords.union(set(range(start, end+1)))
				breakpoints[tg].append(start)
				breakpoints[tg].append(end+1)
				map_handle.write('\t'.join([str(x) for x in [tg, start, end+1, dom_name, score, i_evalue]]) + '\n')		
				dom_start_names[tg + '_' + str(start)] = tg + '|' + dom_name + '|' + str(tg_dom_name_iter[dom_name]+1) 
				tg_dom_name_iter[dom_name] += 1
				dom_start_evalues[tg + '_' + str(start)] = str(i_evalue)
		map_handle.close()

		chopped_features = defaultdict(list)
		with open(prot_file) as ocf:
			for rec in SeqIO.parse(ocf, 'fasta'):
				tg = rec.id
				record = feat_record[tg]
				tg_seq = str(rec.seq)
				prev_end_coord = 1
				tg_interdomain_index = 1
				if not tg in breakpoints and len(tg_seq) >= minimal_length:
					dn = tg + '|full_protein|1'
					de = 'NA'
					chopped_feature = mapChunkProteinCoordsToFeatureCoords(1, len(tg_seq), tg_seq, lt_coord_info[tg], dn, de, eukaryotic_gene_cluster_flag, logObject)
					chopped_features[record].append(chopped_feature)
				else:
					for tg_seq_chunk in split_by_idx(tg_seq, ([0] + sorted(breakpoints[tg]))):
						if tg_seq_chunk.strip() == '': continue
						end_coord = prev_end_coord + len(tg_seq_chunk) - 1
						# limit chosen arbitrarily because avg. alpha-helix is ~10 residues and beta-strand 4-10 residues
						if len(tg_seq_chunk) >= minimal_length:
							dn = dom_start_names[tg + '_' + str(prev_end_coord-1)]
							de = dom_start_evalues[tg + '_' + str(prev_end_coord-1)]
							if dn == 'NA':
								dn = tg + '|inter-domain_region|' + str(tg_interdomain_index)
								tg_interdomain_index += 1
							chopped_feature = mapChunkProteinCoordsToFeatureCoords(prev_end_coord, end_coord, tg_seq_chunk, lt_coord_info[tg], dn, de, eukaryotic_gene_cluster_flag, logObject)
							chopped_features[record].append(chopped_feature)

						prev_end_coord = end_coord + 1

		occgf = open(ccds_gbk_file, 'w')
		with open(gbk) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				seq = Seq(rec.seq)
				gbk_rec = SeqRecord(seq, id=rec.id, name=rec.id, description=rec.description)
				gbk_rec.annotations['molecule_type'] = 'DNA'
				gbk_rec.features = chopped_features[rec.id]
				SeqIO.write(gbk_rec, occgf, 'genbank')
		occgf.close()

	except:
		msg = 'An issue occurred with creating chopped up version of GenBank file %s.' % gbk
		logObject.warning(msg)
		sys.stderr.write(traceback.format_exc())
		sys.stderr.write(msg + '\n')
		sys.exit(1)

def batchCreateChoppedGenbanks(genbanks, minimal_length, dm_scratch_dir, modified_genbank_dir, dom_to_cds_relations_file, logObject, eukaryotic_gene_cluster_flag=False, threads=1):
	"""	
	Description:
	Create chopped CDS GenBank files from regular GenBank input with CDS features.
	********************************************************************************************************************
	Parameters:
	- genbanks: Input GenBank files.
	- minimal_length: The minimum length in amino acids for a domain matching or intra-domain region to be kept and
	                  tagged as a chopped CDS feature.
	- 
	- modified_genbank_dir: The output directory where to write chopped CDS GenBank files. 
	- dom_to_cds_relations_file: The chopped CDS to full CDS relationship file.
	- logObject: A logging object.
	- eukaryotic_gene_cluster_flag: Flag for weather eukaryotic gene clusters are being investigated.
	- threads: The number of threads to use [Default is 1].
	********************************************************************************************************************
	Returns:
	- chopped_genbanks: A python list containing paths to chopped CDS GenBank files.
	********************************************************************************************************************
	"""

	chopped_genbanks = []
	try:
		pfam_db_file = None
		z = None
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
					name, _, db_file, z = line.split('\t')
					if name == 'pfam': 
						pfam_db_file = db_file
						pfam_z = int(z)
	
			assert(pfam_db_file != None and os.path.isfile(pfam_db_file) and z != None)
		except:
			msg = 'Issues validating that the Pfam database - please run setup_annotation_database.py to setup databases - note you might be interested in downloading the minimal database instead of the full one.'
			logObject.warning(msg)
			sys.stderr.write(msg + '\n')

		
		proteome_dir = dm_scratch_dir + 'GenBank_Proteins/'
		mapping_dir = dm_scratch_dir + 'cCDS_to_CDS_Mappings/'
		util.setupReadyDirectory([proteome_dir, mapping_dir])

		gbk_mod_inputs = []
		for gbk in genbanks:
			ccds_gbk_file = modified_genbank_dir + gbk.split('/')[-1]
			prot_file = proteome_dir + '.'.join(gbk.split('/')[-1].split('.')[:-1]) + '.faa'
			mapping_file = mapping_dir + '.'.join(gbk.split('/')[-1].split('.')[:-1]) + '.txt'
			gbk_mod_inputs.append([gbk, prot_file, mapping_file, ccds_gbk_file, pfam_db_file, pfam_z, minimal_length, eukaryotic_gene_cluster_flag, logObject])

		msg = "Creating domain-chopped up version of GenBank files for %d gene clusters" % len(gbk_mod_inputs) 
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(createChoppedGenbank, gbk_mod_inputs), total=len(gbk_mod_inputs)):
			pass
		p.close()

		ccc_handle = open(dom_to_cds_relations_file, 'w')
		for f in os.listdir(modified_genbank_dir):
			gbk_file = modified_genbank_dir + f
			chopped_genbanks.append(gbk_file)
			with open(gbk_file) as ogf:
				for rec in SeqIO.parse(ogf, 'genbank'):
					for feat in rec.features:
						if feat.type != 'cCDS': continue
						name = feat.qualifiers.get('locus_tag')[0]
						tg = name.split('|')[0]
						ccc_handle.write('\t'.join([name, tg, gbk_file]) + '\n')
		ccc_handle.close()
	
		return(chopped_genbanks)
	except:
		msg = 'Issue regenerating GenBank files with CDS features chopped up.'
		logObject.warning(msg)
		sys.stderr.write(traceback.format_exc())
		sys.stderr.write(msg + '\n')
		sys.exit(1) 

def reinflateOrthoGroups(ortho_matrix_file, prot_dir, rog_dir, logObject, cdhit_params="-c 0.98 -aL 0.95 -aS 0.95 -n 5 -M 4000", threads=1):
	"""
	Description:
	This function reinflates a matrix of ortholog groups to include all proteins in a given directory.
	The function first reads the ortholog group matrix and creates a set of protein IDs that are representatives of
	ortholog groups from the representative (dereplicated) set of gene clusters. The function then uses the CD-HIT
	program to cluster all proteins in the prot_dir directory and reads the CD-HIT clustering output to create a
	dictionary that maps non-representative protein IDs to ortholog groups.

	*******************************************************************************************************************
	Parameters:
	- orthogroup_matrix_file: The ortholog group vs sample matrix file, where cells correspond to locus tag identifiers.
	- prot_dir: A directory containing protein FASTA files.
	- rog_dir: A directory to write temporary + result files pertaining to reinflation to.
	- logObject: An object for logging messages.
	- cdhit_params: CD-HIT parameters. Default is: "-c 0.98 -aL 0.95 -aS 0.95 -n 5 -M 4000"
	- threads: The number of threads to use for the CD-HIT clustering step.
	*******************************************************************************************************************
	"""
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
		cdhit_cmd = ['cd-hit', '-i', comp_prot_file, '-o', cdhit_nr_prefix, cdhit_params,
					 '-d', '0', '-T', str(threads)]

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
							sys.stderr.write('Warning: The protein %s is a representative of ortholog group %s, but can potentially belong to multiple. Skipping its incorporation for ortholog group %s.\n' % (pid, rep_to_hg[pid], hg))
							logObject.warning('The protein %s is a representative of a ortholog group %s, but can potentially belong to multiple. Skipping its incorporation for ortholog group %s.' % (pid, rep_to_hg[pid], hg))
						elif pid in accounted:
							sys.stderr.write('Warning: The protein %s has already been clustered into a ortholog group, but can potentially belong to multiple. Skipping its incorporation for ortholog group %s.\n' % (pid, hg))
							logObject.warning('The protein %s has already been clustered into a ortholog group, but can potentially belong to multiple. Skipping its incorporation for ortholog group %s.' % (pid, hg))
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

def dereplicateUsingSkani(genbanks, focal_genbanks, derep_dir, kept_dir, logObject, skani_small_genomes_preset=False, 
						  skani_identiy_threshold=99.0, skani_coverage_threshold=95.0, mcl_inflation=None, threads=1):
	"""
	Description:
	This function dereplicates a set of GenBank files using the skani to estimate pairwise gene-cluster ANI and either
	single-linkage clustering (slclust) or MCL to cluster and select representative gene-clusters. If focal gene
	clusters are noted by the user, they are given priority as representatives otherwise representatives are chosen
	based on length (longest given priority).
	*******************************************************************************************************************
	Parameters:
	- genbanks: A list of paths to GenBank files.
	- focal_genbanks: A list of paths to GenBank files that should be included in the dereplicated set if possible,
	                even if they are not the largest sequences in their cluster.
	- derep_dir: The directory to write the dereplicated GenBank files to.
	- kept_dir: The directory to write the GenBank files that were kept after dereplication to.
	- logObject: An object for logging messages.
	- skani_small_genomes_preset: Use the --small-genomes preset in skani for faster computes.
	- skani_identiy_threshold: The minimum identity threshold for two sequences to be considered similar.
	- skani_coverage_threshold: The minimum coverage threshold for two sequences to be considered similar.
	- mcl_inflation: The inflation factor to use for the MCL clustering algorithm. If not provided (default), single-
	                 linkage clustering (via slclust) will be used instead.
	- threads: The number of threads to use for the skani and MCL clustering algorithms.
	*******************************************************************************************************************
	Returns:
	A tuple of two lists:
		- The first list contains the paths to the dereplicated GenBank files.
		- The second list contains the members of each cluster in the dereplicated set.
	*******************************************************************************************************************
	"""

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
		skani_sketch_cmd = ['skani', 'sketch', '-t', str(threads), '-l', fasta_listing_file, '-o', skani_sketch_db]
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
		skani_dist_cmd = ['skani', 'dist', '-t', str(threads), '-q', skani_sketch_db + '*', '-r', skani_sketch_db + '*',
						  '--min-af', str(skani_coverage_threshold), '-o', skani_result_file]
		if skani_small_genomes_preset:
			skani_dist_cmd += ['--small-genomes']

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
						 '-te', str(threads)]

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
	"""
	Description:
	This function partitions gene cluster gene/protein sequences by ortholog groups.
	*******************************************************************************************************************
	Parameters:
	- ortho_matrix_file: A file containing an orthomatrix.
	- prot_dir: A directory containing protein sequences.
	- nucl_dir: A directory containing DNA sequences.
	- hg_prot_dir: A directory to write protein sequences for each ortholog group.
	- hg_nucl_dir: A directory to write DNA sequences for each ortholog group.
	- logObject: A logging object.
	*******************************************************************************************************************
	"""
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
		sys.stderr.write('Issues with partitioning sequences to ortholog groups.\n')
		logObject.error('Issues with partitioning sequences to ortholog groups.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def partitionAndCreateUpstreamNuclAlignments(ortho_matrix_file, nucl_upstr_dir, hg_upst_dir, upst_algn_dir, logObject,
											 threads=1, use_super5=False):
	"""
	Description:
	This function partitions upstream DNA sequences into ortholog groups and creates alignments for each group.
	*******************************************************************************************************************
	Parameters:
	- ortho_matrix_file: A file containing an orthomatrix.
	- nucl_upstr_dir: A directory containing upstream DNA sequences.
	- hg_upst_dir: A directory to write upstream DNA sequences for each ortholog group.
	- upst_algn_dir: A directory to write alignments for each ortholog group.
	- logObject: A logging object.
	- threads: The number of threads to use for alignment.
	- use_super5: Whether to use the SUPER5 algorithm for MUSCLE alignment.
	"""
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
			align_cmd = ['muscle', '-align', upst_file, '-output', upst_algn_file, '-nt', '-threads', str(threads), '-perturb', '12345']
			if use_super5:
				align_cmd = ['muscle', '-super5', upst_file, '-output', upst_algn_file, '-nt', '-threads', str(threads), '-perturb', '12345']
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
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def createProteinAlignments(prot_dir, prot_algn_dir, logObject, use_super5=False, threads=1):
	"""
	Description:
	This function creates protein alignments from a directory of protein sequences.
	*******************************************************************************************************************
	Parameters:
	- prot_dir: A directory containing protein sequences.
	- prot_algn_dir: A directory to write protein alignments.
	- logObject: A logging object.
	- use_super5: Whether to use the SUPER5 algorithm for MUSCLE alignment.
	- threads: The number of threads to use for alignment.
	"""
	try:
		for pf in os.listdir(prot_dir):
			prefix = '.faa'.join(pf.split('.faa')[:-1])
			prot_file = prot_dir + pf
			prot_algn_file = prot_algn_dir + prefix + '.msa.faa'
			align_cmd = ['muscle', '-align', prot_file, '-output', prot_algn_file, '-amino', '-threads', str(threads), '-perturb', '12345']
			if use_super5:
				align_cmd = ['muscle', '-super5', prot_file, '-output', prot_algn_file, '-amino', '-threads', str(threads), '-perturb', '12345']
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
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def createCodonAlignments(prot_algn_dir, nucl_dir, codo_algn_dir, logObject, threads=1):
	"""
	Description:
	This function creates codon alignments from a directory of protein alignments and a directory of DNA sequences.
	*******************************************************************************************************************
	Parameters:
	- prot_algn_dir: A directory containing protein alignments.
	- nucl_dir: A directory containing DNA sequences.
	- codo_algn_dir: A directory to write codon alignments.
	- logObject: A logging object.
	- threads: The number of threads to use for alignment.
	"""
	try:
		pal2nal_cmds = []
		for paf in os.listdir(prot_algn_dir):
			prefix = '.msa.faa'.join(paf.split('.msa.faa')[:-1])
			prot_algn_file = prot_algn_dir + paf
			nucl_file = nucl_dir + prefix + '.fna'
			codo_algn_file = codo_algn_dir + prefix + '.msa.fna'
			pal2nal_cmds.append(['pal2nal.pl', prot_algn_file, nucl_file, '-output', 'fasta', '>', codo_algn_file, logObject])
		
		msg = "Running PAL2NAL to generate codon alignments for %d ortholog groups" % len(pal2nal_cmds) 
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(util.multiProcess, pal2nal_cmds), total=len(pal2nal_cmds)):
			pass
		p.close()

	except Exception as e:
		sys.stderr.write('Issues with creating codon alignments.\n')
		logObject.error('Issues with creating codon alignments.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def trimAlignments(prot_algn_dir, codo_algn_dir, prot_algn_trim_dir, codo_algn_trim_dir, logObject, threads=1):
	"""
	Description:
	This function trims protein and codon alignments using TrimAl.
	*******************************************************************************************************************
	Parameters:
	- prot_algn_dir: The directory containing the protein alignments.
	- codo_algn_dir: The directory containing the codon alignments.
	- prot_algn_trim_dir: The directory where the trimmed protein alignments will be saved.
	- codo_algn_trim_dir: The directory where the trimmed codon alignments will be saved.
	- logObject: A logging object.
	- threads: The number of threads to use for trimming the alignments.
	*******************************************************************************************************************
	"""
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

		msg = "Running trimal to generate trimmed protein and codon alignments for %d ortholog groups" % (len(trim_cmds)/2) 
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(util.multiProcess, trim_cmds), total=len(trim_cmds)):
			pass
		p.close()

	except Exception as e:
		sys.stderr.write('Issues with trimming protein/codon alignments.\n')
		logObject.error('Issues with trimming protein/codon alignments.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def createGeneTrees(codo_algn_trim_dir, codo_algn_dir, tree_dir, logObject, threads=1):
	"""
	Description:
	This function creates gene trees from trimmed codon alignments using FastTree2 for ortholog groups.
	*******************************************************************************************************************
	Parameters:
	- codo_algn_trim_dir: The directory containing trimmed codon alignments.
	- tree_dir: The directory where trees in Newick format will be saved.
	- logObject: A logging object.
	- threads: The number of threads to use.
	*******************************************************************************************************************
	"""
	try:
		fasttree_cmds = []
		for catf in os.listdir(codo_algn_dir):
			if not catf.endswith('.msa.fna'): continue
			prefix = '.msa.faa'.join(catf.split('.msa.fna')[:-1])
			codo_algn_trim_file = codo_algn_trim_dir + catf
			codo_algn_file = codo_algn_dir + catf
			seqlen = 0
			with open(codo_algn_trim_file) as ocatf:
				for rec in SeqIO.parse(ocatf, 'fasta'):
					seqlen = len(str(rec.seq))
			tree_file = tree_dir + prefix + '.tre'
			if seqlen > 0:
				fasttree_cmds.append(['fasttree', '-nt', codo_algn_trim_file, '>', tree_file, logObject])
			else:
				fasttree_cmds.append(['fasttree', '-nt', codo_algn_file, '>', tree_file, logObject])
				msg = 'Trimmed codon alignment was blank so defaulting to untrimmed for constructing ortholog group tree for %s.' % prefix
				logObject.warning(msg)
				sys.stderr.write('Warning: ' + msg + '\n')
				
		msg = "Running FastTree 2 to generate gene trees (based on trimmed codon alignments) for %d ortholog groups" % len(fasttree_cmds) 
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(util.multiProcess, fasttree_cmds), total=len(fasttree_cmds)):
			pass
		p.close()

	except Exception as e:
		sys.stderr.write('Issues with creating gene-trees.\n')
		logObject.error('Issues with creating gene-trees.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def createProfileHMMsAndConsensusSeqs(prot_algn_dir, phmm_dir, cons_dir, logObject, threads=1):
	"""
	Description:
	This function creates profile HMMs and emits consensus sequences based on protein MSAs using HMMER.
	*******************************************************************************************************************
	Parameters:
	- prot_algn_dir: The directory containing protein alignments.
	- phmm_dir: The directory where profile HMMs in HMMER3 HMM format will be saved.
	- cons_dir: The directory where consensus sequences in FASTA format will be saved.
	- logObject: A logging object.
	- threads: The number of threads to use.
	*******************************************************************************************************************
	"""
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

		msg = "Running HMMER3 hmmbuild to generate profile HMMs for %d ortholog groups" % len(hmmbuild_cmds)
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(util.multiProcess, hmmbuild_cmds), total=len(hmmbuild_cmds)):
			pass
		p.close()

		msg = "Running HMMER3 hmmemit to generate consensus protein sequences for %d ortholog groups" % len(hmmemit_cmds)
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(util.multiProcess, hmmemit_cmds), total=len(hmmemit_cmds)):
			pass
		p.close()

	except Exception as e:
		sys.stderr.write('Issues with creating profile HMMs and consensus sequences.\n')
		logObject.error('Issues with creating profile HMMs and consensus sequences.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def refineGeneCalling(custom_database, hg_prot_dir, hg_nucl_dir, refine_workpace_dir, logObject, use_super5=True, threads=1):
	"""
	Description:
	This function refines protein alignments by DIAMOND BLASTing proteins from a reference gene cluster and filtering
	sites which are not aligned. The intent of the function was to refine alignments when dealing with eukaryotic
	gene clusters and miniprot based mapping of gene calls in genomes had resulted in inclusion of intronic sequences.
	We have since improved the parsing of miniprot based gene mapping to make this function largely unecessary and might
	remove the functionality in future versions.
	*******************************************************************************************************************
	Parameters:
	- custom_database: Proteins in FASTA format from a reference gene cluster.
	- hg_prot_dir: The directory where the original (un-refined) ortholog group protein alignments are stored.
	- hg_nucl_dir: The directory where the origianl (un-refined) ortholog group nucleotide alignments are stored.
	- refine_workspace_dir: The workspace to perform the refinement of protein and nucleotide alignments.
	- logObject: A logging object.
	- use_super5: Whether to use SUPER5 algorithm for MUSCLE based alignment.
	- threads: The number of threads to use.
	*******************************************************************************************************************
	Returns:
	A list consisting of two items:
	- hg_prot_refined_dir (first item): The directory where the refined ortholog group protein alignments are stored.
	- hg_nucl_refined_dir (second item): The directory where the refined ortholog group protein alignments are stored.
	*******************************************************************************************************************
	"""
	hg_prot_refined_dir = refine_workpace_dir + 'OG_Protein_Refined_Sequences/'
	hg_nucl_refined_dir = refine_workpace_dir + 'OG_Nucleotide_Refined_Sequences/'
	try:
		tmp_hg_prot_seq_dir = refine_workpace_dir + 'OG_Prot_Seqs_with_References/'
		tmp_hg_prot_aln_dir = refine_workpace_dir + 'OG_Prot_Alns_with_References/'
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
		search_cmd = ['diamond', 'blastp', '--ignore-warnings', '-p', str(threads), '-d', dmnd_db, '-q', custom_database,
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
			sys.stderr.write(traceback.format_exc())
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
			sys.stderr.write(traceback.format_exc())
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

			align_cmd = ['muscle', '-align', hg_seq_with_que, '-output', hg_aln_with_que, '-amino', '-threads', str(threads), '-perturb', '12345']
			if use_super5:
				align_cmd = ['muscle', '-super5', hg_seq_with_que, '-output', hg_aln_with_que, '-amino', '-threads', str(threads), '-perturb', '12345']
			try:
				subprocess.call(' '.join(align_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
								executable='/bin/bash')
				assert(os.path.isfile(hg_aln_with_que))
				logObject.info('Successfully ran: %s' % ' '.join(align_cmd))
			except Exception as e:
				logObject.error('Had an issue running MUSCLE: %s' % ' '.join(align_cmd))
				sys.stderr.write('Had an issue running MUSCLE: %s\n' % ' '.join(align_cmd))
				logObject.error(e)
				sys.stderr.write(traceback.format_exc())
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

def annotateCustomDatabase(protein_faa, custom_protein_db_faa, annotation_dir, logObject, threads=1,
						   max_annotation_evalue=1e-5):
	"""
	Description:
	This function will annotate consensus sequences for ortholog groups with a custom database of protein sequences
	provided by the user in FASTA format. The best hit per ortholog group is selected based on bitscore provided that
	the E-value threshold is met.
	*******************************************************************************************************************
	Parameters:
	- protein_faa: Consensus protein sequences for ortholog groups in FASTA format.
	- custom_protein_db_faa: Custom database of reference proteins in FASTA format.
	- annotation_dir: Directory where to perform annotation analysis.
	- logObject: A logging object.
	- threads: The number of threads to use.
	- max_annotation_evalue: The maximum E-value by DIAMOND to regard an alignment between a consensus ortholog group
	                         sequence and a reference protein sequence.
	*******************************************************************************************************************
	Returns:
	A dictionary where the key in the ortholog group identifier and the value is a list of two lists, where the first
	list contains identifiers/descriptions of the custom/reference proteins and the second is a list of the respective
	E-values.
	*******************************************************************************************************************
	"""
	custom_annotations = {}
	try:
		custom_annot_dir = annotation_dir + 'Custom_Annotation/'
		util.setupReadyDirectory([custom_annot_dir])
		dmnd_db = custom_annot_dir + 'Custom.dmnd'
		blastp_file = custom_annot_dir + 'Custom.txt'
		makedb_cmd = ['diamond', 'makedb', '--ignore-warnings', '--in', custom_protein_db_faa, '-d', dmnd_db]
		search_cmd = ['diamond', 'blastp', '--ignore-warnings', '-p', str(threads), '-d', dmnd_db, '-q', protein_faa,
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
			sys.stderr.write(traceback.format_exc())
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
			sys.stderr.write(traceback.format_exc())
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

def runPyhmmer(inputs):
	name, db_file, z, protein_faa, annotation_result_file, threads = inputs
	try:
		alphabet = pyhmmer.easel.Alphabet.amino()
		sequences = []
		with pyhmmer.easel.SequenceFile(protein_faa, digital=True, alphabet=alphabet) as seq_file:
			sequences = list(seq_file)

		outf = open(annotation_result_file, 'w')
		if name == 'pfam':
			with pyhmmer.plan7.HMMFile(db_file) as hmm_file:
				for hits in pyhmmer.hmmsearch(hmm_file, sequences, bit_cutoffs='trusted', Z=int(z), cpus=threads):
					for hit in hits:
						accession = 'NA'
						try:
							accession = hits.query.accession.decode()
						except:
							accession = 'NA'
						outf.write('\t'.join([hits.query.name.decode(), accession, 'NA', hit.name.decode(), 'NA', str(hit.evalue), str(hit.score)]) + '\n')

		else:
			with pyhmmer.plan7.HMMFile(db_file) as hmm_file:
				for hits in pyhmmer.hmmsearch(hmm_file, sequences, Z=int(z), cpus=threads):
					for hit in hits:
						accession = 'NA'
						try:
							accession = hits.query.accession.decode()
						except:
							accession = 'NA'

						outf.write('\t'.join([hits.query.name.decode(), accession, 'NA', hit.name.decode(), 'NA', str(hit.evalue),
											  str(hit.score)]) + '\n')
		outf.close()
	except:
		raise RuntimeError('Problem running pyhmmer!')

def annotateConsensusSequences(protein_faa, annotation_dir, logObject, threads=1, max_annotation_evalue=1e-5):
	"""
	Description:
	This function will attempt to annotate consensus sequences for ortholog groups with the default databases supported
	within zol assuming they have been (properly) setup. The best hit per ortholog group per database is selected based
	on bitscore (if FASTA database) or score (if pHMM database) provided that the E-value threshold is met.
	*******************************************************************************************************************
	Parameters:
	- protein_faa: Consensus protein sequences for ortholog groups in FASTA format.
	- annotation_dir: Directory where to perform annotation analysis.
	- logObject: A logging object.
	- threads: The number of threads to use.
	- max_annotation_evalue: The maximum E-value by DIAMOND to regard an alignment between a consensus ortholog group
	                         sequence and a database sequence.
	*******************************************************************************************************************
	Returns:
	A dictionary of dictionaries where the primary key is the name of the database and the secondary key is the homolog
	group identifier and the values are a list of two lists, where the first list contains identifiers/descriptions of
	the best-hit database proteins and the second is a list of the respective E-values.
	*******************************************************************************************************************
	"""
	zol_data_directory = str(os.getenv("ZOL_DATA_PATH")).strip()
	db_locations = None
	if zol_data_directory != 'None':
		try:
			zol_data_directory = os.path.abspath(zol_data_directory) + '/'
			db_locations = zol_data_directory + 'database_location_paths.txt'
		except:
			pass
	if db_locations == None or not os.path.isfile(db_locations):
		sys.stderr.write('Warning: databases do not appear to be setup or setup properly!\n')
		annotations = defaultdict(lambda: defaultdict(lambda: ['NA', 'NA']))
		return (default_to_regular(annotations))

	try:
		dmnd_individual_threads = 1
		dmnd_pool_size = threads
		if threads > 5:
			dmnd_individual_threads = math.floor(threads/5)
			dmnd_pool_size = 5

		hmm_individual_threads = 1
		hmm_pool_size = threads
		if threads > 4:
			hmm_individual_threads = math.floor(threads/4)
			hmm_pool_size = 4

		assert(os.path.isfile(db_locations))
		dmnd_search_cmds = []
		hmm_search_cmds = []
		name_to_info_file = {}
		hmm_based_annotations = set([])
		with open(db_locations) as odls:
			for line in odls:
				line = line.strip()
				name, annot_info_file, db_file, z = line.split('\t')
				name_to_info_file[name] = annot_info_file
				annotation_result_file = annotation_dir + name + '.txt'
				if db_file.endswith('.hmm'):
					hmm_based_annotations.add(name)
					hmm_search_cmds.append([name, db_file, z, protein_faa, annotation_result_file, hmm_individual_threads])
				elif db_file.endswith('.dmnd') and not name in set(['riboprots', 'mobsuite']):
					search_cmd = ['diamond', 'blastp', '--ignore-warnings', '-p', str(dmnd_individual_threads), '-d', db_file,
								  '-q', protein_faa, '-o', annotation_result_file, logObject]
					dmnd_search_cmds.append(search_cmd)

		msg = "Running pyhmmer hmmsearch for functional annotation for %d databases" % len(hmm_search_cmds)
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(hmm_pool_size)
		for _ in tqdm.tqdm(p.imap_unordered(runPyhmmer, hmm_search_cmds), total=len(hmm_search_cmds)):
			pass
		p.close()

		msg = "Running DIAMOND blastp for functional annotation for %d databases" % len(dmnd_search_cmds)
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(dmnd_pool_size)
		for _ in tqdm.tqdm(p.imap_unordered(util.multiProcess, dmnd_search_cmds), total=len(dmnd_search_cmds)):
			pass
		p.close()

		annotations = defaultdict(lambda: defaultdict(lambda: ['NA', 'NA'])) # db -> query -> [hit descriptions, evalue]
		for rf in os.listdir(annotation_dir):
			if not rf.endswith('.txt'): continue
			db_name = rf.split('.txt')[0]
			if db_name in set(['riboprots', 'mobsuite']): continue
			annot_info_file = name_to_info_file[db_name]

			id_to_description = defaultdict(lambda: 'NA')
			with open(annot_info_file) as oaif:
				for line in oaif:
					line = line.strip()
					ls = line.split('\t')
					id_to_description[ls[0]] = ls[1]

			# or by_score if HMMscan - lets avoid a second variable
			best_hits_by_bitscore = defaultdict(lambda: [[], [], 0.0])
			# note second listing here is evalues for most annotation dbs but for Pfam it is accessions (the PF ids)
			if db_name in hmm_based_annotations:
				# parse HMM based results from pyhmmer
				with open(annotation_dir + rf) as oarf:
					for line in oarf:
						line = line.rstrip('\n')
						if line.startswith('#'): continue
						ls = line.split()
						query = ls[3]
						accession = ls[1]
						hit = ls[0]
						evalue = decimal.Decimal(ls[5])
						score = float(ls[6])
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
								best_hits_by_bitscore[query][1].append(accession)
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
	Convert a defaultdict to a regular old dict.

	Function taken from:
	https://stackoverflow.com/questions/26496831/how-to-convert-defaultdict-of-defaultdicts-of-defaultdicts-to-dict-of-dicts-o
	"""
	try:
		if isinstance(d, defaultdict):
			d = {k: default_to_regular(v) for k, v in d.items()}
		return d
	except:
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def determineConsensusOrderOfHGs(genbanks, ortho_matrix_file, markovian_file, consensus_path_file, logObject, domain_mode=False):
	"""
	Description:
	This function determines the consensus order and directionality of ortholog groups in a set of gene cluster
	GenBanks. It is closely based on code from lsaBGC.
	*******************************************************************************************************************
	Parameters:
	- genbanks: A list of gene cluster GenBank files.
	- orthogroup_matrix_file: The ortholog group vs sample matrix file, where cells correspond to locus tag identifiers.
	- logObject: A logging object.
	*******************************************************************************************************************
	Results:
	- hg_order_scores: A dictionary where keys are ortholog group identifiers and the values are a list of two items:
	                   (1) consensus order and (2) consensus direction.
	*******************************************************************************************************************
	"""
	try:
		gc_gene_to_hg = {}
		most_conserved_hgs = set([])
		hg_conservation_values = []
		max_hg_conservation = 0.0
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
				hg_conservation = sample_count/float(len(ls[1:])) 
				hg_conservation_values.append([hg, hg_conservation])
				if hg_conservation >= max_hg_conservation:
					max_hg_conservation = hg_conservation
				if hg_conservation == 1.0:
					core_hgs.add(hg)
					if sc_sample_count/float(len(ls[1:])) == 1.0:
						single_copy_core_hgs.add(hg)

		for hgv in hg_conservation_values:
			if hgv[1] == max_hg_conservation:
				most_conserved_hgs.add(hgv[0])

		gc_gene_counts = defaultdict(int)
		gc_genes = defaultdict(set)
		scaff_relative_gene_locations = {}
		for gbk in genbanks:
			prefix = '.'.join(gbk.split('/')[-1].split('.')[:-1])
			if domain_mode:
				gene_locations = util.parseGbk(gbk, prefix, logObject, feature_type='cCDS')
			else:
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
		hg_pair_score = defaultdict(int)
		hg_preceding_score = defaultdict(lambda: defaultdict(int))
		hg_following_score = defaultdict(lambda: defaultdict(int))
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
					hg_preceding_score[hg][hg_previ] += 1
					hg_following_score[hg_previ][hg] += 1
					hg_pair_score[tuple([hg_previ, hg])] += 1
				try:
					hg_after = hgs[j + 1]
					# make sure you don't get lost with broken/fragmented genes in BGCs that might be
					# in the process being lost.
					if hg != hg_after:
						hg_preceding_score[hg_after][hg] += 1
						hg_following_score[hg][hg_after] += 1
						hg_pair_score[tuple([hg, hg_after])] += 1
				except:
					hg_after = 'end'
					hg_preceding_score[hg_after][hg] += 1
					hg_following_score[hg][hg_after] += 1
					hg_pair_score[tuple([hg, hg_after])] += 1

		markovian_handle = open(markovian_file, 'w')
		markovian_handle.write('og\tog_after\tog_direction\tog_after_direction\tsupport\n')
		for hg in hg_following_score:
			for hg_after in hg_following_score[hg]:
				hg_dir = '-'
				if direction_forward_support[hg] >= direction_reverse_support[hg]: 
					hg_dir = '+'
				hg_after_dir = '-'
				if direction_forward_support[hg_after] >= direction_reverse_support[hg_after]:
					hg_after_dir = '+'
				markovian_handle.write(hg + '\t' + hg_after + '\t' + str(hg_following_score[hg][hg_after]) + '\t' + hg_dir + '\t' + hg_after_dir + '\n')
		markovian_handle.close()

		consensus_handle = open(consensus_path_file, 'w')		
		curr_og = 'start'
		visited_ogs = set([])
		while curr_og != 'end':
			next_og = ""
			max_score = 0 
			for follow_og in hg_following_score[curr_og]:
				if hg_following_score[curr_og][follow_og] >= max_score:
					max_score = hg_following_score[curr_og][follow_og]
					next_og = follow_og
			consensus_handle.write(curr_og + '\t' + next_og + '\n')
			if next_og in visited_ogs:
				break
			visited_ogs.add(curr_og)
			curr_og = next_og
		consensus_handle.close()

		anchor_edge = None
		if len(single_copy_core_hgs) > 0:
			# first attempt to find anchor edge using a single copy core ortholog groups if any exist
			for hps in sorted(hg_pair_score.items(), key=itemgetter(1), reverse=True):
				if hps[0][0] in single_copy_core_hgs and hps[0][1] in single_copy_core_hgs:
					anchor_edge = hps[0]
					break
		if len(core_hgs) > 0 and anchor_edge == None:
			# looks like that failed, now lets use any available core ortholog groups (not necessarily single copy) if any exist 
			for hps in sorted(hg_pair_score.items(), key=itemgetter(1), reverse=True):
				if hps[0][0] in core_hgs and hps[0][1] in core_hgs:
					anchor_edge = hps[0]
					break
			try:
				assert (anchor_edge != None)
			except:
				for hps in sorted(hg_pair_score.items(), key=itemgetter(1), reverse=True):
					if hps[0][0] in core_hgs or hps[0][1] in core_hgs:
						anchor_edge = hps[0]
						break

		if anchor_edge == None:
			# ahh, that also failed welp - lets use the most conserved gene available and write a warning to the log file and console
			stars = '*'*34 + '\n'
			sys.stderr.write(stars + 'WARNING!!! No core ortholog groups were detected across homologous gene cluster\ninstances - the consensus order and direction predictions will likely be lower quality.\n' + stars)
			logObject.warning('No core ortholog groups were detected across homologous gene cluster\ninstances - the quality of the consensus order and direction\npredictions will be lower.\n')
			try:
				for hps in sorted(hg_pair_score.items(), key=itemgetter(1), reverse=True):
					if hps[0][0] in most_conserved_hgs and hps[0][1] in most_conserved_hgs:
						anchor_edge = hps[0]
						break
				try:
					assert (anchor_edge != None)
				except:
					for hps in sorted(hg_pair_score.items(), key=itemgetter(1), reverse=True):
						if hps[0][0] in most_conserved_hgs or hps[0][1] in most_conserved_hgs:
							anchor_edge = hps[0]
							break
				assert(anchor_edge != None)
			except:
				sys.stderr.write(traceback.format_exc())
				sys.stderr.write("\nUnexpected error, no anchor edge found, could be because no protocore ortholog group exists, which shouldn't be the case!\n")
				sys.exit(1)


		# use to keep track of which HGs have been accounted for already at different steps of assigning order
		accounted_hgs = set([anchor_edge[0], anchor_edge[1]])

		# primary expansion left
		curr_hg = anchor_edge[0]
		left_expansion = [curr_hg]
		while not curr_hg == 'start':
			new_hg = None
			for i, hg in enumerate(sorted(hg_preceding_score[curr_hg].items(), key=itemgetter(1), reverse=True)):
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
			for i, hg in enumerate(sorted(hg_following_score[curr_hg].items(), key=itemgetter(1), reverse=True)):
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
				for phg in sorted(hg_preceding_score[hg].items(), key=itemgetter(1), reverse=True):
					if best_score < phg[1] and phg[0] in accounted_hgs:
						best_score = phg[1]
						relative_pos = 'after'
						neighboriest_hg = phg[0]
						break
				for fhg in sorted(hg_following_score[hg].items(), key=itemgetter(1), reverse=True):
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
				consensus_direction = '-'
				if direction_forward_support[hg] >= direction_reverse_support[hg]: consensus_direction = '+'
				hg_order_scores[hg] = [i, consensus_direction]
				i += 1
		return hg_order_scores
	except Exception as e:
		sys.stderr.write('Issues in attempting to calculate order score for each ortholog group.\n')
		logObject.error("Issues in attempting to calculate order score for each ortholog group.")
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def determineHGStats(orthogroup_matrix_file, hg_nucl_dir, logObject, representative_associated_members=None,
					 impute_broad_conservation=False):
	"""
	Description:
	This function determines basic statistics for ortholog groups, including: (i) percentage of gene clusters which
	feature theme, (ii) the median GC%, (iii) the median GC skew, (iv) whether they are found at max a copy-count of
	one per gene cluster, and (v) the median legnth of nucleotide sequences.
	*******************************************************************************************************************
	Parameters:
	- orthogroup_matrix_file: The ortholog group vs sample matrix file, where cells correspond to locus tag identifiers.
	- hg_nucl_dir: The directory of FASTA files with nucleotide sequence for each ortholog group.
	- logObject: A logging object.
	- representative_associated_members:
	- impute_broad_conservation:
	*******************************************************************************************************************
	Results:
	- hg_order_scores: a dictionary where the key is the ortholog group identifier
	*******************************************************************************************************************
	"""
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
		logObject.error('Issues with determining basic stats for ortholog groups.')
		sys.stderr.write('Issues with determining basic stats for ortholog groups.\n')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def individualHyphyRun(inputs):
	"""
	Description:
	This functions run HyPhy based analyses (GARD + FUBAR) for a single ortholog group.
	*******************************************************************************************************************
	Parameters:
	- inputs: a list which can be expanded to the following items:
		- hg: ortholog group identifier.
		- hg_codo_algn_file: ortholog group codon alignment file in FASTA format.
		- hg_codo_tree_file: ortholog group (approximate) phylogeny file in Newick format.
		- gard_output: output *.json file from GARD analysis.
		- best_gard_output: output *.best file from GARD analysis.
		- fubar_outdir: results directory for FUBAR analysis.
		- skip_gard: boolean flag for whether to skip GARD analysis.
		- skip_busted: boolean flag for whether to skip BUSTED analysis.
		- gard_mode: analysis mode for GARD - either "Faster" or "Normal".
		- logObject: a logging object.
	*******************************************************************************************************************
	"""
	hg, hg_codo_algn_file, hg_full_codo_tree_file, gard_output, best_gard_output, fubar_outdir, busted_outdir, skip_gard, skip_busted, gard_mode, gard_timeout, logObject = inputs
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
				sys.stderr.write(traceback.format_exc())
				sys.exit(1)
			if not skip_busted:
				busted_cmd = ['hyphy', 'CPU=1', 'busted', '--alignment', hg_codo_algn_file, '--tree', hg_full_codo_tree_file]
				try:
					subprocess.call(' '.join(busted_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
									executable='/bin/bash')
					assert(os.path.isfile(hg_codo_algn_file + '.BUSTED.json'))
					os.system('mv %s %s' % (hg_codo_algn_file + '.BUSTED.json', busted_outdir))
					logObject.info('Successfully ran: %s' % ' '.join(busted_cmd))
				except Exception as e:
					logObject.error('Had an issue running BUSTED: %s' % ' '.join(busted_cmd))
					sys.stderr.write('Had an issue running BUSTED: %s\n' % ' '.join(busted_cmd))
					logObject.error(e)
					sys.stderr.write(traceback.format_exc())
					sys.exit(1)
		else:
			gard_cmd = ['hyphy', 'CPU=1', 'gard', '--mode', gard_mode, '--alignment', hg_codo_algn_file,
							  '--output', gard_output, '--output-lf', best_gard_output]

			try:
				subprocess.call(' '.join(gard_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
								executable='/bin/bash', timeout=gard_timeout)
				assert(os.path.isfile(best_gard_output))
				logObject.info('Successfully ran: %s' % ' '.join(gard_cmd))
			except subprocess.TimeoutExpired:
				logObject.error('Timed out running GARD: %s, defaulting to using original alignment in downstream selection analyses.' % ' '.join(gard_cmd))
				sys.stderr.write('Timed out running GARD: %s, defaulting to using original alignment in downstream selection analyses.\n' % ' '.join(gard_cmd))
				logObject.error(traceback.format_exc())
				sys.stderr.write(traceback.format_exc())
				best_gard_output = hg_codo_algn_file
			except Exception as e:
				logObject.error('Had an issue running GARD: %s, defaulting to using original alignment in downstream selection analyses.' % ' '.join(gard_cmd))
				sys.stderr.write('Had an issue running GARD: %s, defaulting to using original alignment in downstream selection analyses.\n' % ' '.join(gard_cmd))
				sys.stderr.write(traceback.format_exc())
				logObject.error(traceback.format_exc())
				best_gard_output = hg_codo_algn_file

			fubar_cmd = ['hyphy', 'CPU=1', 'fubar', '--alignment', best_gard_output]
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
				sys.stderr.write(traceback.format_exc())
				sys.exit(1)

			if not skip_busted:
				busted_cmd = ['hyphy', 'CPU=1', 'busted', '--alignment', best_gard_output]
				try:
					subprocess.call(' '.join(busted_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
									executable='/bin/bash')
					assert(os.path.isfile(best_gard_output + '.BUSTED.json'))
					os.system('mv %s %s' % (best_gard_output + '.BUSTED.json', busted_outdir))
					logObject.info('Successfully ran: %s' % ' '.join(busted_cmd))
				except Exception as e:
					logObject.error('Had an issue running BUSTED: %s' % ' '.join(busted_cmd))
					sys.stderr.write('Had an issue running BUSTED: %s\n' % ' '.join(busted_cmd))
					logObject.error(e)
					sys.stderr.write(traceback.format_exc())
					sys.exit(1)

	except Exception as e:
		sys.stderr.write('Issues with running HYPHY based analyses for ortholog group %s\n' % hg)
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def runHyphyAnalyses(codo_algn_dir, tree_dir, gard_results_dir, fubar_results_dir, busted_results_dir, logObject, 
					 skip_gard=False, skip_busted=False, gard_mode="Faster", gard_timeout=60, threads=1):
	"""f
	Description:
	This function oversees running of HyPhy based analyses (GARD + FUBAR) for ortholog groups and parses resulting
	statistics from resulting JSON files to include in the consolidated report created at the end of zol.
	*******************************************************************************************************************
	Parameters:
	- codo_algn_dir: The directory with codon alignments for each ortholog group.
	- tree_dir: The directory with gene-trees for each ortholog group (recall - these are made using FastTree2 on
			    trimmed codon alignments.
	- gard_results_dir: The directory where GARD result files should be saved.
	- fubar_results_dir: The directory where FUBAR result files should be saved.
	- busted_results_dir: The directory where BUSTED result files should be saved.
	- logObject: A logging object.
	- skip_gard: Boolean indicating whether user has requested to skip GARD analsyis.
	- skip_busted: Boolean indicating whether user has requested to skip BUSTED analysis.
	- gard_mode: Which mode to run GARD analysis using, can either be "Faster" or "Normal".
	- threads: The number of threads to use.
	*******************************************************************************************************************
	"""
	try:
		hyphy_inputs = []
		for caf in os.listdir(codo_algn_dir):
			if not caf.endswith('.msa.fna'): continue
			hg = caf.split('.msa.fna')[0]
			hg_codo_algn_file = codo_algn_dir + caf
			hg_codo_tree_file = tree_dir + hg + '.tre'
			gard_output = gard_results_dir + hg + '.json'
			best_gard_output = gard_results_dir + hg + '.best'
			hyphy_inputs.append([hg, hg_codo_algn_file, hg_codo_tree_file, gard_output, best_gard_output, fubar_results_dir,
							     busted_results_dir, skip_gard, skip_busted, gard_mode, gard_timeout, logObject])

		msg = "Running HyPhy recombination/selection analyses for %d ortholog groups" % len(hyphy_inputs)
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(individualHyphyRun, hyphy_inputs), total=len(hyphy_inputs)):
			pass
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

		busted_pval = {}
		for f in os.listdir(busted_results_dir):
			if f.endswith('.json'):
				hg = f.split('.msa.fna.BUSTED.json')[0]
				if f.endswith('.best.BUSTED.json'):
					hg = f.split('.best.BUSTED.json')[0]
				busted_json_file = busted_results_dir + f
				try:
					with open(busted_json_file) as objr:
						busted_results = json.load(objr)
					pval = float(busted_results['test results']['p-value'])
					busted_pval[hg] = pval
				except:
					busted_pval[hg] = 'NA'	

		return([gard_partitions, fubar_sel_props, fubar_sel_sites, fubar_deba, busted_pval])
	except Exception as e:
		sys.stderr.write('Issues with running HyPhy GARD, BUSTED, or FUBAR analyses.\n')
		logObject.error('Issues with running HyPhy GARD, BUSTED, or FUBAR analyses.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def determineSeqSimProteinAlignment(inputs):
	"""
	Description:
	This function computes the sequence similarity/identity between proteins in a MSA.
	*******************************************************************************************************************
	Parameters:
	- input: A list which can be expanded to the following items:
		- hg: The ortholog group identifier.
		- protein_alignment_file: The protein multiple sequence alignment file for the ortholog group in FASTA format.
		- outf: The output file where to write pairwise sequence similarities.
		- logObject: A logging object.
	*******************************************************************************************************************
	"""
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

def computeBetaRDgc(prot_algn_dir, evo_dir, logObject, threads=1):
	"""
	Description:
	This function computes the BetaRD-gc statistic for each ortholog group - which is an estimate of how sequence
	similarity varies for the ortholog group in relation to other ortholog groups from the gene cluster.
	Note, Beta-RD gene-cluster statistic here is being computed in a different manner than what we did in lsaBGC, it is
	a different statistic.
	*******************************************************************************************************************
	Parameters:
	- prot_algn_dir: The directory with protein alignments for ortholog groups.
	- evo_dir: The workspace/directory where evolutionary analyses are to be performed under.
	- logObject: A logging object.
	- threads: The number of threads to use.
	*******************************************************************************************************************
	Returns:
	- A list of two dictionaries:
		- hg_med_beta_rd: A dictionary mapping ortholog groups to the median BetaRD-gc statistic.
		- hg_max_beta_rd: A dictionary mapping ortholog groups to the max BetaRD-gc statistic.
	*******************************************************************************************************************
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

		msg = "Determining Beta-RDgc statistic for %d ortholog groups" % len(inputs)
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(determineSeqSimProteinAlignment, inputs), total=len(inputs)):
			pass
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
		sys.stderr.write('Issues with calculating Beta-RD gene-cluster for ortholog groups.\n')
		logObject.error('Issues with calculating Beta-RD gene-cluster for ortholog groups.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return([hg_med_beta_rd, hg_max_beta_rd])

def calculateMSAEntropy(inputs):
	"""
	Description:
	This function computes the average entropy statistic for a MSA of proteins for an individual ortholog group.
	*******************************************************************************************************************
	Parameters:
	- input: A list which can be expanded to the following items:
		- hg: The ortholog group identifier.
		- nucl_algn_fasta: The nucleotide/codon multiple sequence alignment file for the ortholog group in FASTA format.
		- outf: The output file where to write the average entropy calculated.
		- logObject: A logging object.
	*******************************************************************************************************************
	"""
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
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def runEntropyAnalysis(codo_algn_trim_dir, upst_algn_dir, evo_dir, logObject, threads=1):
	"""
	Description:
	This function computes the average entropy statistic for ortholog groups.
	*******************************************************************************************************************
	Parameters:
	- codon_algn_trim_dir: The directory with trimmed codon alignments for ortholog groups.
	- upst_algn_dir: The directory with upstream nucleotide alignments for ortholog groups.
	- evo_dir: The workspace/directory where evolutionary analyses are to be performed under.
	- logObject: A logging object.
	- threads: The number of threads to use.
	*******************************************************************************************************************
	Returns:
	- A list of two dictionaries:
		- hg_entropy: A dictionary mapping ortholog groups to the average entropy statistic.
		- hg_upst_entropy: A dictionary mapping ortholog groups to the upstream sequence average entropy statistic.
	*******************************************************************************************************************
	"""
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

		msg = "Computing sequence and upstream sequence entropy for ortholog groups"
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(calculateMSAEntropy, inputs), total=len(inputs)):
			pass
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
		sys.stderr.write('Issues with calculating entropy for ortholog groups or their upstream regions.\n')
		logObject.error('Issues with calculating entropy for ortholog groups or their upstream regions.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def calculateAmbiguity(codo_algn_dir, codo_algn_trim_dir, logObject):
	"""
	Description:
	This function computes the proportion of ambiguous sites (>10% gaps) in full and trimmed codon alignments.
	*******************************************************************************************************************
	Parameters:
	- codon_algn_dir: The directory with codon alignments for ortholog groups.
	- codon_algn_trim_dir: The directory with trimmed codon alignments for ortholog groups.
	- logObject: A logging object.
	*******************************************************************************************************************
	Returns:
	- A list of two dictionaries:
		- full_amb_prop: A dictionary mapping ortholog groups to the proportion of ambiguous sites in the full codon
		                 alignment.
		- trim_amb_prop: A dictionary mapping ortholog groups to the proportion of ambiguous sites in the trimmed codon
		                 alignment.
	*******************************************************************************************************************
	"""
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
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return([full_amb_prop, trim_amb_prop])

def runTajimasDAnalysisPerHG(inputs):
	"""
	Description:
	This function oversees the calculation of Tajima's D and proportion of segregating sites per ortholog group.
	*******************************************************************************************************************
	Parameters:
	- input: A list which can be expanded to the following items:
		- hg: The ortholog group identifier.
		- trim_codon_align: The trimmed codon multiple sequence alignment file for the ortholog group in FASTA format.
		- outf: The output file where to write the Tajima's D statistic and the proportion of segregating sites for the
		        ortholog group.
		- logObject: A logging object.
	*******************************************************************************************************************
	"""
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
		sys.stderr.write('Issues with calculating Tajima\'s D for ortholog group %s.\n' % hg)
		logObject.error('Issues with calculating Tajima\'s D for ortholog group %s.' % hg)
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def determineBGCAndViralScores(pfam_annotations, logObject):
	"""
	Description:
	Maps Pfam annotations per homolog group to BGC and viral scores from GECCO weights and v-Scores, respecitvely.
	********************************************************************************************************************
	Parameters:
	- pfam_annotations: dictionary of pfam annotations for homolog groups.
	- logObject: a Python logging object.
	********************************************************************************************************************
	Returns:
	- gecco_weights: 
	- vscores: 
	********************************************************************************************************************
	"""

	gecco_weights = defaultdict(lambda: 'NA')
	vscores = defaultdict(lambda: 'NA')
	zol_data_directory = str(os.getenv("ZOL_DATA_PATH")).strip()
	gecco_weights_file = None
	vscore_file = None
	if zol_data_directory != 'None':
		try:
			zol_data_directory = os.path.abspath(zol_data_directory) + '/'
			gecco_weights_file = zol_data_directory + 'GECCO_Weights.txt'
			vscore_file = zol_data_directory + 'VScoreDataNormalized.csv'
		except:
			pass
	if gecco_weights_file == None or not os.path.isfile(gecco_weights_file) or vscore_file == None or not os.path.isfile(vscore_file):
		sys.stderr.write('Warning: weight/score files do not appear to be setup or setup properly!\n')
		return([gecco_weights, vscores])
	else:
		pf_vscores = {}
		pf_gecco_weights = {}
		try:
			with open(vscore_file) as ovsf:
				csv_handle = list(csv.reader(ovsf))
				for i, ls in enumerate(csv_handle):
					if i == 0: continue
					if not ls[0].startswith('PF'): continue
					pf = ls[0].strip().split('.')[0]
					pf_vscores[pf] = float(ls[2])

			with open(gecco_weights_file) as ogwf:
				for line in ogwf:
					line = line.strip() 
					pf, weight = line.split('\t')
					pf = pf.strip().split('.')[0]
					pf_gecco_weights[pf] = float(weight)
		
			for hg in pfam_annotations:
				max_vscore = -1.0
				max_bgc = -7.0
				for pf in pfam_annotations[hg][1]:
					if not pf.startswith('PF'): continue
					pf = pf.split('.')[0]

					if pf in pf_vscores:
						if pf_vscores[pf] > max_vscore:
							max_vscore = pf_vscores[pf]
					if pf in pf_gecco_weights:
						if pf_gecco_weights[pf] > max_bgc:
							max_bgc = pf_gecco_weights[pf]
								
				if max_bgc == -7:
					max_bgc = 'NA'
				if max_vscore == -1.0:
					max_vscore = 'NA'
				
				gecco_weights[hg] = max_bgc
				vscores[hg] = max_vscore
			return([gecco_weights, vscores])
		
		except Exception as e:
			sys.stderr.write('Issues with mapping pfam domain annotations to BGC and viral scores from GECCO and vScore.\n')
			logObject.error('Issues with calculating Tajima\'s D for ortholog groups.')
			sys.stderr.write(str(e) + '\n')
			sys.stderr.write(traceback.format_exc())
			sys.exit(1)

def runTajimasDAnalysis(codo_algn_trim_dir, evo_dir, logObject, threads=1):
	"""
	Description:
	This function runs Tajima's D analysis for all ortholog groups.
	********************************************************************************************************************
	Parameters:
	- codo_algn_trim_dir: The directory where trimmed codon alignments are stored for ortholog groups.
	- evo_dir: The workspace/directory where Tajima's D analyses should be performed under.
	- logObject: A logging object.
	- threads: The number of threads to use.
	********************************************************************************************************************
	Returns:
	- A list with two items:
		- hg_tajimas_d: A dictionary mapping ortholog groups to their Tajima's D estimates.
		- hg_seg_sites: A dictionary mapping ortholog groups to the proportion of segregating sites in their trimmed
		                codon alignments.
	********************************************************************************************************************
	"""
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

		msg = "Computing Tajima's D statistic (using trimmed codon alignments) for %d ortholog groups" % len(inputs)
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(runTajimasDAnalysisPerHG, inputs), total=len(inputs)):
			pass
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
		sys.stderr.write('Issues with calculating Tajima\'s D for ortholog groups.\n')
		logObject.error('Issues with calculating Tajima\'s D for ortholog groups.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def calculateTajimasD(sequences):
	"""
	Description:
	Takes a list of aligned sequences (trimmed codon alignments) and computes Tajima's D statistic. The code for this
	functionality was largely taken from Tom Whalley's Tajima's D implementation in Python and further modified /
	corrected to better match the calculation of the statistic as described by Tajima 1989.
	********************************************************************************************************************
	Parameters:
	- sequences: a list of aligned nucleotide sequences.
	********************************************************************************************************************
	Returns:
	- An estimate of Tajima's D for the input.
	********************************************************************************************************************
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


def compareFocalAndComparatorGeneClusters(focal_genbank_ids, comparator_genbank_ids, codo_algn_trim_dir, upst_algn_dir,
										  logObject, representative_associated_members=None,
										  impute_broad_conservation=False):
	"""
	Description:
	This function performs comparative analyses between focal and comparator/complementary gene clusters if requested by
	the user. Will compute conservation percentages between the two sets and also FST for each ortholog group for the
	focal gene cluster and upstream regions.
	*******************************************************************************************************************
	Parameters:
	- focal_genbank_ids: A set of gene cluster identifiers which correspond to the focal set delineated by the user.
	- comparator_genbank_ids: A set of gene cluster identifiers which correspond to the comparator or complementary set
	                          either specified by the user or automatically determined.
	- codo_algn_trim_dir: The directory where the trimmed codon alignments for ortholog groups are stored.
	- upst_algn_dir: The directory where the alignments of the upstream sequences for ortholog groups are stored.
	- logObject: A logging object.
	- representative_associated_members: A mapping of gene clusters (including those removed due to redundancy from
	                                     dereplication) to representative gene clusters.
	- impute_broad_conservation: Whether to impute conservation comprehensively, regarding members of a set of
	                             similar gene clusters as featuring an ortholog group if their respetive representative
	                             gene cluster had the ortholog group.
	*******************************************************************************************************************
	Returns:
	- comp_stats: A dictionary which contains four dictionaries:
		- prop_foc_with: A dictionary which maps ortholog groups to the proportion of focal gene clusters with them.
		- prop_com_with: A dictionary which maps ortholog groups to the proportion of comparator gene clusters with them.
		- fst: A dictionary which maps ortholog groups to the FST for the focal gene clusters.
		- upst_fst:
	comp_stats[hg] = {'prop_foc_with': prop_foc_with, 'prop_com_with': prop_com_with, 'fst': fst, 'fst_upst': upst_fst}
	*******************************************************************************************************************
	"""
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

def createOrthoGroupMatrixFromPrecomputedFile(precomputed_orthogroups_file, fo_prot_dir, ortho_matrix_file, logObject):
	"""
	Description:
	Create an orthogroup matrix from a file listing ortholog group designations for locus tags.
	********************************************************************************************************************
	Parameters:
	- precomputed_orthogroups_file: The file listing the precomputed designations of orthogroups per locus tag.
	- fo_prot_dir: The directory of proteins extracted from CDS features in gene cluster GenBank files.
	- ortho_matrix_file: The output file to write the final orthogroup vs. sample matrix used by zol in subsequent steps.
	- logObject: A logging object.
	********************************************************************************************************************
	"""
	try:
		lt_to_og = {}
		try:
			with open(precomputed_orthogroups_file) as opof:
				for line in opof:
					line = line.strip()
					ls = line.split('\t')
					lt_to_og[ls[0]] = ls[1]
		except:
			msg = 'Issues processing precomputed orthogroups designations in the file %s' % precomputed_orthogroups_file
			sys.stderr.write(msg + '\n')
			logObject.error(msg)
			raise RuntimeError()

		all_gcs = set([])
		all_ogs = set([])
		og_gc_lts = defaultdict(lambda: defaultdict(set))
		for prot_faa in os.listdir(fo_prot_dir):
			prot_faa_file = fo_prot_dir + prot_faa
			gc = '.faa'.join(prot_faa.split('.faa')[:-1])
			all_gcs.add(gc)
			with open(prot_faa_file) as opff:
				for rec in SeqIO.parse(opff, 'fasta'):
					assert(rec.id.startswith(gc + '|'))
					lt = rec.id.split(gc + '|')[1]
					og = None
					try:
						og = lt_to_og[lt]
					except:
						msg = 'Issues finding corresponding orthogroup designation for the locus tag' % lt
						sys.stderr.write(msg + '\n')
						logObject.error(msg)
						raise RuntimeError()
					all_ogs.add(og)
					og_gc_lts[og][gc].add(rec.id)

		outf = open(ortho_matrix_file, 'w')
		outf.write('\t'.join(['Sample'] + sorted(all_gcs)) + '\n')
		for og in sorted(all_ogs):
			printlist = [og]
			for gc in sorted(list(all_gcs)):
				printlist.append(', '.join(og_gc_lts[og][gc]))
			outf.write('\t'.join(printlist) + '\n')
		outf.close()

	except Exception as e:
		sys.stderr.write('Issues with creating a sample vs. ortholog group matrix file from pre-computed locus tag to orthogroup designations.\n')
		logObject.error('Issues with creating a sample vs. ortholog group matrix file from pre-computed locus tag to orthogroup designations')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def consolidateReport(consensus_prot_seqs_faa, comp_stats, hg_stats, annotations, evo_stats, type_weights, final_report_xlsx,
					  final_report_tsv, full_prot_clusters, logObject, domain_mode=False, run_hyphy=False, ces=False):
	"""
	Description:
	This function creates the final consolidated TSV and XLSX report for zol where each row corresponds to an ortholog
	group.
	********************************************************************************************************************
	Parameters:
	- consensus_prot_seqs_faa: FASTA file with consensus protein sequences for each ortholog group.
	- comp_stats: A dictionary with information from comparative analysis of gene-clusters in focal vs. rest/comparator
	              set of gene-cluster for each ortholog group.
	- hg_stats: A dictionary with general information/statistics for each ortholog group.
	- annotations: A dictionary containing annotation information for each ortholog group.
	- evo_stats: A dictionary containing evolutionary statistics for each ortholog group.
	- type_weights: A dictionary containing orthogroup specific weights for BGC or viral related.
	- final_report_xlsx: The path to the final XLSX report spreadsheet.
	- final_report_tsv: The path to the final TSV report table.
	- logObject: A logging object.
	- run_hyphy: Whether HyPhy analysis was requested by the user.
	- ces: Whether comprehensive reporting of stats for all ortholog groups, regardless of whether they are found in
	       multiple copies in some gene clusters (some statistics are filtered by default for such ortholog groups).
	********************************************************************************************************************
	"""

	try:
		# Note to self, eventually conditionally display all columns (e.g. FUBAR columns) when requested by user
		# to avoid having columns with NA values.
		header = ['Ortholog Group (OG) ID', 'OG is Single Copy?', 'Proportion of Total Gene Clusters with OG',
				  'OG Median Length (bp)', 'OG Consensus Order', 'OG Consensus Direction']
		if domain_mode != False:
			header = ['Ortholog Group (OG) ID', 'OG is Single Copy?', 'Proportion of Total Gene Clusters with OG',
				      'OG Median Length (bp)', 'Single-Linkage Full Protein Cluster', 'OG Consensus Order', 
				      'OG Consensus Direction']
		if comp_stats != None:
			header += ['Proportion of Focal Gene Clusters with OG', 'Proportion of Comparator Gene Clusters with OG',
					   'Fixation Index', 'Upstream Region Fixation Index']
		header += ['Tajima\'s D', 'Proportion of Filtered Codon Alignment is Segregating Sites', 'Entropy',
				   'Upstream Region Entropy', 'Median Beta-RD-gc', 'Max Beta-RD-gc',
				   'Proportion of sites which are highly ambiguous in codon alignment',
				   'Proportion of sites which are highly ambiguous in trimmed codon alignment', 'Median GC',
				   'Median GC Skew', 'BGC score (GECCO weights)', 'Viral score (V-Score)']
		if run_hyphy:
			header += ['GARD Partitions Based on Recombination Breakpoints',
			           'Number of Sites Identified as Under Positive or Negative Selection by FUBAR',
				       'Average delta(Beta, Alpha) by FUBAR across sites',
				       'Proportion of Sites Under Selection which are Positive',
					   "P-value for gene-wide episodic selection by BUSTED"]
		header += ['Custom Annotation (E-value)', 'KO Annotation (E-value)', 'PGAP Annotation (E-value)',
				   'PaperBLAST Annotation (E-value)', 'CARD Annotation (E-value)', 'IS Finder (E-value)',
				   'MIBiG Annotation (E-value)', 'VOG Annotation (E-value)', 'VFDB Annotation (E-value)',
				   'Pfam Domains', 'CDS Locus Tags', 'OG Consensus Sequence']

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
			hg_scs = util.gatherValueFromDictForHomologGroup(hg, hg_stats['hg_single_copy_status'])
			hg_cons = util.gatherValueFromDictForHomologGroup(hg, hg_stats['hg_prop_samples'])
			hg_mlen = util.gatherValueFromDictForHomologGroup(hg, hg_stats['hg_median_lengths'])
			hg_lts = 'NA'
			if hg in hg_stats['hg_locus_tags']: hg_lts = '; '.join(sorted(hg_stats['hg_locus_tags'][hg]))
			hg_full_amb = util.gatherValueFromDictForHomologGroup(hg, hg_stats['hg_full_ambiguity'])
			hg_trim_amb = util.gatherValueFromDictForHomologGroup(hg, hg_stats['hg_trim_ambiguity'])
			hg_gc = util.gatherValueFromDictForHomologGroup(hg, hg_stats['hg_median_gc'])
			hg_gcs = util.gatherValueFromDictForHomologGroup(hg, hg_stats['hg_median_gcskew'])
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
				hg_bpva = util.gatherValueFromDictForHomologGroup(hg, evo_stats['busted_pvals'])
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
			if domain_mode:
				fp_clusters = set([])
				for dog in hg_stats['hg_locus_tags'][hg]:
					fp = '|'.join(dog.split('|')[:-2])
					fpc = full_prot_clusters[fp]
					fp_clusters.add(fpc)
				try:
					assert(len(fp_clusters) == 1)
				except:
					msg = 'Issue with validating all domain instances belong to one full protein cluster - shouldn\'t happen - please report on Github issues!'
					sys.stderr.write(msg)
					logObject.error(msg + '\n')

				fp_cluster = list(fp_clusters)[0]
				row = ['D' + hg, hg_scs, hg_cons, hg_mlen, fp_cluster, hg_ordr, hg_dire]
			if comp_stats != None:
				fp = comp_stats[hg]['prop_foc_with']
				cp = comp_stats[hg]['prop_com_with']
				if hg_scs == True or ces:
					fst = comp_stats[hg]['fst']
					fst_upst = comp_stats[hg]['fst_upst']
				row += [fp, cp, fst, fst_upst]
			hg_gw = type_weights['bgc'][hg]
			hg_vs = type_weights['viral'][hg]
			row += [hg_tajd, hg_segs, hg_entr, hg_upst_entr, hg_med_brdgc, hg_max_brdgc, hg_full_amb, hg_trim_amb,
					hg_gc, hg_gcs, hg_gw, hg_vs]
			if run_hyphy:
				row += [hg_gpar, hg_ssit, hg_deba, hg_spro, hg_bpva]

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
		dd_sheet.write(0, 0, 'Data Dictionary describing columns of "ZoL Results" spreadsheet can be found below and on zol\'s Wiki page at:')
		dd_sheet.write(1, 0, 'https://github.com/Kalan-Lab/zol/wiki/3.-more-info-on-zol#explanation-of-report')


		wrap_format = workbook.add_format({'text_wrap': True, 'valign': 'vcenter', 'align': 'center', 'border': 1})
		header_format = workbook.add_format({'bold': True, 'text_wrap': True, 'valign': 'top', 'fg_color': '#FFFFFF', 'border': 1})

		data_dict_zol = data_dictionary.zol_dd()
		data_dict_zol_df = util.loadTableInPandaDataFrameFromString(data_dict_zol)
		worksheet_dd =  writer.sheets['Data Dictionary']
		worksheet_dd.set_column(1, 3, 50)

		for col_num, value in enumerate(data_dict_zol_df.columns.values):
			worksheet_dd.write(3, col_num + 1, value, header_format)

		colnames = ['Column', 'Description', 'Notes']
		for index, row in data_dict_zol_df.iterrows():
			row_ind = index + 4
			format = wrap_format
			for col_ind in range(0,3):
				col_name = colnames[col_ind]
				worksheet_dd.write(row_ind, col_ind+1, row[col_name], format)

		numeric_columns = {'Proportion of Total Gene Clusters with OG', 'Proportion of Focal Gene Clusters with OG',
						   'Proportion of Comparator Gene Clusters with OG', 'Fixation Index',
						   'Upstream Region Fixation Index', 'OG Median Length (bp)', 'OG Consensus Order',
						   'Tajima\'s D', 'Entropy', 'Upstream Region Entropy',
						   'GARD Partitions Based on Recombination Breakpoints',
						   'Number of Sites Identified as Under Positive or Negative Selection by FUBAR',
						   'Proportion of Sites Under Selection which are Positive', 
						   'P-value for gene-wide episodic selection by BUSTED', 'Median Beta-RD-gc',
						   'Max Beta-RD-gc', 'Proportion of Filtered Codon Alignment is Segregating Sites',
						   'Proportion of sites which are highly ambiguous in codon alignment',
						   'Proportion of sites which are highly ambiguous in trimmed codon alignment',
						   'Average delta(Beta, Alpha) by FUBAR across sites', 'Median GC', 'Median GC Skew', 
						   'BGC score (GECCO weights)', 'Viral score (V-Score)'}

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
			hfo = ['G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V']
			if domain_mode:
				hfo = ['H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W']
			# prop focal gene-clusters with hg
			worksheet.conditional_format(hfo[0] + '2:' + hfo[0] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#f7de99", 'max_color': "#c29006",'min_type': 'num', 'max_type': 'num',
										  "min_value": 0.0, "max_value": 1.0})
			# prop comparator gene-clusters with hg
			worksheet.conditional_format(hfo[1] + '2:' + hfo[1] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#f7de99", 'max_color': "#c29006",'min_type': 'num', 'max_type': 'num',
										  "min_value": 0.0, "max_value": 1.0})
			# fst
			worksheet.conditional_format(hfo[2] + '2:' + hfo[2] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#aecaf5", 'max_color': "#6198ed",'min_type': 'num', 'max_type': 'num',
										  "min_value": 0.0, "max_value": 1.0})

			# upstream region fst
			worksheet.conditional_format(hfo[3] + '2:' + hfo[3] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#aecaf5", 'max_color': "#6198ed",'min_type': 'num', 'max_type': 'num',
										  "min_value": 0.0, "max_value": 1.0})

			# taj-d
			worksheet.conditional_format(hfo[4] + '2:' + hfo[4] + str(num_rows),
										 {'type': '3_color_scale', 'min_color': "#f7a09c", "mid_color": "#e0e0e0", 'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										  'max_color': "#87cefa", "min_value": -2.0, "mid_value": 0.0, "max_value": 2.0})
			# prop seg sites
			worksheet.conditional_format(hfo[5] + '2:' + hfo[5] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#eab3f2", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#a37ba8", "min_value": 0.0, "max_value": 1.0})

			# entropy
			worksheet.conditional_format(hfo[6] + '2:' + hfo[6]+ str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#f7a8bc", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#fa6188", "min_value": 0.0, "max_value": 1.0})

			# upstream region entropy
			worksheet.conditional_format(hfo[7] + '2:' + hfo[7] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#f7a8bc", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#fa6188", "min_value": 0.0, "max_value": 1.0})

			# beta-rd gc
			worksheet.conditional_format(hfo[8] + '2:' + hfo[8] + str(num_rows),
										 {'type': '3_color_scale', 'min_color': "#fac087", "mid_color": "#e0e0e0", 'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										 'max_color': "#9eb888", "min_value": 0.75, "mid_value": 1.0, "max_value": 1.25})
			# max beta-rd gc
			worksheet.conditional_format(hfo[9] + '2:' + hfo[9] + str(num_rows),
										 {'type': '3_color_scale', 'min_color': "#fac087", "mid_color": "#e0e0e0", 'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										  'max_color': "#9eb888", "min_value": 0.75, "mid_value": 1.0, "max_value": 1.25})

			# ambiguity full ca
			worksheet.conditional_format(hfo[10] + '2:' + hfo[10] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#ed8c8c", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#ab1616", "min_value": 0.0, "max_value": 1.0})

			# ambiguity trim ca
			worksheet.conditional_format(hfo[11] + '2:' + hfo[11] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#ed8c8c", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#ab1616", "min_value": 0.0, "max_value": 1.0})

			# GC
			worksheet.conditional_format(hfo[12] + '2:' + hfo[12] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#abffb7", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#43bf55", "min_value": 0.0, "max_value": 1.0})

			# GC Skew
			worksheet.conditional_format(hfo[13] + '2:' + hfo[13] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#c7afb4", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#965663", "min_value": -2.0, "max_value": 2.0})
			
			# BGC score
			worksheet.conditional_format(hfo[14] + '2:' + hfo[14] + str(num_rows),
									{'type': '2_color_scale', 'min_color': "#f5aca4", 'min_type': 'num', 'max_type': 'num',
									'max_color': "#c75246", "min_value": -7.0, "max_value": 13.0})

			# viral score
			worksheet.conditional_format(hfo[15] + '2:' + hfo[15] + str(num_rows),
									{'type': '2_color_scale', 'min_color': "#dfccff", 'min_type': 'num', 'max_type': 'num',
									'max_color': "#715a96", "min_value": 0.0, "max_value": 10.0})
			
		else:
			hfo = ['G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R']
			if domain_mode:
				hfo = ['H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S']
			# taj-d
			worksheet.conditional_format(hfo[0] + '2:' + hfo[0] + str(num_rows),
										 {'type': '3_color_scale', 'min_color': "#f7a09c", "mid_color": "#e0e0e0",'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										  'max_color': "#87cefa", "min_value": -2.0, "mid_value": 0.0, "max_value": 2.0})

			# prop seg sites
			worksheet.conditional_format(hfo[1] + '2:' + hfo[1] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#eab3f2", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#a37ba8", "min_value": 0.0, "max_value": 1.0})

			# entropy
			worksheet.conditional_format(hfo[2] + '2:' + hfo[2] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#f7a8bc", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#fa6188", "min_value": 0.0, "max_value": 1.0})

			# upstream region entropy
			worksheet.conditional_format(hfo[3] + '2:' + hfo[3] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#f7a8bc", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#fa6188", "min_value": 0.0, "max_value": 1.0})

			# median beta-rd gc
			worksheet.conditional_format(hfo[4] + '2:' + hfo[4] + str(num_rows),
										 {'type': '3_color_scale', 'min_color': "#fac087", "mid_color": "#e0e0e0",'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										  'max_color': "#9eb888", "min_value": 0.75, "mid_value": 1.0, "max_value": 1.25})
			# max beta-rd gc
			worksheet.conditional_format(hfo[5] + '2:' + hfo[5] + str(num_rows),
										 {'type': '3_color_scale', 'min_color': "#fac087", "mid_color": "#e0e0e0",'min_type': 'num', 'max_type': 'num', 'mid_type': 'num',
										  'max_color': "#9eb888", "min_value": 0.75, "mid_value": 1.0, "max_value": 1.25})

			# ambiguity full ca
			worksheet.conditional_format(hfo[6] + '2:' + hfo[6] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#ed8c8c", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#ab1616", "min_value": 0.0, "max_value": 1.0})

			# ambiguity trim ca
			worksheet.conditional_format(hfo[7] + '2:' + hfo[7] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#ed8c8c", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#ab1616", "min_value": 0.0, "max_value": 1.0})

			# GC
			worksheet.conditional_format(hfo[8] + '2:' + hfo[8] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#abffb7", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#43bf55", "min_value": 0.0, "max_value": 1.0})

			# GC Skew
			worksheet.conditional_format(hfo[9] + '2:' + hfo[9] + str(num_rows),
										 {'type': '2_color_scale', 'min_color': "#c7afb4", 'min_type': 'num', 'max_type': 'num',
										  'max_color': "#965663", "min_value": -2.0, "max_value": 2.0})

			# BGC score
			worksheet.conditional_format(hfo[10] + '2:' + hfo[10] + str(num_rows),
									{'type': '2_color_scale', 'min_color': "#f5aca4", 'min_type': 'num', 'max_type': 'num',
									'max_color': "#c75246", "min_value": -7.0, "max_value": 13.0})

			# viral score
			worksheet.conditional_format(hfo[11] + '2:' + hfo[11] + str(num_rows),
									{'type': '2_color_scale', 'min_color': "#dfccff", 'min_type': 'num', 'max_type': 'num',
									'max_color': "#715a96", "min_value": 0.0, "max_value": 10.0})

		worksheet.autofilter('A1:BA' + str(num_rows))
		worksheet.filter_column(2, 'x >= 0.1')
		workbook.close()

	except Exception as e:
		sys.stderr.write('Issues creating consolidated results files.\n')
		logObject.error('Issues creating consolidated results files.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def plotHeatmap(hg_stats, genbanks, plot_result_pdf, work_dir, logObject, height=7, width=10,
				full_genbank_labels=False):
	"""
	Description:
	This function create a heatmap figure showing the presence of ortholog groups across gene clusters.
	********************************************************************************************************************
	Parameters:
	- hg_stats: A dictionary containing general statistics for each ortholog group (including which gene clusters have
	            them).
	- genbanks: The set of gene-clusters to account for in the heatmap.
	- plot_result_pdf: The resulting heatmap plot PDF path.
	- work_dir: The workspace/directory where to store intermediate files for plotting.
	- logObject: A logging object.
	- height: The height of the plot in inches.
	- width: The width of the plot in inches.
	- full_genbank_labels: Whether to use the full gene cluster identifiers as labels.
	********************************************************************************************************************
	"""
	try:
		representative_genbanks = set([])
		for gbk in genbanks:
			gbk_prefix = None
			if gbk.endswith('.gbk') or gbk.endswith('.gbff') or gbk.endswith('.genbank'):
				gbk_prefix = '.'.join(gbk.split('/')[-1].split('.')[:-1])
			assert (gbk_prefix != None)
			representative_genbanks.add(gbk_prefix)

		# create input tracks
		ml_track_file = work_dir + 'OG_Median_Length_Info.txt'
		hm_track_file = work_dir + 'OG_Heatmap_Info.txt'
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

		rscript_path = work_dir + 'clusterHeatmap.R'
		util.clusterHeatmapR(ml_track_file, hm_track_file, plot_result_pdf, height, width, rscript_path, logObject)
		plot_cmd = ['Rscript', rscript_path]
		try:
			subprocess.call(' '.join(plot_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
							executable='/bin/bash')
			assert (os.path.isfile(plot_result_pdf))
			logObject.info('Successfully ran: %s' % ' '.join(plot_cmd))
		except Exception as e:
			logObject.error('Had an issue running R based plotting - potentially because of R setup issues in conda: %s' % ' '.join(plot_cmd))
			sys.stderr.write('Had an issue running R based plotting - potentially because of R setup issues in conda: %s\n' % ' '.join(plot_cmd))
			logObject.error(e)
			sys.stderr.write(traceback.format_exc())
			sys.exit(1)
	except Exception as e:
		sys.stderr.write('Issues creating heatmap visualization in zol.\n')
		logObject.error('Issues creating heatmap visualization in zol.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
