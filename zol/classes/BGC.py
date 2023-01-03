import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from operator import itemgetter
import traceback
import copy
import _pickle as cPickle

lsaBGC_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-3])
gecco_pickle_weights_file_file = lsaBGC_main_directory + '/db/GECCO_PF_Weights.pkl'

class BGC:
	def __init__(self, bgc_genbank, bgc_id, is_expansion_bgc, prediction_method='ANTISMASH'):
		self.bgc_genbank = bgc_genbank
		self.bgc_id = bgc_id
		self.is_expansion_bgc = is_expansion_bgc
		self.prediction_method = prediction_method.upper()
		self.gene_information = None
		self.cluster_information = None

	def parseGECCO(self, comprehensive_parsing=True, flank_size=2000):
		""" Function to parse BGC Genbank produced by GECCO BGC."""
		domains = []
		domain_weights = {}

		gecco_pfam_weights_pickle_handle = open(gecco_pickle_weights_file_file, "rb")
		gecco_pfam_weights = cPickle.load(gecco_pfam_weights_pickle_handle)
		rec = SeqIO.read(self.bgc_genbank, 'genbank')
		full_sequence = str(rec.seq)
		for feature in rec.features:
			if feature.type == 'misc_feature':
				start = feature.location.start + 1
				end = feature.location.end
				aSDomain = "NA"
				description = "NA"
				dom_weight = -7
				try:
					aSDomain = feature.qualifiers['standard_name'][0]
				except:
					pass
				try:
					description = feature.qualifiers['function'][0]
				except:
					pass
				try:
					dom_weight = gecco_pfam_weights[aSDomain]
				except:
					pass
				domain_weights[aSDomain + '|' + str(start+1) + '|' + str(end)] = dom_weight
				domains.append({'start': start + 1, 'end': end, 'type': feature.type, 'aSDomain': aSDomain, 'description': description, 'is_multi_part': False})

		product = 'NA'
		try:
			product = rec.annotations['structured_comment']['GECCO-Data']['biosyn_class']
		except:
			pass
		bgc_info = [{'prediction_method': self.prediction_method, 'detection_rule': 'NA', 'product': product, 'contig_edge': 'NA', 'full_sequence': full_sequence}]

		# determine top 10% of domains with highest GECCO CRF weights (as recommended by Martin Larralde)
		num_total_domains = len(domain_weights)
		core_domains = set([])
		for i, d in enumerate(sorted(domain_weights.items(), key=itemgetter(1), reverse=True)):
			if i <= num_total_domains*0.1:
				core_domains.add(d[0])

		# sys.stderr.write('Processing %s\n' % self.bgc_genbank)
		genes = {}
		core_genes = set([])
		gene_order = {}

		for feature in rec.features:
			if feature.type == "CDS":
				lt = feature.qualifiers.get('locus_tag')[0]
				start = feature.location.start + 1
				end = feature.location.end
				direction = "-" if feature.location.strand == -1 else "+"

				try:
					product = feature.qualifiers.get('product')[0]
				except:
					product = "hypothetical protein"

				grange = set(range(start, end + 1))

				gene_domains = []
				core_overlap = False
				for d in domains:
					drange = set(range(d['start'], d['end'] + 1))
					if len(drange.intersection(grange)) > 0:
						gene_domains.append(d)
						if (d['aSDomain'] + '|' + str(d['start']) + '|' + str(d['end'])) in core_domains:
							core_overlap = True
							core_genes.add(lt)

				gene_order[lt] = start

				prot_seq, nucl_seq, nucl_seq_with_flanks, relative_start, relative_end = [None] * 5
				if comprehensive_parsing:
					prot_seq = feature.qualifiers.get('translation')[0]

					flank_start = start - flank_size
					flank_end = end + flank_size

					if flank_start < 1: flank_start = 1

					if flank_end >= len(full_sequence): flank_end = None
					if end >= len(full_sequence): end = None

					if end:
						nucl_seq = full_sequence[start - 1:end]
					else:
						nucl_seq = full_sequence[start - 1:]
						end = len(full_sequence)

					if flank_end:
						nucl_seq_with_flanks = full_sequence[flank_start - 1:flank_end]
					else:
						nucl_seq_with_flanks = full_sequence[flank_start - 1:]

					gene_length = end - start

					relative_start = nucl_seq_with_flanks.find(nucl_seq)
					relative_end = relative_start + gene_length

					if direction == '-':
						nucl_seq = str(Seq(nucl_seq).reverse_complement())
						nucl_seq_with_flanks = str(Seq(nucl_seq_with_flanks).reverse_complement())
						relative_start = nucl_seq_with_flanks.find(nucl_seq)
						relative_end = relative_start + gene_length

				genes[lt] = {'bgc_name': self.bgc_id, 'start': start, 'end': end, 'direction': direction,
							 'product': product, 'prot_seq': prot_seq, 'nucl_seq': nucl_seq,
							 'nucl_seq_with_flanks': nucl_seq_with_flanks, 'gene_domains': gene_domains,
							 'core_overlap': core_overlap, 'relative_start': relative_start,
							 'relative_end': relative_end, 'is_expansion_bgc': self.is_expansion_bgc,
							 'is_multi_part': False}

		number_of_core_gene_groups = 0
		tmp = []
		for lt in sorted(gene_order.items(), key=itemgetter(1), reverse=True):
			if lt[0] in core_genes:
				tmp.append(lt[0])
			elif len(tmp) > 0:
				number_of_core_gene_groups += 1
				tmp = []
		if len(tmp) > 0:
			number_of_core_gene_groups += 1

		for i, pc in enumerate(bgc_info):
			bgc_info[i]['count_core_gene_groups'] = number_of_core_gene_groups

		self.gene_information = genes
		self.cluster_information = bgc_info


	def parseDeepBGC(self, comprehensive_parsing=True, flank_size=2000):
		""" Function to parse BGC Genbank produced by DeepBGC. """
		domains = []
		full_sequence = ""
		domain_score = {}
		product = "NA|NA"
		with open(self.bgc_genbank) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				full_sequence = str(rec.seq)
				for feature in rec.features:
					if feature.type == 'PFAM_domain':
						start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
						end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
						aSDomain = "NA"
						description = "NA"
						deepbgc_score = 0.0
						try:
							aSDomain = feature.qualifiers.get('db_xref')[0]
						except:
							pass
						try:
							description = feature.qualifiers.get('description')[0]
						except:
							pass
						try:
							deepbgc_score = float(feature.qualifiers.get('deepbgc_score')[0])
						except:
							pass
						domain_score[aSDomain + '|' + str(start+1) + '|' + str(end)] = deepbgc_score
						domains.append({'start': start + 1, 'end': end, 'type': feature.type, 'aSDomain': aSDomain, 'description': description, 'is_multi_part': False})
					elif feature.type == 'cluster':
						product_class = "NA"
						product_activity = "NA"
						try:
							product_activity = feature.qualifiers.get('product_activity')[0]
						except:
							pass
						try:
							product_classes = []
							product_class_scpus = feature.qualifiers.get('product_class_score')[0]
							for pci in product_class_scpus.split(','):
								pc, pcs = pci.split('=')
								if float(pcs) >= 0.5:
									product_classes.append(pc)
							product_class = '-'.join(sorted(product_classes))
						except:
							pass
						product = 'PC:' + product_class + '_PA:' + product_activity

		bgc_info = [{'prediction_method': self.prediction_method, 'detection_rule': 'NA', 'product': product, 'contig_edge': 'NA', 'full_sequence': full_sequence}]

		# determine top 10% of domains with lowest e-values
		num_total_domains = len(domain_score)
		core_domains = set([])
		for i, d in enumerate(sorted(domain_score.items(), key=itemgetter(1), reverse=True)):
			if i <= num_total_domains*0.1:
				core_domains.add(d[0])

		# sys.stderr.write('Processing %s\n' % self.bgc_genbank)
		genes = {}
		core_genes = set([])
		gene_order = {}
		with open(self.bgc_genbank) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type == "CDS":
						lt = feature.qualifiers.get('locus_tag')[0]
						start = min([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
						end = max([int(x) for x in str(feature.location)[1:].split(']')[0].split(':')])
						direction = str(feature.location).split('(')[1].split(')')[0]

						try:
							product = feature.qualifiers.get('product')[0]
						except:
							product = "hypothetical protein"

						grange = set(range(start, end + 1))

						gene_domains = []
						core_overlap = False
						for d in domains:
							drange = set(range(d['start'], d['end'] + 1))
							if len(drange.intersection(grange)) > 0:
								gene_domains.append(d)
								if (d['aSDomain'] + '|' + str(d['start']) + '|' + str(d['end'])) in core_domains:
									core_overlap = True
									core_genes.add(lt)

						gene_order[lt] = start

						prot_seq, nucl_seq, nucl_seq_with_flanks, relative_start, relative_end = [None] * 5
						if comprehensive_parsing:
							flank_start = start - flank_size
							flank_end = end + flank_size

							if flank_start < 1: flank_start = 1

							if flank_end >= len(full_sequence): flank_end = None
							if end >= len(full_sequence): end = None

							if end:
								nucl_seq = full_sequence[start - 1:end]
							else:
								nucl_seq = full_sequence[start - 1:]
								end = len(full_sequence)

							if flank_end:
								nucl_seq_with_flanks = full_sequence[flank_start - 1:flank_end]
							else:
								nucl_seq_with_flanks = full_sequence[flank_start - 1:]

							gene_length = end - start

							relative_start = nucl_seq_with_flanks.find(nucl_seq)
							relative_end = relative_start + gene_length

							if direction == '-':
								nucl_seq = str(Seq(nucl_seq).reverse_complement())
								nucl_seq_with_flanks = str(Seq(nucl_seq_with_flanks).reverse_complement())
								relative_start = nucl_seq_with_flanks.find(nucl_seq)
								relative_end = relative_start + gene_length

							try:
								prot_seq = feature.qualifiers.get('translation')[0]
							except:
								prot_seq = Seq(nucl_seq).translate()

						genes[lt] = {'bgc_name': self.bgc_id, 'start': start, 'end': end, 'direction': direction,
									 'product': product, 'prot_seq': prot_seq, 'nucl_seq': nucl_seq,
									 'nucl_seq_with_flanks': nucl_seq_with_flanks, 'gene_domains': gene_domains,
									 'core_overlap': core_overlap, 'relative_start': relative_start,
									 'relative_end': relative_end, 'is_expansion_bgc': self.is_expansion_bgc,
									 'is_multi_part': False}

		number_of_core_gene_groups = 0
		tmp = []
		for lt in sorted(gene_order.items(), key=itemgetter(1), reverse=True):
			if lt[0] in core_genes:
				tmp.append(lt[0])
			elif len(tmp) > 0:
				number_of_core_gene_groups += 1
				tmp = []
		if len(tmp) > 0:
			number_of_core_gene_groups += 1

		for i, pc in enumerate(bgc_info):
			bgc_info[i]['count_core_gene_groups'] = number_of_core_gene_groups

		self.gene_information = genes
		self.cluster_information = bgc_info

	def parseAntiSMASH(self, comprehensive_parsing=True, flank_size=2000):
		""" Functoin to parse BGC Genbank produced by antiSMASH """
		bgc_info = []
		domains = []
		core_positions = set([])
		full_sequence = ""
		with open(self.bgc_genbank) as ogbk:
			domain_feature_types = ['PFAM_domain', 'CDS_motif', 'aSDomain']
			for rec in SeqIO.parse(ogbk, 'genbank'):
				full_sequence = str(rec.seq)
				for feature in rec.features:
					if comprehensive_parsing and feature.type in domain_feature_types:
						start = None
						end = None
						is_multi_part = False
						if not 'join' in str(feature.location):
							start = min([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
							end = max([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
						else:
							is_multi_part = True
							all_starts = []
							all_ends = []
							for exon_coord in str(feature.location)[5:-1].split(', '):
								start = min(
									[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
								end = max([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
								all_starts.append(start);
								all_ends.append(end);
							start = min(all_starts)
							end = max(all_ends)

						aSDomain = "NA"
						description = "NA"
						try:
							aSDomain = feature.qualifiers.get('aSDomain')[0]
						except:
							pass
						try:
							description = feature.qualifiers.get('description')[0]
						except:
							pass
						domains.append({'start': start, 'end': end, 'type': feature.type, 'aSDomain': aSDomain,
										'description': description, 'is_multi_part': is_multi_part})
					elif feature.type == 'protocluster':
						detection_rule = feature.qualifiers.get('detection_rule')[0]
						try:
							product = feature.qualifiers.get('product')[0]
						except:
							product = "NA"
						contig_edge = feature.qualifiers.get('contig_edge')[0]
						bgc_info.append(
							{'detection_rule': detection_rule, 'product': product, 'contig_edge': contig_edge,
							 'full_sequence': str(rec.seq)})
					elif feature.type == 'proto_core':
						if not 'join' in str(feature.location):
							core_start = min([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
							core_end = max([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
							core_positions = core_positions.union(set(range(core_start + 1, core_end + 1)))
						else:
							core_starts = []
							core_ends = []
							for exon_coord in str(feature.location)[5:-1].split(', '):
								start = min([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
								end = max([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
								core_starts.append(start); core_ends.append(end)
							core_positions = core_positions.union(set(range(min(core_starts)+1, max(core_ends)+1)))

		if len(bgc_info) == 0:
			bgc_info = [
				{'detection_rule': 'NA', 'product': 'NA', 'contig_edge': 'NA', 'full_sequence': full_sequence}]

		# sys.stderr.write('Processing %s\n' % self.bgc_genbank)
		genes = {}
		core_genes = set([])
		gene_order = {}
		with open(self.bgc_genbank) as ogbk:
			for rec in SeqIO.parse(ogbk, 'genbank'):
				for feature in rec.features:
					if feature.type == "CDS":
						lt = feature.qualifiers.get('locus_tag')[0]
						start = None
						end = None
						direction = None
						all_coords = []
						is_multi_part = False
						if not 'join' in str(feature.location):
							start = min([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
							end = max([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
							direction = str(feature.location).split('(')[1].split(')')[0]
							all_coords.append([start, end, direction])
						else:
							is_multi_part = True
							all_starts = []
							all_ends = []
							all_directions = []
							for exon_coord in str(feature.location)[5:-1].split(', '):
								start = min([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
								end = max([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
								direction = exon_coord.split('(')[1].split(')')[0]
								all_starts.append(start); all_ends.append(end); all_directions.append(direction)
								all_coords.append([start, end, direction])
							assert(len(set(all_directions)) == 1)
							start = min(all_starts)
							end = max(all_ends)
							direction = all_directions[0]
						try:
							product = feature.qualifiers.get('product')[0]
						except:
							product = "hypothetical protein"

						rule_based_bgc_cds = False
						try:
							if 'rule-based-clusters' in feature.qualifiers.get('gene_functions')[0]:
								rule_based_bgc_cds = True
						except:
							pass

						grange = set(range(start, end + 1))
						core_overlap = False
						if len(grange.intersection(core_positions)) > 0 and rule_based_bgc_cds:
							core_overlap = True
							core_genes.add(lt)

						gene_order[lt] = start

						prot_seq, nucl_seq, nucl_seq_with_flanks, relative_start, relative_end, gene_domains = [None] * 6
						if comprehensive_parsing:
							prot_seq = feature.qualifiers.get('translation')[0]
							gene_domains = []
							for d in domains:
								drange = set(range(d['start'], d['end'] + 1))
								if len(drange.intersection(grange)) > 0:
									gene_domains.append(d)

							flank_start = start - flank_size
							flank_end = end + flank_size

							if flank_start < 1: flank_start = 1

							if flank_end >= len(full_sequence): flank_end = None
							if end >= len(full_sequence): end = len(full_sequence)

							nucl_seq = ''
							for sc, ec, dc in sorted(all_coords, key=itemgetter(0), reverse=False):
								if ec >= len(full_sequence):
									nucl_seq += full_sequence[sc - 1:]
								else:
									nucl_seq += full_sequence[sc - 1:ec]

							if flank_end:
								nucl_seq_with_flanks = full_sequence[flank_start - 1:flank_end]
							else:
								nucl_seq_with_flanks = full_sequence[flank_start - 1:]

							gene_length = end - start

							relative_start = nucl_seq_with_flanks.find(nucl_seq)
							relative_end = relative_start + gene_length

							if direction == '-':
								nucl_seq = str(Seq(nucl_seq).reverse_complement())
								nucl_seq_with_flanks = str(Seq(nucl_seq_with_flanks).reverse_complement())
								relative_start = nucl_seq_with_flanks.find(nucl_seq)
								relative_end = relative_start + gene_length

						genes[lt] = {'bgc_name': self.bgc_id, 'start': start, 'end': end, 'direction': direction,
									 'product': product, 'prot_seq': prot_seq, 'nucl_seq': nucl_seq,
									 'nucl_seq_with_flanks': nucl_seq_with_flanks, 'gene_domains': gene_domains,
									 'core_overlap': core_overlap, 'relative_start': relative_start,
									 'relative_end': relative_end, 'is_multi_part': is_multi_part,
									 'is_expansion_bgc': self.is_expansion_bgc}

		number_of_core_gene_groups = 0
		tmp = []
		for lt in sorted(gene_order.items(), key=itemgetter(1), reverse=True):
			if lt[0] in core_genes:
				tmp.append(lt[0])
			elif len(tmp) > 0:
				number_of_core_gene_groups += 1
				tmp = []
		if len(tmp) > 0:
			number_of_core_gene_groups += 1

		for i, pc in enumerate(bgc_info):
			bgc_info[i]['count_core_gene_groups'] = number_of_core_gene_groups

		self.gene_information = genes
		self.cluster_information = bgc_info

	def parseGenbanks(self, comprehensive_parsing=True, flank_size=2000):
		"""
		Function to determine whether AntiSMASH, DeepBGC, or GECCO Genbank processing is appropriate.
		If comprehensive parsing is disabled, only minimal info from the BGC will be extracted into the BGC object.
		Gene flanks are not used currently in the software, so might not work as intended, and are an artifact of usage
		in earlier versions of zol-DiscoVary.py
		"""
		if self.prediction_method.upper() == 'ANTISMASH':
			self.parseAntiSMASH(comprehensive_parsing=comprehensive_parsing, flank_size=flank_size)
		elif self.prediction_method.upper() == 'DEEPBGC':
			self.parseDeepBGC(comprehensive_parsing=comprehensive_parsing, flank_size=flank_size)
		elif self.prediction_method.upper() == 'GECCO':
			self.parseGECCO(comprehensive_parsing=comprehensive_parsing, flank_size=flank_size)
		else:
			raise RuntimeError("Unable to parse file because BGC prediction method is not an accepted option!")

	def refineGenbank(self, refined_genbank_file, first_bg, second_bg):
		"""
		Function to prune and update coordinates of BGC Genbank - main functoin of zol-Refiner
		"""
		try:
			rgf_handle = open(refined_genbank_file, 'w')
			start_coord = min([self.gene_information[first_bg]['start'], self.gene_information[first_bg]['end'], self.gene_information[second_bg]['start'], self.gene_information[second_bg]['end']])
			end_coord = max([self.gene_information[first_bg]['start'], self.gene_information[first_bg]['end'], self.gene_information[second_bg]['start'], self.gene_information[second_bg]['end']])
			pruned_coords = set(range(start_coord, end_coord+1))
			with open(self.bgc_genbank) as ogbk:
				recs = [x for x in SeqIO.parse(ogbk, 'genbank')]
				try:
					assert(len(recs) == 1)
				except Exception as e:
					raise RuntimeError(traceback.format_exc())
				for rec in recs:
					original_seq = str(rec.seq)
					filtered_seq = ""
					if end_coord == len(original_seq):
						filtered_seq = original_seq[start_coord-1:]
					else:
						filtered_seq = original_seq[start_coord-1:end_coord]

					new_seq_object = Seq(filtered_seq)

					updated_rec = copy.deepcopy(rec)
					updated_rec.seq = new_seq_object

					updated_features = []
					for feature in rec.features:
						start = None
						end = None
						direction = None
						all_coords = []
						if not 'join' in str(feature.location) and not 'order' in str(feature.location):
							start = min([int(x.strip('>').strip('<')) for x in
										 str(feature.location)[1:].split(']')[0].split(':')]) + 1
							end = max([int(x.strip('>').strip('<')) for x in
									   str(feature.location)[1:].split(']')[0].split(':')])
							direction = str(feature.location).split('(')[1].split(')')[0]
							all_coords.append([start, end, direction])
						elif 'order' in str(feature.location):
							all_starts = []
							all_ends = []
							all_directions = []
							for exon_coord in str(feature.location)[6:-1].split(', '):
								start = min(
									[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
								end = max(
									[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
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
								end = max(
									[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
								direction = exon_coord.split('(')[1].split(')')[0]
								all_starts.append(start)
								all_ends.append(end)
								all_directions.append(direction)
								all_coords.append([start, end, direction])
							start = min(all_starts)
							end = max(all_ends)

						feature_coords = set(range(start, end+1))
						if len(feature_coords.intersection(pruned_coords)) > 0:
							fls = []
							for sc, ec, dc in all_coords:
								updated_start = sc - start_coord + 1
								updated_end = ec - start_coord + 1
								if ec > end_coord:
									if feature.type == 'CDS':
										continue
									else:
										updated_end = end_coord - start_coord + 1  # ; flag1 = True
								if sc < start_coord:
									if feature.type == 'CDS':
										continue
									else:
										updated_start = 1  # ; flag2 = True

								strand = 1
								if dc == '-':
									strand = -1
								fls.append(FeatureLocation(updated_start - 1, updated_end, strand=strand))
							if len(fls) > 0:
								updated_location = fls[0]
								if len(fls) > 1:
									updated_location = sum(fls)
								feature.location = updated_location
								updated_features.append(feature)
					updated_rec.features = updated_features
					SeqIO.write(updated_rec, rgf_handle, 'genbank')
			rgf_handle.close()
		except Exception as e:
      			raise RuntimeError(traceback.format_exc())