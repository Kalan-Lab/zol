#!/usr/bin/env python3

"""
### Program: processNCBIGenBank.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology
"""

# BSD 3-Clause License
#
# Copyright (c) 2023-2025, Kalan-Lab
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
import os
import sys
import argparse
from Bio import SeqIO
import gzip
from zol import util
import traceback
import copy

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: processNCBIGenBank.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Process NCBI Genbanks to create proteome + genbanks with specific locus tag.

	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--input_ncbi_genbank', help='Path to genomic assembly in GenBank format.', required=True)
	parser.add_argument('-g', '--genbank_outdir',
						help='Path to output directory where genbank file should be written. Should already be created!',
						required=True)
	parser.add_argument('-p', '--proteome_outdir',
						help='Path to output directory where proteome file should be written. Should already be created!',
						required=True)
	parser.add_argument('-n', '--name_mapping_outdir',
						help='Path to output directory where gene old-to-new mapping text files should be written. Should already be created!',
						required=True)
	parser.add_argument('-e', '--error_outdir', help='Path to output directory where error files should be written. Should already be created!', required=False, default=None)
	parser.add_argument('-s', '--sample_name', help='Sample name', default='Sample', required=False)
	parser.add_argument('-l', '--locus_tag', help='Locus tag', default="AAAA", required=False)
	parser.add_argument('-r', '--rename_all_lts', action='store_true', help='Only rename locus_tags when needed - e.g. they are missing.', default=False, required=False)
	parser.add_argument('-enl', '--error-no-lt', action='store_true', help='Do not rename locus tags if not found/other issue\n- will result in skipping inclusion of entire genome.', default=False, required=False)
	parser.add_argument('-ent', '--error-no-translation', action='store_true', help='Do not skip CDS without translation\n- will result in skipping inclusion of entire genome.', default=False, required=False)
	args = parser.parse_args()
	return args


def processAndReformatNCBIGenbanks():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE REQUIRED INPUTS
	"""
	myargs = create_parser()

	input_genbank_file = os.path.abspath(myargs.input_ncbi_genbank)
	gbk_outdir = os.path.abspath(myargs.genbank_outdir) + '/'
	pro_outdir = os.path.abspath(myargs.proteome_outdir) + '/'
	map_outdir = os.path.abspath(myargs.name_mapping_outdir) + '/'

	try:
		assert (util.is_genbank(input_genbank_file))
	except Exception as e:
		raise RuntimeError('Issue with input Genbank file from NCBI.')

	outdirs = [gbk_outdir, pro_outdir, map_outdir]
	for outdir in outdirs:
		if not os.path.isdir(outdir):
			sys.stderr.write(f"Output directory {outdir} does not exist! Please create and retry program.")

	"""
	PARSE OPTIONAL INPUTS
	"""

	sample_name = myargs.sample_name
	locus_tag = myargs.locus_tag
	rename_all_flag = myargs.rename_all_lts
	error_no_translation = myargs.error_no_translation
	error_no_lt = myargs.error_no_lt
	error_outdir = myargs.error_outdir

	error_outfile_handle = None
	if error_outdir is not None:
		error_outfile = error_outdir + sample_name + '.txt'
		error_outfile_handle = open(error_outfile, 'w')

	"""
	START WORKFLOW
	"""

	# Step 1: Process NCBI Genbank and (re)create genbank/proteome with updated locus tags.
	gbk_outfile = gbk_outdir + sample_name + '.gbk'
	pro_outfile = pro_outdir + sample_name + '.faa'
	map_outfile = map_outdir + sample_name + '.txt'
	tot_cds_features = 0
	accounted_features = 0
	try:
		gbk_outfile_handle = open(gbk_outfile, 'w')
		pro_outfile_handle = open(pro_outfile, 'w')
		map_outfile_handle = open(map_outfile, 'w')

		locus_tag_iterator = 1		
		oigf = None
		if input_genbank_file.endswith('.gz'):
			oigf = gzip.open(input_genbank_file, 'rt')
		else:
			oigf = open(input_genbank_file)

		previously_accounted_lts = set([])

		for rec in SeqIO.parse(oigf, 'genbank'):
			updated_features = []
			for feature in rec.features:
				if feature.type == "CDS":
					tot_cds_features += 1
					
					start, end, direction, all_coords = util.process_location_string(str(feature.location))	# type: ignore

					old_locus_tag = None
					prot_seq = None
					try:
						old_locus_tag = feature.qualifiers.get('locus_tag')[0]
					except Exception as e:
						if error_no_lt:
							msg = f'Error: A CDS does not have a locus tag. Please check the Genbank file for sample {sample_name}.'
							if error_outfile_handle is not None:
								error_outfile_handle.write(msg + '\n') # type: ignore
							raise RuntimeError(msg)
						else:
							if not rename_all_flag:
								msg = "Warning: A CDS does not have a locus tag. Will give it an arbitary one."
								sys.stderr.write(msg + '\n')

					try:
						prot_seq = str(feature.qualifiers.get('translation')[0]).replace('*', '')
					except Exception as e:
						if error_no_translation:
							msg = f'Error: A CDS does not have a translation. Please check the Genbank file for sample {sample_name}.'
							if error_outfile_handle is not None:
								error_outfile_handle.write(msg + '\n') # type: ignore
							raise RuntimeError(msg)
						else:
							msg = "Warning: A CDS does not have a translation. Skipping it."
							sys.stderr.write(msg + '\n')
							continue

					new_locus_tag = None
					if old_locus_tag != None and not rename_all_flag:
						new_locus_tag = old_locus_tag
					elif rename_all_flag or old_locus_tag == None:
						new_locus_tag = locus_tag + '_'
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
						locus_tag_iterator += 1
						feature.qualifiers['locus_tag'] = new_locus_tag
					
					if new_locus_tag in previously_accounted_lts:
						msg = f'Error: Locus tag {new_locus_tag} exists multiple times in the Genbank file for sample {sample_name}.'
						if error_outfile_handle is not None:
							error_outfile_handle.write(msg + '\n') # type: ignore
						raise RuntimeError(msg)
					
					pro_outfile_handle.write('>' + str(new_locus_tag) + ' ' + rec.id + ' ' + str(start) + ' ' + str(end) + ' ' + str(direction) + '\n' + prot_seq + '\n')
					if old_locus_tag == None:
						old_locus_tag = 'NA'
					map_outfile_handle.write(str(old_locus_tag) + '\t' + str(new_locus_tag) + '\n')
					accounted_features += 1
					previously_accounted_lts.add(new_locus_tag)
					updated_features.append(feature)
				else:
					updated_features.append(feature)
			rec.features = updated_features
			SeqIO.write(rec, gbk_outfile_handle, 'genbank')
		oigf.close()
		gbk_outfile_handle.close()
		pro_outfile_handle.close()
		map_outfile_handle.close()

	except Exception as e:

		try:
			os.remove(gbk_outfile)
		except OSError as error: 
			sys.stderr.write(f"Error removing files: {error}")
		try:
			os.remove(pro_outfile)
		except OSError as error: 
			sys.stderr.write(f"Error removing files: {error}")
		try:
			os.remove(map_outfile)
		except OSError as error: 
			sys.stderr.write(f"Error removing files: {error}")

		msg = f"Issue processing NCBI Genbank file for sample {sample_name} for some non-evident reason -\n"
		msg += f"please consider reporting and sharing your GenBank file with us on GitHub issues so we can\n"
		msg += f"improve.\n"
		if error_outfile_handle is not None:
			error_outfile_handle.write(msg + '\n') # type: ignore
		else:
			sys.stderr.write(traceback.format_exc() + '\n')
		raise RuntimeError(msg + '\n')

	if accounted_features/tot_cds_features < 0.9:
		try:
			os.remove(gbk_outfile)
		except OSError as error: 
			sys.stderr.write(f"Error removing files: {error}")
		try:
			os.remove(pro_outfile)
		except OSError as error: 
			sys.stderr.write(f"Error removing files: {error}")
		try:
			os.remove(map_outfile)
		except OSError as error: 
			sys.stderr.write(f"Error removing files: {error}")

		msg = f'Processed only {accounted_features} features of {tot_cds_features} total features for sample\n'
		msg += f'{sample_name} - so removing this file from database inclusion.\n'
		if error_outfile_handle is not None:
			error_outfile_handle.write(msg + '\n') # type: ignore

		raise RuntimeError(msg)

	if error_outfile_handle is not None:
		error_outfile_handle.close()

	sys.exit(0)

if __name__ == '__main__':
	processAndReformatNCBIGenbanks()