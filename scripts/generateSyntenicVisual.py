#!/usr/bin/env python3

"""
Program: generateSyntenicVisual.py
Author: Rauf Salamzade
Kalan Lab
UW Madison, Department of Medical Microbiology and Immunology
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
import pandas
import argparse
import subprocess
from time import sleep
from zol import util
import traceback

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: generateSyntenicVisual.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--zol_report_tsv', help='Path to zol tsv report.', required=True)
	parser.add_argument('-o', '--output_dir', help='Path to output directory.', required=True)
	parser.add_argument('-m', '--metric', help='Metric to use for coloring in plot. Valid options are headers for evolutionary statistics.\nPlease surround by quotes if there is a space in the header of the column. Default is "Tajima\'s D"', required=False, default="Tajima's D")
	parser.add_argument('-l', '--length', type=int, help='Specify the height/length of the plot. Default is 7.', required=False, default=3)
	parser.add_argument('-w', '--width', type=int, help='Specify the width of the plot. Default is 10.', required=False, default=7)

	args = parser.parse_args()
	return args

def genSynVis():
	"""
	PARSE INPUTS
	"""
	myargs = create_parser()

	input_zol_report = myargs.zol_report_tsv
	metric_name = myargs.metric
	outdir = os.path.abspath(myargs.output_dir) + '/'
	plot_height = myargs.length
	plot_width = myargs.width

	try:
		assert(os.path.isfile(input_zol_report))
	except Exception as e:
		sys.stderr.write(f'Error validating input file {input_zol_report} exists!')

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n")
		sleep(5)
	else:
		os.mkdir(outdir)

	df = pandas.read_csv(input_zol_report, sep='\t', header=0)
	"""
	0	ortholog group (OG) ID
	1	OG is Single Copy?
	2	Proportion of Total Gene Clusters with OG
	3	OG Median Length (bp)
	4	OG Consensus Order
	5	OG Consensus Direction
	6	Proportion of Focal Gene Clusters with OG
	7	Proportion of Comparator Gene Clusters with OG
	8	Fixation Index
	9	Upstream Region Fixation Index
	10	Tajima's D
	11	Proportion of Filtered Codon Alignment is Segregating Sites
	12	Entropy
	13	Upstream Region Entropy
	14	Median Beta-RD-gc
	15	Max Beta-RD-gc
	16	Proportion of sites which are highly ambiguous in codon alignment
	17	Proportion of sites which are highly ambiguous in trimmed codon alignment
	18	GARD Partitions Based on Recombination Breakpoints
	19	Number of Sites Identified as Under Positive or Negative Selection by FUBAR
	20	Average delta(Beta, Alpha) by FUBAR across sites
	21	Proportion of Sites Under Selection which are Positive
	22	Custom Annotation (E-value)
	23	KO Annotation (E-value)
	24	PGAP Annotation (E-value)
	25	PaperBLAST Annotation (E-value)
	26	CARD Annotation (E-value)
	27	IS Finder (E-value)
	28	MIBiG Annotation (E-value)
	29	VOG Annotation (E-value)
	30	VFDB Annotation (E-value)
	31	Pfam Domains
	32	CDS Locus Tags
	33	OG Consensus Sequence
	"""

	df.sort_values(by="OG Consensus Order", ascending=True, inplace=True)

	prev_end = 1
	plot_input_file = outdir + 'Track.txt'
	pif_handle = open(plot_input_file, 'w')
	pif_handle.write('\t'.join(['OG', 'Start', 'End', 'Direction', 'SC', 'Metric']) + '\n')
	for index, row in df.iterrows():
		og = row['Ortholog Group (OG) ID']
		og_cons = row['Proportion of Total Gene Clusters with OG']
		#if float(og_cons) < 0.25: continue
		og_mlen = float(row['OG Median Length (bp)'])
		og_dir = row['OG Consensus Direction']
		sc_flag = row['OG is Single Copy?']
		sc_mark = ''
		if sc_flag == False:
			sc_mark = '//'

		start = prev_end
		end = prev_end + og_mlen
		dir = 1
		if og_dir == '-':
			dir = 0
		metric_val = row[metric_name]
		print_row = [og, start, end, dir, sc_mark, metric_val]
		pif_handle.write('\t'.join([str(x) for x in print_row]) + '\n')
		prev_end = end + 200
	pif_handle.close()

	plot_result_pdf = outdir + 'Plot.pdf'
	rscript_path = outdir + 'generateSyntenicVisual.R'
	util.generate_syntenic_visual_r(plot_input_file, plot_result_pdf, plot_height, plot_width, rscript_path)
	plot_cmd = ['Rscript', rscript_path]
	try:
		subprocess.call(' '.join(plot_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		assert (os.path.isfile(plot_result_pdf))
	except Exception as e:
		sys.stderr.write(f"Had an issue running R based plotting - potentially because of R setup issues in conda: {' '.join(plot_cmd)}\n")
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
		
if __name__ == '__main__':
	genSynVis()
