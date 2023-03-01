import os
import sys
import pandas
import argparse
from time import sleep

zol_main_directory = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
plot_prog = zol_main_directory + 'scripts/generateSyntenicVisual.R'

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
	except:
		sys.stderr.write('Error validating input file %s exists!' % input_zol_report)

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n")
		sleep(5)
	else:
		os.mkdir(outdir)

	df = pandas.read_csv(input_zol_report, sep='\t', header=0)
	"""
	0	Homolog Group (HG) ID
	1	HG is Single Copy?
	2	Proportion of Total Gene Clusters with HG
	3	HG Median Length (bp)
	4	HG Consensus Order
	5	HG Consensus Direction
	6	Proportion of Focal Gene Clusters with HG
	7	Proportion of Comparator Gene Clusters with HG
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
	28	MI-BiG Annotation (E-value)
	29	VOG Annotation (E-value)
	30	VFDB Annotation (E-value)
	31	Pfam Domains
	32	CDS Locus Tags
	33	HG Consensus Sequence
	"""

	df.sort_values(by="HG Consensus Order", ascending=True, inplace=True)

	prev_end = 1
	plot_input_file = outdir + 'Track.txt'
	pif_handle = open(plot_input_file, 'w')
	pif_handle.write('\t'.join(['HG', 'Start', 'End', 'Direction', 'SC', 'Metric']) + '\n')
	for index, row in df.iterrows():
		hg = row['Homolog Group (HG) ID']
		hg_cons = row['Proportion of Total Gene Clusters with HG']
		if float(hg_cons) < 0.25: continue
		hg_mlen = float(row['HG Median Length (bp)'])
		hg_dir = row['HG Consensus Direction']
		sc_flag = row['HG is Single Copy?']
		sc_mark = ''
		if sc_flag == False:
			sc_mark = '//'

		start = prev_end
		end = prev_end + hg_mlen
		dir = 1
		if hg_dir == '-':
			dir = 0
		metric_val = row[metric_name]
		print_row = [hg, start, end, dir, sc_mark, metric_val]
		pif_handle.write('\t'.join([str(x) for x in print_row]) + '\n')
		prev_end = end + 200
	pif_handle.close()

	plot_result_pdf = outdir + 'Plot.pdf'
	plot_cmd = ['Rscript', plot_prog, plot_input_file, str(plot_height), str(plot_width), plot_result_pdf]
	try:
		subprocess.call(' '.join(plot_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		assert (os.path.isfile(plot_result_pdf))
	except Exception as e:
		sys.stderr.write('Had an issue running R based plotting: %s\n' % ' '.join(plot_cmd))
		sys.exit(1)

if __name__ == '__main__':
	genSynVis()
