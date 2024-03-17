import os
import sys
from Bio import SeqIO
import argparse

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: extractProteinsFromGenBank.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
    
    A simple script 
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--gene_cluster_genbank', help='Path to gene cluster genbank.', required=True)
	args = parser.parse_args()
	return args

myargs = create_parser()
gbk_file = myargs.gene_cluster_genbank

counter = 0
with open(gbk_file) as ogbf:
    for rec in SeqIO.parse(ogbf, 'genbank'):
        for feat in rec.features:
            if not feat.type == 'CDS': continue
            protein_id = None
            try:
                protein_id = feat.qualifiers.get('locus_tag')[0]
            except:
                pass
            
            if protein_id == None:
                try:
                    protein_id = feat.qualifiers.get('protein_id')[0]
                except:
                    pass
            if protein_id == None:
                protein_id = 'CDS_' + str(counter)
                counter += 1
            print('>' + protein_id + '\n' + str(feat.qualifiers.get('translation')[0]))