#!/usr/bin/env python3

"""
Program: extractProteinsFromGenBank.py
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
