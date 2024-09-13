#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO

for rec in SeqIO.parse(sys.argv[1], 'genbank'):
    print('>' + rec.id + '\n' + str(rec.seq))
