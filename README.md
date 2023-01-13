# *zol*

### Zoom on Locus (zol) and Find Additional Instances (fai)

**zol** is a program to create table reports showing homolog group conservation, annotation, and evolutionary stats for any gene-cluster or locus of interest (works for eukaryotes, but designed for bacteria).

**fai** is a program to search for additional instances of a gene-cluster or genome locus in some set of genomes (bacteria specific). Inspired by cblaster (in concept) and ClusterFinder (in algorithm). Works similar to our lsaBGC-Expansion.py program described in the lsaBGC manuscript.

### Installation:



### Usage:


### fai (finding homologous instances)



### zol (generating table reports)

```commandline
zol.py -i Genbanks_Directory/ -g Genomes_Directory/ -o Results/
```

Please also cite:

* MUSCLE5 for performing multiple sequence alignments and PAL2NAL for converting to codon alignments.
* clinker and treemmer for visualization
* DIAMOND for alignments in determining homolog groups and FastTree2 for subsequent phylogeny construction.
* HyPhy and FASTME for selection analyses.
* STAG for consensus tree construction used in gene tree congruence statistic (and in the near future visualizations).
* antiSMASH, GECCO, DeepBGC, VIBRANT, or ICEfinder if you used to identify a BGC, phage, or ICEs.
* PFAM, KEGG, NCBI's PGAP, MIBIG, VOG, VFDB, CARD, and ISFinder databases used for annotation. 
* lsaBGC, BiG-SCAPE/CORASON, cblaster, or BiG-SLICE studies if you used them to identify homologous BGCs.

## License:

```
BSD 3-Clause License

Copyright (c) 2023, Kalan-Lab
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## Major Updates 
