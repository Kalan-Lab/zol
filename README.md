# *zol (& fai)*

zol (& fai) are tools to zoom in on a locus and perform comparative genomics (uh genetics) between homologous instances of gene clusters (not just BGCs, but phages and ICEs too!)

zol and fai are similar to functionalities offered in lsaBGC but are designed to broadly look at similarities between a variety of gene clusters.

The tools are also designed to explore a single set of homologous gene-clusters (e.g. a single GCF), not the BGC-ome of an entire genus or species like lsaBGC.

Some features in zol reports are more up to date than lsaBGC (but we plan to incorporate these in future versions of lsaBGC). This includes HyPhy-based selection inference and gene tree to gene-cluster tree concordance statistics.

### Zoom on Locus (zol) and Find Additional Instances (fai)

**zol** is a program to create table reports showing homolog group conservation, annotation, and evolutionary stats for any gene-cluster or locus of interest (works for eukaryotes, but designed for bacteria).

**fai** is a program to search for additional instances of a gene-cluster or genome locus in some set of genomes (bacteria specific, soon to work for eukaryotes). Inspired by cblaster (in concept) and ClusterFinder (in algorithm). Works similar to our lsaBGC-Expansion.py program described in the lsaBGC manuscript. Users would need to prepare a set of genomes for searching using **prepTG**. It aims to find really similar instances of gene-clusters in a rapid and automated fashion, for more divergent searches, we recommend users check out [cblaster](https://github.com/gamcil/cblaster), which offers a more interactive experience for users to select desired gene-clusters.

### Installation:

```bash
# 1. clone Git repo and cd into it!
git clone https://github.com/Kalan-Lab/zol
cd zol/

# 2. create conda environment using yaml file and activate it!
conda env create -f zol_env.yml -p /path/to/zol_conda_env/
conda activate /path/to/zol_conda_env/

# 3. complete python installation with the following commands:
python setup.py install
pip install -e .

# 4. depending on internet speed, this can take 20-30 minutes
# within zol Git repo with conda environment activated, run:
setup_annotation_dbs.py
```

### Usage:

### prepTG (preparing to run fai)

prepTG formats and parses information in provided GenBank files or can run prodigal (for bacteria only!) for gene-calling if provided FASTA files and subsequently create GenBank files.

```bash
prepTG -i Folder_with_Target_Genomes/ -o prepTG_DB/
```

For additoinal details on prepTG (e.g. how to download genomes from NCBI), please check out the [1. more info on prepTG]() wiki page.

### fai (finding homologous instances)

fai uses an HMM-based approach to identify homologous instances of a gene-cluster or known set of homologous gene-clusters using an approach very much analogous to lsaBGC-Expansion.py. It is more general, has flexibility in input types, and offers additional parameters/conditions for user adjustment:

1. Provide GenBank(s) of known instance(s) of gene cluster in an input directory

```bash
fai -i Known_GeneCluster_Genbanks/ -sg Search_Genomes/ -tg prepTG_Database/ -o fai_Results/
```

2. Provide gene-cluster coordinates along a FASTA reference genome 

```bash
fai -r Reference.fasta -rc scaffold01 -rs 40201 -re 45043 -tg prepTG_Database/ -o fai_Results/
```

3. Provide proteins gene-cluster using set of proteins that should be co-clustered (similar to cblaster!)

```bash
fai -pq Gene-Cluster_Query_Proteins.faa -tg prepTG_Database/ -o fai_Results/
```

For additional details on fai (e.g. how it relates to cblaster and lsaBGC-Expansion, plots it can create to assess homologous gene-clusters detected), please check out the [2. more info on fai']() wiki page.

### zol (generating table reports)

```bash
zol.py -i Genbanks_Directory/ -g Genomes_Directory/ -o Results/
```

zol produces an xlsx spreadsheet report similar to lsaBGC-PopGene.py where rows correspond to each individual orthogroup/homolog-group and columns provide basic stats, consensus order, annotation information using multiple databases, and evolutionary/selection-inference statistics. Coloring is automatically applied on select quantitative field for users to more easily assess trends.

Annotation databases include: KEGG, NCBI's PGAP, PaperBLAST, VOGs (phage gene clusters), MI-BiG (characterized BGCs), VFDB (virulence factors), CARD (antibiotic resistance), ISfinder (transposons/insertion-sequences).

zol also produces a heatmap and can also be run in dereplication mode to obtain a diverse and representative set of GenBanks/gene-clusters for visual exploration with the wonderful [clinker](https://github.com/gamcil/clinker) tool. 

For details on the stats/annotations zol infers, please refer to the [3. more info on zol](https://github.com/Kalan-Lab/zol/wiki/3.-more-info-on-zol/) wiki page.


![image](https://user-images.githubusercontent.com/4260723/214199267-81546d36-c98d-4394-a146-f9679e0628fe.png)


## Dependencies and Citation

#### Manuscript in preparation! Please cite this GitRepo in the meantime if you find it useful!

**Please consider citing the following dependencies!**
* **MUSCLE5** for performing multiple sequence alignments and PAL2NAL for converting to codon alignments.
* **DIAMOND** for alignments in determining homolog groups and **FastTree2** for subsequent phylogeny construction.
* **HyPhy** and **FASTME** for selection analyses.
* **FastANI** for dereplication of gene-clusters/GenBanks.
* **antiSMASH, GECCO, DeepBGC, VIBRANT**, or **ICEfinder** if you used to identify a BGC, phage, or ICEs.
* **PFAM, KEGG, NCBI's PGAP, MIBIG, VOG, VFDB, CARD,** and **ISFinder** databases used for annotation. 
* **lsaBGC, BiG-SCAPE/CORASON, cblaster**, or **BiG-SLICE** studies if you used them to identify homologous BGCs.

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
