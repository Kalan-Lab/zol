# *zol (& fai)*

zol (& fai) are tools to zoom in on a locus and perform comparative genomics (uh genetics) between homologous instances of gene clusters (not just BGCs, but phages and ICEs too!)

zol and fai are similar to functionalities offered in [lsaBGC](https://github.com/Kalan-Lab/lsaBGC) but are designed to broadly look at similarities between a variety of gene clusters.

The tools are also designed to explore a single set of homologous gene-clusters (e.g. a single GCF), not the BGC-ome of an entire genus or species like [lsaBGC](https://github.com/Kalan-Lab/lsaBGC).

zol produces a basic heatmap, but for visualizations of gene-clusters we recommend other tools such as [clinker](https://github.com/gamcil/clinker) and [gggenomes](https://github.com/thackl/gggenomes), which we think the in-depth spreadsheet complements nicely.

Some features in zol reports are more up to date than lsaBGC (but we plan to incorporate these in future versions of lsaBGC), including HyPhy-based site-specific selection inferences.

![overview_for_github](https://user-images.githubusercontent.com/4260723/229955061-64d32132-f016-485d-9cd9-8493f2de5dd5.png)

### Zoom on Locus (zol) and Find Additional Instances (fai)

**`zol`** is a program to create table reports showing homolog group conservation, annotation, and evolutionary stats for any gene-cluster or locus of interest (works for eukaryotes, but designed for bacteria).

**`fai`** is a program to search for additional instances of a gene-cluster or genome locus in some set of genomes (bacteria specific, soon to work for eukaryotes). Inspired by cblaster (in concept) and ClusterFinder (in algorithm). Works similar to our lsaBGC-Expansion.py program described in the lsaBGC manuscript. Users would need to prepare a set of genomes for searching using **`prepTG`**. It aims to find really similar instances of gene-clusters in a rapid and automated fashion, for more divergent searches, we recommend users check out [cblaster](https://github.com/gamcil/cblaster), which offers a more interactive experience for users to select desired gene-clusters. 
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
# end product will be 28 GB! but can also run in minimal mode
# (which will only download PGAP HMM models < 5 GB) using -m.
# within zol Git repo with conda environment activated, run:
setup_annotation_dbs.py
```

### Test case:

Following installation, you can run a provided test case focused on a subset of Enterococcal polysaccharide antigen instances in *E. faecalis* and *E. faecium* as such:

```bash
bash run_tests.sh
```

### Usage:

### prepTG (preparing to run fai)

prepTG formats and parses information in provided GenBank files or can run prodigal (for bacteria only!) for gene-calling if provided FASTA files and subsequently create GenBank files.

```bash
prepTG -i Folder_with_Target_Genomes/ -o prepTG_DB/
```

For additoinal details on prepTG (e.g. how to download genomes from NCBI), please check out the [1. more info on prepTG](https://github.com/Kalan-Lab/zol/wiki/1.-more-info-on-prepTG) wiki page.

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

For additional details on fai (e.g. how it relates to cblaster and lsaBGC-Expansion, plots it can create to assess homologous gene-clusters detected), please check out the [2. more info on fai](https://github.com/Kalan-Lab/zol/wiki/2.-more-info-on-fai) wiki page.

### zol (generating table reports)

```bash
zol.py -i Genbanks_Directory/ -g Genomes_Directory/ -o Results/
```

zol produces an xlsx spreadsheet report where rows correspond to each individual orthogroup/homolog-group and columns provide basic stats, consensus order, annotation information using multiple databases, and evolutionary/selection-inference statistics. Coloring is automatically applied on select quantitative field for users to more easily assess trends.

Annotation databases include: KEGG, NCBI's PGAP, PaperBLAST, VOGs (phage related genes), MI-BiG (genes from characterized BGCs), VFDB (virulence factors), CARD (antibiotic resistance), ISfinder (transposons/insertion-sequences).

For details on the stats/annotations zol infers, please refer to the [zol](https://github.com/Kalan-Lab/zol/wiki/3.-more-info-on-zol/) wiki page.

![image](https://user-images.githubusercontent.com/4260723/229951285-787042d3-d93b-43d8-b897-63c10d3d9a1a.png)

## Dependencies and Citation

#### Manuscript in preparation! Please cite this GitRepo in the meantime if you find it useful!

**Please consider citing the following dependencies!**
* **pyrodigal**, **prodigal**, and **miniprot** for gene-calling/mapping.
* **MUSCLE5** for performing multiple sequence alignments and PAL2NAL for converting to codon alignments.
* **DIAMOND** for alignments in determining homolog groups and **FastTree2** for subsequent phylogeny construction.
* **CD-HIT** for query protein clustering in fai and 're-inflation' approach in zol.
* **HyPhy** and **FASTME** for selection analyses.
* **skani** for dereplication of gene-clusters/GenBanks.
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

### Version 1.02

- prepTG sample to GenBank relations now specified locally so creation of database is not locked into one location.
- Individual pickle files produced by prepTG per genome/metagenome for lower memory use with fai.
- New "Gene-Clumper" mode for gene-cluster discovery in fai, which is now the default.
- Fixed bug pertaining to overlap between merged gene-clusters based on `--max_gene_disconnect` parameter when using "HMM" mode.
- Improved filtering and retention of GenBanks in zol.
- Fixed bug in re-inflation method in zol.

### Version 1.01

- Remove unused individual proteome files in prepTG database directory.
- Store only gene-location information for scaffolds with hits by query proteins in fai to keep memory use low.
- Introduce parallelization to HMM step of fai and use global variables to access common data without duplicating in memory.
- Improve parsing of different input formats for fai and generate new PDF at end mapping individual protein names to non-redundantified protein queries.
- Declare "< 3 segrating sites found" as reason for inability to calculate Tajima's D instead of just "NA", which could also arise from not enough sequences or the sequence length threshold being met. 
