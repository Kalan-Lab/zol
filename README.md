# *zol (& fai)*

zol (& fai) are tools to zoom in on a locus and perform comparative genomics (uh genetics) between homologous instances of gene clusters (not just BGCs, but phages and ICEs too!)

zol and fai are similar to functionalities offered in lsaBGC but are designed to broadly look at similarities between a variety of gene clusters.

The tools are also designed to explore a single set of homologous gene-clusters (e.g. a single GCF), not the BGC-ome of an entire genus or species like lsaBGC.

Some features in zol reports are more up to date than lsaBGC (but we plan to incorporate these in future versions of lsaBGC). This includes HyPhy-based selection inference and gene tree to gene-cluster tree concordance statistics.

### Zoom on Locus (zol) and Find Additional Instances (fai)

**zol** is a program to create table reports showing homolog group conservation, annotation, and evolutionary stats for any gene-cluster or locus of interest (works for eukaryotes, but designed for bacteria).

**fai** is a program to search for additional instances of a gene-cluster or genome locus in some set of genomes (bacteria specific, soon to work for eukaryotes). Inspired by cblaster (in concept) and ClusterFinder (in algorithm). Works similar to our lsaBGC-Expansion.py program described in the lsaBGC manuscript.

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

# 4. depending on internet speed can take 20-30 minutes
# within zol Git repo with conda environment activated, run:
setup_annotation_dbs.py
```

### Usage:

### fai (finding homologous instances; works for Prokaryotes **currently, stay tuned!**)

fai uses an HMM based approach to identify homologous instances of a gene-cluster or known set of homologous gene-clusters using an approach very much analogous to lsaBGC-Expansion.py. It is more general and has flexibility in input types:

1. Provide GenBank(s) of known instance(s) of gene cluster in an input directory

```commandline
fai -i Known_GeneCluster_Genbanks/ -sg Search_Genomes/ -o fai_Results/
```

2. Provide gene-cluster coordinates along a FASTA reference genome 

```commandline
fai -r Reference.fasta -rc scaffold01 -rs 40201 -re 45043 -o fai_Results/
```

3. Provide proteins gene-cluster using set of proteins that should be co-clustered (similar to cblaster!)

```commandline
fai -pq Gene-Cluster_Query_Proteins.faa -o fai_Results/
```


#### Remember: ncbi-genome-download is a great way to download genomes (is included in conda-environment for zol)

Github link: https://github.com/kblin/ncbi-genome-download

```bash
# Example for downloading all Streptomyces genomes in RefSeq
ncbi-genome-download -F fasta -s refseq -g "Streptomyces" --flat-output -o RefSeq_Streptomyces_Genomes/ bacteria -p 4
```

#### Relation to cblaster and lsaBGC-Expansion

Similar to lsaBGC-Expansion it is inspired by cblaster and ClusterFinder and similar methods. 

Similar to lsaBGC-Expansion - fai aims to work on fragmented genomes and potentially MAGs and can report multiple segments, provided they are on scaffold/contig edges. This is not the default for fai, where `--draft-mode` needs to be specified to initiate it. 

Unlike cblaster, it can also filter for syntenic similarity with regards to gene order and direction in the reference gene cluster. It does this through comparing candidate homologous gene-cluster segments to known reference gene-clusters provide as the queries via Pearson correlation of gene-coordinate midpoints that are facing in the same relative directions (this is the same approach used in lsaBGC-Expansion).

cblaster uniquely offers the ability to perform gene cluster searches remotely, in case you lack computational resources and are interested in searching a large set of genomes.

#### Using cblaster instead of fai to generate input for zol

It is possible to use cblaster (Gilchrist et al. 2021) instead of fai to generate the inputs required for zol.

For instance:

```bash
# use cblaster to search for homologous co-clusters in NCBI genomes
cblaster search -qf queries.faa -s cblaster_results.json

# use cblaster to extract GenBannks of homologous gene-clusters detected
cblaster extract_clusters session.json -o example_directory/
```

### zol (generating table reports; works for Eukaryotes + Prokaryotes)

```commandline
zol.py -i Genbanks_Directory/ -g Genomes_Directory/ -o Results/
```

zol produces a table report similar to lsaBGC-PopGene.py where rows correspond to each individual orthogroup/homolog-group and columns provide basic stats, consensus order, annotation information using multiple databases, and evolutionary/selection-inference statistics.

Annotation databases include KEGG, NCBI's PGAP, PaperBLAST, VOGs (phage gene clusters), MI-BiG (characterized BGCs), VFDB (virulence factors), CARD (antibiotic resistance), ISfinder (transposons/insertion-sequences)

zol also produces a heatmap PDF as well, because visuals are always nice too!

## Showcase Example(s)

### Finding new auxiliary genes of the cyphomycin synthesis encoding BGC in ***Streptomyces***

At the Natural Products (Meta)Genome Mining symposium in Hellerup, Denmark in 2022, Marnix Medema made a great suggestion for the utility of the lsaBGC-Expansion algorithm for discovering new genes when discovering new instances of GCFs in genomes. This was difficult to implement in lsaBGC because the homolog group set is hard-locked at this stage in the lsaBGC workflow; however, the combination of fai & zol make this possible to do and we show such an application here by first identifying homologous instances of the BGC encoding for the synthesis of cyphomycin described by Chevrette et al. 2019.

Note, `examples/cyphomycin/` can be found in the main zol GitHub directory. Please change directories there.

```bash
# change directory to example folder
cd /path/to/zol/examples/cyphomycin/
```

**1. Run fai to find homologous BGCs to the cyphomycin encoding BGC in other *Streptomyces* genomes**

```bash
fai -i Reference_GBKs/ -tg Target_Genomes/ -o fai_Results/ -e 1e-10 -m 30 -kpq Key_Proteins.faa -kpe 1e-20 -kpm 5 -c 4
```

So some parameters are being used here. First we are providing a directory with the reference BGC for cyphomycin downloaded from MI-BiG with the `-i` argument. 

Next, we are specifying a directory of genomes to search, where genomes are in FASTA format, using the `-tg` argument. Note, to prevent putting all Streptomyces genomes in this folder, I preselected the genomes where I found homologous instances of the reference BGC. Some manual investigation was done to ensure final sets were to the desired level of similarity.

`-o` simply specifies the output directory.

`-e` is the E-value threshold to use for finding homologs with DIAMOND in the target genomes. 

`-m` is the minimum number of homolog groups/genes in a segment the HMM identifies as in a "homologous" state to consider the BGC as present.

`-kpq` is a FASTA file containing proteins for key proteins in the gene cluster. Here we provide sequences for 7 of the largest genes in the reference BGC.

`-kpm` is the minimum number of key proteins in a candidate segment that are needed to be reported.

`-kpe` is the minimum E-value threshold needed for the key proteins. Because there is a smaller set, you might want to increase stringency for them to prevent false positives.

`-c` is the number of cpus/threads to use.

After running, homologous gene-clusters to the reference can be found in the subdirectory of the results at: `fai_Results/Additional_Gene_Cluster_Instances_Identified/`.

**2. Perform comparative analytics using zol**

```bash
# first copy over the reference and additional homologous instances into a single directory:
mkdir zol_Input/
cp examples/Reference_GBKs/* fai_Results/Additional_Gene_Cluster_Instances_Identified/*  zol_Input/

# now run zol
zol -i zol_Input/ -o zol_Results/ -c 4 
```
After running, final results can be found at: `zol_Results/`

Including the following static PDF: 



**3. Filter for genes which are missing in cyphomycin reference BGC but found in homologous BGCs**

**4. Filter for genes which are unique to the cyphomycin reference BGC but missing in homologous BGCs**

**5. Investigate cluster further using clinker**

clinker is an amazing program by Gilchirst & Chooi 2021 which produces interactive visualizations. We integrate it here and can be further used to explore homologous gene-clusters including syntenic conservation.


### Application to finding and investigating homologous phages and ICEs (coming soon!)

## Dependencies and Citation

#### Manuscript in preparation! Please cite this GitRepo in the meantime if you find it useful!

**Please consider citing the following dependencies!**
* **MUSCLE5** for performing multiple sequence alignments and PAL2NAL for converting to codon alignments.
* **Treemmer** and **clinker** for visualization
* **DIAMOND** for alignments in determining homolog groups and **FastTree2** for subsequent phylogeny construction.
* **HyPhy** and **FASTME** for selection analyses.
* **STAG** for consensus tree construction used in gene tree congruence statistic.
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
