# *zol (& fai)*

[![Preprint](https://img.shields.io/badge/Preprint-bioRxiv-darkblue?style=flat-square&maxAge=2678400)](https://www.biorxiv.org/content/10.1101/2023.06.07.544063v3)
[![Documentation](https://img.shields.io/badge/Documentation-Wiki-darkgreen?style=flat-square&maxAge=2678400)](https://github.com/Kalan-Lab/zol/wiki)
[![Docker](https://img.shields.io/badge/Docker-DockerHub-darkred?style=flat-square&maxAge=2678400)](https://hub.docker.com/r/raufs/zol)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/zol/README.html) [![Conda](https://img.shields.io/conda/dn/bioconda/zol.svg)](https://anaconda.org/bioconda/zol/files)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/zol/badges/latest_release_date.svg)](https://anaconda.org/bioconda/zol)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/zol/badges/platforms.svg)](https://anaconda.org/bioconda/zol)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/zol/badges/license.svg)](https://anaconda.org/bioconda/zol)

***zol (& fai)*: tools for targeted searching and evolutionary investigations of gene clusters (sets of co-located genes - e.g. biosynthetic gene clusters, viruses/phages, operons, etc.).**

First, fai allows users to search for homologous/orthologous instances of a query gene cluster in a database of (meta-)genomes. There are some other similar tools, including convenient webservers, to fai (which we highlight and recommend as altneratives on [this documentation page](https://github.com/Kalan-Lab/zol/wiki/5.1-tutorial-for-using-zol-with-output-from-fast.genomics-and-CAGECAT)); but, fai also has some unique/rarer options. Mainly, fai pays special attention to see whether gene cluster hits in target (meta-)genomes are on scaffold/contig edges and takes consideration of this, during both detection and downstream assessment. E.g. fai will mark individual coding genes and gene cluster instances if they are on the edge of a scaffold/contig, which can then be used as a filter in zol. *This is important for calculation of conservation of genes across homologous gene clusters!* 

After finding homologous instances of a gene cluster - using fai or other software - users often wish to investigate the similarity between instances. This is often performed using pairwise similarity assessment via visualization with tools such as clinker, gggenomes, etc. While these tools are great, **if you found 100s or 1000s of gene cluster instances** such visualizations can get overwhelming and computationally expensive to render. To simplify the identification of interesting functional, evolutionary, and conservation patterns across 100s to 1000s of homologous gene cluster instances, we developed zol to perform *de novo* ortholog group predictions and create detailed color-formatted XLSX spreadsheets summarizing information. More recently, we have also introduced scalable visualization tools (*cgc & cgcg*) that allow for simpler assessment of information represented across thousands of homologous gene cluster instances.

<p align="center">
<img src="https://github.com/user-attachments/assets/b0ec16bf-f302-4018-a7eb-91ff8a8b7817" width="600">
</p>

### Citation:
> [zol & fai: large-scale targeted detection and evolutionary investigation of gene clusters](https://www.biorxiv.org/content/10.1101/2023.06.07.544063v3). *bioRxiv 2023.* Rauf Salamzade, Patricia Q Tran, Cody Martin, Abigail L Manson, Michael S Gilmore, Ashlee M Earl, Karthik Anantharaman, Lindsay R Kalan

*In addition, please cite important [dependency software or databases](https://github.com/Kalan-Lab/zol/wiki/6.-dependencies) for your specific analysis accordingly.*

> [!CAUTION]
> Please avoid using versions 1.5.1 to 1.5.3 in which zol has the possibility to get stuck in a while loop and write a large file. This issue is resolved in v1.5.4. 

## Main Contents:

1. [Documetation](https://github.com/Kalan-Lab/zol/wiki)
2. [Overview of Major Results](https://github.com/Kalan-Lab/zol/wiki/0.-overview-of-major-result-files)
3. [Short note on resource requirements](#short-note-on-resource-requirements)
4. [Installation](#installation)
5. [Test Case](#test-case)
6. [Example Usages](https://github.com/Kalan-Lab/zol/wiki/4.-basic-usage-examples)
7. [Tutorial with Tips and Tricks](https://github.com/Kalan-Lab/zol/wiki/5.-tutorial-%E2%80%90-a-detailed-walkthrough)

### Auxiliary tools within the suite:
* [abon, atpoc, and apos: Assessing the conservation of a focal sample's BGC-ome, phage-ome, and plasmid-ome](https://github.com/Kalan-Lab/zol/wiki/0.-overview-of-major-result-files#abon-atpoc-and-apos-results)
* [(***New***) cgc: Summary visualization of 1000s of gene clusters](https://github.com/Kalan-Lab/zol/wiki/5.3-visualization-of-1000s-of-gene-clusters-using-cgc)
* [(***New***) cgcg: Network visualization of ortholog groups across 1000s of gene clusters](https://github.com/Kalan-Lab/zol/wiki/5.3-visualization-of-1000s-of-gene-clusters-using-cgc)
* [(***New***) salt: Assessing support for lateral gene transfer](https://github.com/Kalan-Lab/zol/wiki/5.4-horizontal-or-lateral-transfer-assessment-of-gene-clusters-using-salt)


## Short Note on Resource Requirements:

Different programs in the zol suite have different resource requirements. Moving forward, the default settings in the `zol` program itself should usually allow for low memory usage and faster runtime. For thousands of gene cluster instances, we recommend to either use the dereplication/reinflation approach (see manuscript for comparison on evolutionary statistics between this approach and a full processing) or using CD-HIT clustering (a greedy incremental clustering approach - which is nicely illustrated/explained on the [MMSeqs2 wiki](https://github.com/soedinglab/MMseqs2/wiki#clustering-modes)) to determine protein clusters/families (not true ortholog groups). Disk space is generally not a huge concern for zol analysis, but if working with thousands of gene clusters things can temporarily get large. 

Available disk space is the primary concern however for `fai` and `prepTG`. This is mostly the case for users interested in the construction and searching of large databases (containing over a thousand genomes). Generally, `prepTG` and `fai` are designed to work on metagenomic as well as genomic datasets and do not have a high memory usage, but genomic files stack up in space and DIAMOND alignment files can quite get large as well.

## Installation:

#### Bioconda (Recommended):

Note, (for some setups at least) ***it is critical to specify the conda-forge channel before the bioconda channel to properly configure priority and lead to a successful installation.***
 
**Recommended**: For a significantly faster installation process, use `mamba` in place of `conda` in the below commands, by installing `mamba` in your base conda environment.

```bash
# 1. install and activate zol
conda create -n zol_env -c conda-forge -c bioconda zol
conda activate zol_env

# 2. depending on internet speed, this can take 20-30 minutes
# end product will be ~40 GB! You can also run in minimal mode
# (which will only download Pfam & PGAP HMM models ~8.5 GB)
# using the -m argument. 
setup_annotation_dbs.py [-m]
```

> [!NOTE]
> When you create a conda environment using `-n`, the environment will typically be stored in your home directory. However, because the databases can be large, you might prefer to instead setup the conda environment somewhere else with more space on your system using `-p`. For instance, `conda create -p /path/to/drive_with_more_space/zol_conda_env/ -c conda-forge -c bioconda zol`. Then, next time around you would simply activate this environment by providing the path to it: `conda activate /path/to/drive_with_more_space/zol_conda_env/`

#### Docker:

___Requires docker to be installed on your system!___

To keep the Docker image size relatively low (currently ~13 GB), only the Pfam and PGAP HMMs/databases are included.

```bash
# get wrapper script from GitHub
wget https://raw.githubusercontent.com/Kalan-Lab/zol/main/docker/run_ZOL.sh

# change permissions to allow execution
chmod a+x ./run_ZOL.sh

# run script
./run_ZOL.sh
```

## Test case:

Following installation, you can run a provided test case focused on a subset of Enterococcal polysaccharide antigen instances in *E. faecalis* and *E. faecium* as such:

#### Bioconda:

```bash
# download test data tar.gz and bash script for running tests
wget https://github.com/Kalan-Lab/zol/raw/main/test_case.tar.gz
wget https://raw.githubusercontent.com/Kalan-Lab/zol/main/run_tests.sh

# run bash-based testing script
bash run_tests.sh
```

#### Docker:

```bash
# download test scripts from (bash script which you can reference for learning how to run zol).
wget https://raw.githubusercontent.com/Kalan-Lab/zol/main/docker/test_docker.sh

# change permissions to allow execution
chmod a+x ./test_docker.sh

# run tests
./test_docker.sh
```

Note, the script `test_docker.sh` must be run in the same folder as run_ZOL.sh!

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
