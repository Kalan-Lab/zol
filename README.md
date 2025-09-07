# *zol (& fai)*

[![Manuscript](https://img.shields.io/badge/Manuscript-NAR-darkblue?style=flat-square&maxAge=2678400)](https://academic.oup.com/nar/article/53/3/gkaf045/8001966)
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
> [zol & fai: large-scale targeted detection and evolutionary investigation of gene clusters](https://academic.oup.com/nar/article/53/3/gkaf045/8001966). *Nucleic Acids Research* 2025. Rauf Salamzade, Patricia Q Tran, Cody Martin, Abigail L Manson, Michael S Gilmore, Ashlee M Earl, Karthik Anantharaman, Lindsay R Kalan

*In addition, please cite important [dependency software or databases](https://github.com/Kalan-Lab/zol/wiki/6.-dependencies) for your specific analysis accordingly.*

```
Usage: zol-suite [-h] [--list-programs] [--version] <program> ...

The zol suite - a comprehensive bioinformatics toolkit for gene cluster analysis.

/================================\
|| ________      ________       ||
|||\_____  \    |\   ____\      ||
|| \|___/  /|   \ \  \___|_     ||
||     /  / /    \ \_____  \    ||
||    /  /_/__    \|____|\  \   ||
||   |\________\    ____\_\  \  ||
||    \|_______|   |\_________\ ||
||                 \|_________| ||
\================================/

Author: Rauf Salamzade
Lab: Kalan Lab; University of Wisconsin - Madison; McMaster University

This interface provides access to all ZOL tools through a single command-line interface.
Each tool has its own specific arguments and functionality.

Typical order of operations:

1.) Run prepTG to prepare a database of target genomes for searches using fai.
2.) Run fai to find additional instances of a gene cluster of interest in the prepTG database.
3.) Run zol to perform comparative gene cluster analysis on the results from fai.
4.) Run cgc and cgcg to visualize the results from zol.

For help with a specific program, use: zol-suite <program> --help

Positional Arguments:
  <program>        ZOL program to run
    abon           Automated analysis of conservation/novelty for a sample's biosynthetic
                   gene clusters.
    apos           Automated analysis of conservation/novelty for a sample's plasmids.
    atpoc          Automated analysis of conservation/novelty for a sample's prophages.
    cgc            Visualization of zol results along a consensus gene cluster sequence.
    cgcg           Visualization of zol results as a graphical network.
    fai            Find additional instances of gene clusters in a genome database using
                   flexible alignment and synteny criteria.
    prepTG         Prepare a database of target genomes for searches using fai.
    salt           Support assessment for lateral transfer of gene clusters (experimental).
    zol            Perform comparative gene cluster analysis.
    zol-scape      Run zol analysis on BiG-SCAPE results.

Options:
  -h, --help       show this help message and exit
  --list-programs  List all available programs and exit
  --version        show program's version number and exit
  ```

> [!CAUTION]
> Please avoid using versions 1.5.1 to 1.5.3 in which zol has the possibility to get stuck in a while loop and write a large file. This issue is resolved in v1.5.4. 

> [!IMPORTANT]
> We recently updated zol to v1.6.0 - which introduces several key updates, including a unified interface that can be issued as `zol-suite`, improved PEP8 compliance for backend code, and lighter databases constructed using prepTG.

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
* [(***New***) cgc and cgcg: Scalable visualization of 1000s of homologous gene clusters](https://github.com/Kalan-Lab/zol/wiki/5.3-visualization-of-1000s-of-gene-clusters-using-cgc-or-cgcg)

## Short Note on Resource Requirements:

Different programs in the zol suite have different resource requirements. Moving forward, the default settings in the `zol` program itself should usually allow for low memory usage and faster runtime. For thousands of gene cluster instances, we recommend to either use the dereplication/reinflation approach (see manuscript for comparison on evolutionary statistics between this approach and a full processing) or using DIAMOND linclust clustering to determine protein clusters/families (not true ortholog groups). Disk space is generally not a huge concern for zol analysis, but if working with thousands of gene clusters things can temporarily get large. 

Available disk space is the primary concern however for `fai` and `prepTG`. This is mostly the case for users interested in the construction and searching of large databases (containing over a thousand genomes). Generally, `prepTG` and `fai` are designed to work on metagenomic as well as genomic datasets and do not have a high memory usage, but genomic files stack up in space and DIAMOND alignment files can quite get large as well.

## Installation:

#### Bioconda (Recommended):

Note, (for some setups at least) ***it is critical to specify the conda-forge channel before the bioconda channel to properly configure priority and lead to a successful installation.***
 
**Recommended**: For a significantly faster installation process, use `mamba` in place of `conda` in the below commands, by installing `mamba` in your base conda environment.

```bash
# 1. install and activate zol

# On Linux:
conda create -n zol_env -c conda-forge -c bioconda zol
conda activate zol_env

# 2. depending on internet speed, this can take 20-30 minutes
# end product will be ~40 GB! You can also run in minimal mode
# (which will only download Pfam & PGAP HMM models ~8.5 GB)
# using the -m argument. 
setup_annotation_dbs.py [-m]

# 3. run interface program
zol-suite [-h]
```

> [!TIP]
> When you create a conda environment using `-n`, the environment will typically be stored in your home directory. However, because the databases can be large, you might prefer to instead setup the conda environment somewhere else with more space on your system using `-p`. For instance, `conda create -p /path/to/drive_with_more_space/zol_conda_env/ -c conda-forge -c bioconda zol`. Then, next time around you would simply activate this environment by providing the path to it: `conda activate /path/to/drive_with_more_space/zol_conda_env/`

> [!CAUTION]
> If you choose to manually define a database directory for `setup_annotation_dbs.py` by setting the `ZOL_DATA_PATH` environmental variable, make sure that it is a unique directory to zol. This directory will be deleted and recreated when you run the script. You don't have to worry about this if using bioconda where the default directory is located within the conda environment space.

> [!NOTE]
> 🍎 For Mac users with Apple Silicon chips, you might need to specify `CONDA_SUBDIR=osx-64` prior to `conda create` as described [here](https://github.com/bioconda/bioconda-recipes/issues/41702#issuecomment-1614173361). So you would issue: `CONDA_SUBDIR=osx-64 conda create -n zol_env -c conda-forge -c bioconda zol`. You might get warnings still related to some Intel related messages - but these should be ok to ignore mostly - if you have concerns - please feel free to just reach out. 

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

## Developer Installation

A YAML file is included which can be used to setup a local environment.

```bash
# clone git repo (from your own fork!)
git clone https://github.com/Kalan-Lab/zol
cd zol

# create conda from yaml file and activate it
conda env create -f zol_env.yaml -p ../conda_env/
conda activate ../conda_env/

# within cloned git repo, run pip install
pip install .

# deactivate and re-activate to ensure proper
# setting of environment variables
conda deactivate
conda activate ../conda_env/
```

Note, this will be using bleeding-edge code that might not yet be released.

You can now make changes, test them, commit, and finally make a pull request.

Please also consider downloading the full testing dataset and running more comprehensive testing:

```bash
wget https://zenodo.org/records/15934591/files/testing_cases.tar.gz?download=1
mv testing_cases* testing_cases.tar.gz 

# run comprehensive_tests.sh shell script found in zol's git repo:
bash comprehensive_tests.sh
```

## License:

```
BSD 3-Clause License

Copyright (c) 2023-2025, Kalan-Lab
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
