=====
*lsa*BGC Suite Overview
=====

*lsa*BGC consists of several individual programs which provide a broad suite of functions for comparative analysis of 
biosynthetic gene clusters across a single focal lineage or taxa (recommended/tested at species or genus levels), to 
understand the allelic variability observed for BGC genes, and mine for novel SNVs within such genes representative
of previously unidentified allelic variants.

![](./images/lsaBGC_Overview.png)

Installation
--------

To learn more about the installation of *lsa*BGC and its dependencies, please take a look at the [Installation](https://github.com/Kalan-Lab/lsaBGC/wiki/01.-Installation) wiki page.

Background / Introduction
--------

What functionalities does *lsa*BGC offer to users? Learn more about the suite's intended usages and where it should not be used, along with recommendations to other great software for exploring and wrangling comparative analysis of secondary metabolite genetic architectures [Background](https://github.com/Kalan-Lab/lsaBGC/wiki/00.-Background) wiki page!

Detailed Walkthrough and Test Cases
--------

A detailed walkthrough for using the lsaBGC suite as intended can be found on the [third Wiki page](https://github.com/Kalan-Lab/lsaBGC/wiki/[03.-Detailed-Walkthrough](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Quick-Start-and-Detailed-Walkthrough)).

We found that the *Corynebacterium kefirresidentii* is a common species complex of the skin microbiome and harbor several BGCs across their compact genome. We use the publicly available genomes from the complex as a small and simple test set to demonstrate the exploratory power of *lsa*BGC. Please have a look at the [lsaBGC_Ckefir_Testing_Cases Github repo](https://github.com/Kalan-Lab/lsaBGC_Ckefir_Testing_Cases) for further details.

Main Programs
--------

*lsa*BGC comprises of 8 primary programs:

Many of the main programs utilize an object oriented infrastructure for processing and analysis. More information on this infrastructure can be found on the wiki page [OOP Framework](https://github.com/Kalan-Lab/lsaBGC/wiki/02.-The-Object-Oriented-Core-of-lsaBGC).

+------------+------------+-----------+-----------+
| Program | Description | Input | Output |
+============+============+============+============+
| [lsaBGC-Ready.py](https://github.com/Kalan-Lab/lsaBGC/wiki/03.-Quick-Start-and-Detailed-Walkthrough) |  Takes existing antiSMASH results (and optionally BiG-SCAPE) and creates inputs necessary to run downstream lsaBGC analyses. | - antiSMASH Results Directory
|                                                                                                      | (optional) BiG-SCAPE Results Directory</li></ul> | <ul><li>OrthoFinder Homolog Group vs. Sample Matrix</li><li>Listing of antiSMASH BGCs</li><li>Listing of Sample Predicted Proteomes/Genbanks - (if BiG-SCAPE results provided) GCF Listings Directory |
|
|
+------------+------------+-----------+-----------+
| [lsaBGC-Cluster.py](https://github.com/Kalan-Lab/lsaBGC/wiki/05.-Clustering-BGCs-into-GCFs) |  Takes the comprehensive list of BGCs and clusters them using MCL into GCFs | - Comprehensive listing of AntiSMASH BGC predictions in Genbank format (from completed/high-quality genomes)</li><li>OrthoFinder Homolog Group vs. Sample Matrix</li></ul> | <ul><li>Summary of GCFs</li><li>Automated report to inform on best clustering parameter choices (if requested)</li><li>List for each GCF of BGC members</li><ul> |
| [lsaBGC-Refiner.py](https://github.com/Kalan-Lab/lsaBGC/wiki/06.-Refinement-of-BGCs-Belonging--to-GCF) | Refines boundaries of BGCs belonging to a single GCF according to user specifications. | <ul><li>BGC instances for focal GCF in Genbank format</li><li>OrthoFinder Homolog Group vs. Sample Matrix</li><li>Boundary Homolog Group ID #1</li><li>Boundary Homolog Group ID #2</li></ul>| <ul><li>BGC instances for focal GCF in Genbank format edited for requested refinement.</li></ul> |
| [lsaBGC-Expansion.py](https://github.com/Kalan-Lab/lsaBGC/wiki/08.-High-throughput-Detection-of-New-GCF-Instances-Across-Draft-Genome-Assemblies) | Uses an HMM based approach to quickly find homologous instances of GCF in draft-quality genomes. | <ul><li>BGC instances for focal GCF in Genbank format</li><li>Additional genomic assemblies listing (post gene-calling)</li></ul> | <ul><li>Expanded list of BGCs belonging to GCF</li><li>Expanded OrthoFinder Homolog Group vs Sample Matrix</li></ul>|
| [lsaBGC-See.py](https://github.com/Kalan-Lab/lsaBGC/wiki/07.-Visualizing-GCFs-Across-Phylogenies) | Visualizes BGC instances of a GCF across a phylogeny |  <ul><li>BGC instances for focal GCF in Genbank format</li><li>(Optional) Species phylogeny</li></ul> | <ul><li>Modified species phylogeny to expand samples which feature multiple BGCs for the GCF (if species phylogeny was provided)</li><li>(Optional) Single-copy-core phylogeny of GCF</li><li>Automated visualization of BGC gene architectures across species or BGC phylogeny in PDF format</li><li>Track file for visualization of gene architecture for BGCs in GCF to be input into iTol.</li></ul>|
| [lsaBGC-Divergence.py](https://github.com/Kalan-Lab/lsaBGC/wiki/09.-Assessing-Evolutionary-Linkage-of-BGCs-with-their-Genome-wide-Contexts) | Determines 𝜷-RT statistic for assessing BGC divergence relative to genome-wide divergence between isolate pairs. | <ul><li>BGC instances for focal GCF in Genbank format</li><li>Pairwise ANI or AAI estimates between samples/genomes with GCF</li></ul> | <ul><li>Report with the 𝜷-RT statistic showcasing the ratio of the genome-wide similarity to the GCF-specific similarity between pairs of isolates with the GCF. </li></ul>|
| [lsaBGC-PopGene.py](https://github.com/Kalan-Lab/lsaBGC/wiki/10.-Population-Genetics-Analysis-of-Genes-Found-in-a-GCF) | Looks at sequence conservation and performs population genetic analyses for each homolog group found in GCF. | <ul><li>BGC instances for focal GCF in Genbank format</li><li>Expanded OrthoFinder Homolog Group vs Sample Matrix</li></ul> | <ul><li>Report with conservation and population-genetic relevant statistic for each homolog group associated with the GCF.</li><li>Automated visualization of genetic variability present in the lineage for each homolog group in PDF format.</li><li>Codon alignment for each homolog group in GCF</li></ul>|
| [lsaBGC-DiscoVary.py](https://github.com/Kalan-Lab/lsaBGC/wiki/11.-Discovering-Novel-Variations-in-GCF-Genes-from-Raw-Sequencing-Reads) | Identifies GCF instances in metagenomes and looks for base-resolution novelty within genes from raw sequencing data not observed in genomic assemblies for the taxonomy. | <ul><li>BGC instances for focal GCF in Genbank format</li><li>Metagenomic/sequencing readsets</li><li>Codon alignments for homolog groups in GCF</li></ul> | <ul><li>Listing of which metagenomic/sequencing readsets are predicted to contain the GCF<li>Table report with novel variants never previously observed in genomic assemblies</li><li>(Optional) Phased homolog group alleles found in metagenomic/sequencing data. [uses DESMAN] </li></ul> | 

Also provided are three workflow/pipeline programs, lsaBGC-AutoProcess.py, lsaBGC-AutoExpansion.py, and lsaBGC-AutoAnalyze.py, which simplify the generation of inputs necessary for the lsaBGC framework and allow for the automatic processing of each GCF post-clustering through standard analysis:

| Program | Description | Input | Output |
| ----------------------------------| ----------- | ------------ |------|
| [lsaBGC-AutoProcess.py](https://github.com/Kalan-Lab/lsaBGC/wiki/04.-Generating-Required-Inputs-for-lsaBGC) | Automatically runs Prokka, AntiSMASH, and OrthoFinder | <ul><li>Genomic assemblies</li></ul> | <ul><li>AntiSMASH BGC predictions in Genbank format</li><li>OrthoFinder Homolog Group vs. Sample Matrix</li></ul> |
| [lsaBGC-AutoExpansion.py](https://github.com/Kalan-Lab/lsaBGC/wiki/08.-High-throughput-Detection-of-New-GCF-Instances-Across-Draft-Genome-Assemblies) | Automatically runs lsaBGC-Expansion for all GCFs and resolves conflicts (e.g. overlapping BGCs for different GCFs) | <ul><li>Directory with BGC listings for each GCF</li><li>Additional genomic assemblies listing (post gene-calling)</li><li>OrthoFinder Homolog Group vs. Sample Matrix</li></ul> | <ul><li>Expanded list of BGCs belonging to GCF</li><li>Expanded OrthoFinder Homolog Group vs Sample Matrix</li></ul>|
| [lsaBGC-AutoAnalyze.py](https://github.com/Kalan-Lab/lsaBGC/wiki/13.-The-lsaBGC-AutoAnalyze-Workflow) | Automatically runs lsaBGC-See.py, lsaBGC-PopGene.py, lsaBGC-Divergence.py, and lsaBGC-DiscoVary for each GCF. | <ul><li>Genomic listing file</li><li>Directory with BGC listings for each GCF</li><li>Additional options</li></ul> | <ul><li>Consolidated reports for lsaBGC-PopGene and lsaBGC-Divergence results</li><li>Visualizations providing overview of lsaBGC analyses</li></ul> |

Future to-do's involve getting these workflows re-written in a DSL framework such as NextFlow.

## Additional Programs / Scripts

Several additional programs and scripts are included in the lsaBGC suite. Major scripts of potential interest are described [here](https://github.com/Kalan-Lab/lsaBGC/wiki/12.-Additional-Programs-and-Scripts).