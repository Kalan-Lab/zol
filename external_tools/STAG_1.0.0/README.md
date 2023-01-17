# STAG: Species Tree inference from All Genes

![STAG logo](stag_logo.png)

STAG is a method for inferring a species tree from a set of gene trees. It can use gene trees inferred from any gene family that contains all the species under consideration. These gene familes can contain paralogues, co-orthologues, or any number of gene duplcation events, and thus STAG is not limited to analysis of gene trees inferred from sets of one-to-one orthologues. This ability to use paralogous gene families enables the use of substantially more data from which to base species tree inference. 

### Additional Advantages of Using STAG
1. STAG trees have realistic branch lengths that can be used for downstream analyses.
2. STAG trees have support values that summarise the fraction of trees that support each bipartition.
3. STAG is fast and memory efficient.

### Installation:
STAG is written in python and can be downloaded from the 'Releases' tab. It is designed to work on Linux operating systems and can be run directly without needing to be installed. STAG requires the FastME program to be installed and in your system path.

#### FastME
FastME can be obtained from http://www.atgc-montpellier.fr/fastme/binaries.php. The package contains a 'binaries/' directory. Choose the appropriate one for your system and copy it to somewhere in the system path e.g. '/usr/local/bin'** and name it 'fastme'. E.g.:

- `sudo cp fastme-2.1.5-linux64 /usr/local/bin/fastme`

### Usage:
An example dataset is supplied, it is the complete set of gene trees for a set of 22 Fungal species and takes less than a minute to run.

To infer the species tree for these species run:
`python stag.py ExampleData/Fungi_SpeciesMap.txt ExampleData/Fungi_Gene_Trees/`
This will analyse the gene trees and write the unrooted species tree to a results directory with todays date. E.g. at the end of the run you should see something like this:

```
Examined 10784 trees
1583 trees had all species present and will be used by STAG to infer the species tree

STAG species tree: /home/david/workspace/git/STAG/stag/ExampleData/STAG_ResultsFeb02/SpeciesTree.tre
```
### Rooting STAG Trees
The tree produced by STAG is unrooted. However, following STAG tree inference the same set of trees used by STAG can be used to root the STAG tree using STRIDE. https://github.com/davidemms/STRIDE: Using the example dataset above the you would root the STAG tree using STRIDE using the following command.

`python ../STRIDE/stride/stride.py -s second_dash -S /home/david/workspace/git/STAG/stag/STAG_ResultsFeb02_1/SpeciesTree.tre -d stag/Fungi_Gene_Trees/`

### Generating Input Datasets for STAG:
The best way to generate a set of gene trees for your species is with OrthoFinder: https://github.com/davidemms/OrthoFinder/

STAG is designed to work with OrthoFinder output. Instead of preparing at species map file you can jsut use the SpeciesIDs.txt file found in the OrthoFinder WorkingDirectory. E.g.
`python stag.py ExampleData/Fungi_SpeciesIDs.txt ExampleData/Fungi_Gene_Trees/`
