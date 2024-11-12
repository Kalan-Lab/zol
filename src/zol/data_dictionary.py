import os
import sys

def zol_dd():
    dd = """
Column\tDescription\tNotes
Ortholog Group (OG) ID\tThe identifier of the ortholog/homolog group.\t
OG is Single Copy?\tWhether the ortholog/homolog group is single copy in the context of the gene-cluster.\tEvolutionary statistics should be evaluated carefully if False or multiple gene-clusters are from the same genome.
Proportion of Total Gene Clusters with OG\tThe proportion of input gene-clusters/GenBanks which feature the homolog group.\t
OG Median Length (bp)\tThe median length of the homolog group in basepairs.\t
OG Consensus Order\tThe consensus order of the homolog group across all gene clusters.\t
OG Consensus Direction\tThe consensus direction of the homolog group across all gene clusters.\t
Proportion of Focal Gene Clusters with OG\t\tOnly produced if comparative analysis is requested by user.
Proportion of Comparator Gene Clusters with OG\t\tOnly produced if comparative analysis is requested by user.
Fixation Index\tFst estimate based on measuring pairwise differences in codon alignments and the statistic developed by Hudson, Slatkin, and Maddison 1992: https://academic.oup.com/genetics/article/132/2/583/6009124?login=true\tOnly produced if comparative analysis is requested by user.
Upstream Region Fixation Index\tFst estimate based on upstream 100 bp nucleotide alignments.\t
Tajima's D\tTajima's D calculated using implementation described in lsaBGC based on statistic developed by Tajima 1989: https://academic.oup.com/genetics/article/123/3/585/5998755?login=true\tInterpret with care and consideration of divergence of genomes from which gene clusters were extracted. Calculation of statistic modified to account for the presence of gaps in alignments. Filtering of codon alignments in zol is currently different than what is applied in lsaBGC.
Proportion of Filtered Codon Alignment is Segregating Sites\tProportion of sites in filtered codon alignment which correspond to segregating sites.\tNote, segregating sites require two different non-gap alleles - gaps are not counted as a distinct allele.
Entropy\tAverage entropy over largely non-ambiguous sites (<10\% ambiguity) in codon alignments.\t
Upstream Region Entropy\tAverage entropy over largely non-ambiguous sites (<10\% ambiguity) in nucleotide alignments of upstream regions.\t
Median Beta-RD-gc\tThe median Beta-RD statistic for ortholog group relative to the full gene-cluster.\tCalculation is similar to what was described in the lsaBGC study, but expected divergence for ortholog group sequence between pair of gene-clusters/samples is not based on genome-wide divergence but gene-cluster divergence.
Max Beta-RD-gc\tThe max Beta-RD statistic observed for the ortholog group between two gene-clusters.\t
Proportion of sites which are highly ambiguous in codon alignment\tThe proportion of sites which are ambiguous (e.g. feature gaps) in greater than 10\% of the sequences of a codon alignments (before trimming/filtering).\t
Proportion of sites which are highly ambiguous in trimmed codon alignment\tThe proportion of sites which are ambiguous (e.g. feature gaps) in greater than 10\% of the sequences of trimmed codon alignments (via trimal).\t
Median GC\tThe median GC\% of genes belonging to the ortholog group.\t
Median GC Skew\tThe median GC skew (G-C)/(G+C) of genes belonging to the ortholog group.\t
GARD Partitions Based on Distinct Segments based on Recombination Breakpoints\tNumber of recombination segments detected by HyPhy's GARD method by Kosakovsky Pond et al. 2006: https://academic.oup.com/bioinformatics/article/22/24/3096/208339 .\tNot run by default due to time requirements.
Number of Sites Identified as Under Positive or Negative Selection\tThe number of sites inferred as under positive Prob[α<β] or negative selection Prob[α>β] based on FUBAR method: Not run by default due to time requirements. Uses HyPhy's FUBAR method: Murrell et al. 2013: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3670733/\t
Average delta(Beta, Alpha) by FUBAR across sites\tThe average difference of β-α across sites in the codon alignment as calculated by FUBAR.\tMore negative values imply greater purifying selection whereas more positive values imply greater positive selection.
Proportion of Sites Under Selection which are Positive\tProportion of the number of sites identified as under either positive or negative selection by FUBAR analysis which are under positive selection.\t
Custom Annotation (E-value)\tCustom annotation based on user providing custom protein database.\t
KO Annotation (E-value)\tBest KEGG ortholog annotation(s) (the HMMER3 E-value associated with the best score)\t
PGAP Annotation (E-value)\tBest PGAP annotation(s) (the HMMER3 E-value associated with the best score)\t
PaperBLAST Annotation (E-value)\tBest PaperBLAST annotation(s) (the DIAMOND E-value associated with the best bitscore). For associated papers BLAST the consensus sequence or the ID here to on the PaperBLAST webpage: https://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi.\t
CARD Annotation (E-value)\tBest CARD annotation(s) of antimicrobial resistance genes (the DIAMOND E-value associated with the best bitscore)\t
IS Finder (E-value)\tBest ISFinder annotation(s) of IS elements / transposons (the DIAMOND E-value associated with the best bitscore)\t
MIBiG Annotation (E-value)\tBest MIBiG annotation(s) for genes in characterized BGCs (the DIAMOND E-value associated with the best bitscore)\t
VOG Annotation (E-value)\tBest VOG annotation(s) for viral/phage ortholog groups (the HMMER3 E-value associated with the best score)\t
Pfam Domains\tPfam domains with E-value < 1e-5 and meeting the "trusted" score thresholds.\t
CDS Locus Tags\tLocus tag identifiers of genes belonging to the ortholog group.\t
Consensus Sequence\tThe consensus sequence for the ortholog group.\t
    """
    return(dd)

def fai_dd():
    dd = """
Column\tDescription\tNotes
sample\tThe identifier of the target genome.\t
gene-cluster-id\tThe identifier of a discrete neighborhood of genes identified as homologous to the query gene cluster.\tOnly in the "Gene Cluster Instance - Report" tab.
aggregate-bitscore\tThe aggregate bitscore of hits to the query gene cluster genes.\tOnly the best hit for each query gene/ortholog-group is retained (based on bitscore).
aai-to-query\tThe average amino-acid identity of the proteins in the target genome to the query gene cluster genes.\tOnly the best hit for each query gene/ortholog-group is retained (based on bitscore).
mean-sequence-to-query-ratio\tThe average sequence-to-query ratio of the proteins.\tOnly the best hit for each query gene/ortholog-group is retained (based on bitscore).
proportion-query-genes-found\tThe proportion of query genes/ortholog-groups found in homologous gene clusters across the target genome ("Genome Wide - Report") or in a specific discrete neighborhood ("Gene Cluster Instance - Report")\t
avg-syntenic-correlation\tPearson product-moment correlation coefficient for global syntenic similarity of a specific discrete neighborhood to the query gene cluster or the average of these values across all discrete neighborhoods which meet user defined filters.\t
number-background-genes\tThe number of background genes in the delineated region within query gene hits which are not represented in the query.\t
number-gene-clusters\tThe number of discrete gene neighborhoods which individually meet reporting criteria.\tOnly in the "Genome Wide - Report" tab
copy-counts\tA string separated by commas listing the copy-count of individual query genes/ortholog-groups.\t
remaining-columns ...\tThe percent identity & sequence-to-query ratios for best hits of each query gene/ortholog-group. Key proteins defined by the user will have the prefix KEY-PROTEIN.\t
    """
    return(dd)
