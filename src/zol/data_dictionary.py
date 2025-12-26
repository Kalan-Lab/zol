import os
import sys

def zol_dd():
    dd = """
Column\tDescription\tNotes
Ortholog Group (OG) ID\tThe identifier of the ortholog/homolog group.\t\t
OG is Single Copy?\tWhether the ortholog/homolog group is single copy in the context of the gene-cluster.\tEvolutionary statistics should be evaluated carefully if False or multiple gene-clusters are from the same genome.
Proportion of Total Gene Clusters with OG\tThe proportion of input gene-clusters/GenBanks which feature the homolog group.\t\t
OG Median Length (bp)\tThe median length of the homolog group in basepairs.\t\t
OG Consensus Order\tThe consensus order of the homolog group across all gene clusters.\t\t
OG Consensus Direction\tThe consensus direction of the homolog group across all gene clusters.\t\t
Custom Annotation (E-value)\tCustom annotation based on user providing custom protein database.\t\t
KO Annotation (E-value)\tBest KEGG ortholog annotation(s) (the HMMER3 E-value associated with the best score)\tHMM-based: HMM profiles are queries, consensus OG sequences are targets.
Pfam Domains\tPfam domains with E-value < 1e-5 and meeting the "trusted" score thresholds.\tHMM-based: HMM profiles are queries, consensus OG sequences are targets.
Proportion of Focal Gene Clusters with OG\t\tOnly produced if comparative analysis is requested by user.
Proportion of Comparator Gene Clusters with OG\t\tOnly produced if comparative analysis is requested by user.
Fixation Index\tFst estimate based on measuring pairwise differences in codon alignments and the statistic developed by Hudson, Slatkin, and Maddison 1992: https://academic.oup.com/genetics/article/132/2/583/6009124?login=true\tOnly produced if comparative analysis is requested by user.
Upstream Region Fixation Index\tFst estimate based on upstream 100 bp nucleotide alignments.\t\t
Tajima's D\tTajima's D calculated using implementation described in lsaBGC based on statistic developed by Tajima 1989: https://academic.oup.com/genetics/article/123/3/585/5998755?login=true\tInterpret with care and consideration of divergence of genomes from which gene clusters were extracted. Calculation of statistic modified to account for the presence of gaps in alignments. Filtering of codon alignments in zol is currently different than what is applied in lsaBGC.
Proportion of Filtered Codon Alignment is Segregating Sites\tProportion of sites in filtered codon alignment which correspond to segregating sites.\tNote, segregating sites require two different non-gap alleles - gaps are not counted as a distinct allele.
Entropy\tAverage entropy over largely non-ambiguous sites (<10% ambiguity) in codon alignments.\t\t
Upstream Region Entropy\tAverage entropy over largely non-ambiguous sites (<10% ambiguity) in nucleotide alignments of upstream regions.\t\t
Median Beta-RD-gc\tThe median Beta-RD statistic for ortholog group relative to the full gene-cluster.\tCalculation is similar to what was described in the lsaBGC study, but expected divergence for ortholog group sequence between pair of gene-clusters/samples is not based on genome-wide divergence but gene-cluster divergence.
Max Beta-RD-gc\tThe max Beta-RD statistic observed for the ortholog group between two gene-clusters.\t\t
Proportion of sites which are highly ambiguous in codon alignment\tThe proportion of sites which are ambiguous (e.g. feature gaps) in greater than 10% of the sequences of a codon alignments (before trimming/filtering).\t\t
Proportion of sites which are highly ambiguous in trimmed codon alignment\tThe proportion of sites which are ambiguous (e.g. feature gaps) in greater than 10% of the sequences of trimmed codon alignments (via trimal).\t\t
Median GC\tThe median GC% of genes belonging to the ortholog group.\t\t
Median GC Skew\tThe median GC skew (G-C)/(G+C) of genes belonging to the ortholog group.\t\t
BGC score (GECCO weights)\tA score from ~ -7 (less BGC like) to ~13 (more BGC like) - based on maximum weight of Pfam domains in GECCO BGC detection (Carroll, Larralde, Fleck, et al. 2021). Primarily relevant for bacterial BGCs. A histogram of values can be found on the wiki documentation link above.\tGECCO paper: https://www.biorxiv.org/content/10.1101/2021.05.03.442509v1.full
Viral score (V-Score)\tA score from 0 (less viral) to 10 (more viral) - based on maximum V-Score (Zhou et al. 2024) for Pfam domains. A histogram of values can be found on the wiki documentation link above.\tV-Score paper: https://www.biorxiv.org/content/10.1101/2024.10.24.619987v1
Hydrophobicity Mean\tThe mean of Hydrophobicity values for proteins belonging to the ortholog group. More info can be found on the [peptides.py documentation](https://peptides.readthedocs.io/en/stable/api/peptide.html#peptides.Peptide.hydrophobicity)\tNote, this metric is computed via the extremely useful [peptides.py](https://github.com/althonos/peptides.py) library. Sequences detected as outliers based on prior distributions from Swiss-Prot are skipped.
Hydrophobicity Std Dev\tThe standard deviation of Hydrophobicity values.\tNote, this metric is computed via the extremely useful [peptides.py](https://github.com/althonos/peptides.py) library. Sequences detected as outliers based on prior distributions from Swiss-Prot are skipped.
Aliphatic Index Mean\tThe mean of Aliphatic Index values for proteins belonging to the ortholog group. Higher Aliphatic Index values coincide with greater thermostability for globular proteins. More info can be found on the [peptides.py documentation](https://peptides.readthedocs.io/en/stable/api/peptide.html#peptides.Peptide.aliphatic_index)\tNote, this metric is computed via the extremely useful [peptides.py](https://github.com/althonos/peptides.py) library. Sequences detected as outliers based on prior distributions from Swiss-Prot are skipped.
Aliphatic Index Std Dev\tThe standard deviation of Aliphatic Index values.\tNote, this metric is computed via the extremely useful [peptides.py](https://github.com/althonos/peptides.py) library. Sequences detected as outliers based on prior distributions from Swiss-Prot are skipped.
m/z Mean\tThe mean of mass over charge (m/z) values for proteins belonging to the ortholog group. More info can be found on the [peptides.py documentation](https://peptides.readthedocs.io/en/stable/api/peptide.html#peptides.Peptide.mz)\tNote, this metric is computed via the extremely useful [peptides.py](https://github.com/althonos/peptides.py) library. Sequences detected as outliers based on prior distributions from Swiss-Prot are skipped.
m/z Std Dev\tThe standard deviation of mean of mass over charge (m/z) values.\tNote, this metric is computed via the extremely useful [peptides.py](https://github.com/althonos/peptides.py) library. Sequences detected as outliers based on prior distributions from Swiss-Prot are skipped.
GARD Partitions Based on Distinct Segments based on Recombination Breakpoints\tNumber of recombination segments detected by HyPhy's GARD method by Kosakovsky Pond et al. 2006: https://academic.oup.com/bioinformatics/article/22/24/3096/208339 .\tNot run by default due to time requirements.
Number of Sites Identified as Under Positive or Negative Selection\tThe number of sites inferred as under positive Prob[α<β] or negative selection Prob[α>β] based on Fubar method: Not run by default due to time requirements. Uses HyPhy's FUBAR method: Murrell et al. 2013: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3670733/\t\t
Average delta(Beta, Alpha) by FUBAR across sites\tThe average difference of β-α across sites in the codon alignment as calculated by FUBAR.\tMore negative values imply greater purifying selection whereas more positive values imply greater positive selection.
Proportion of Sites Under Selection which are Positive\tProportion of the number of sites identified as under either positive or negative selection by FUBAR analysis which are under positive selection.\t\t
P-value for gene-wide episodic selection by BUSTED\tP-value indicating if gene-wide episodic selection has occurred. For more information please checkout the manuscript by Murrell et al. 2015: https://pmc.ncbi.nlm.nih.gov/articles/PMC4408417/.\t\t
PGAP Annotation (E-value)\tBest PGAP annotation(s) (the HMMER3 E-value associated with the best score)\tHMM-based: HMM profiles are queries, consensus OG sequences are targets.
PaperBLAST Annotation (E-value)\tBest PaperBLAST annotation(s) (the DIAMOND E-value associated with the best bitscore). For associated papers BLAST the consensus sequence or the ID here to on the PaperBLAST webpage: https://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi.\tFASTA-based: Consensus OG sequences are queries, database is target.
CARD Annotation (E-value)\tBest CARD annotation(s) of antimicrobial resistance genes (the DIAMOND E-value associated with the best bitscore)\tFASTA-based: Consensus OG sequences are queries, database is target.
TnCentral (E-value)\tBest TnCentral+ISFinder annotation(s) of IS elements / transposon machinery (the DIAMOND E-value associated with the best bitscore)\tFASTA-based: Consensus OG sequences are queries, database is target.
MIBiG Annotation (E-value)\tBest MIBiG annotation(s) for genes in characterized BGCs (the DIAMOND E-value associated with the best bitscore)\tFASTA-based: Consensus OG sequences are queries, database is target.
VOG Annotation (E-value)\tBest VOG annotation(s) for viral/phage ortholog groups (the HMMER3 E-value associated with the best score)\tHMM-based: HMM profiles are queries, consensus OG sequences are targets.
VFDB Annotation (E-value)\tBest VFDB annotation(s) for virulence factor genes (the DIAMOND E-value associated with the best bitscore)\tFASTA-based: Consensus OG sequences are queries, database is target.
CDS Locus Tags\tLocus tag identifiers of genes belonging to the ortholog group.\t\t
Consensus Sequence\tThe consensus sequence for the ortholog group.\t\t
    """
    return(dd)

def fai_dd():
    dd = """
Column\tDescription\tNotes
sample\tThe identifier of the target genome.\t\t
gene-cluster-id\tThe identifier of a discrete neighborhood of genes identified as homologous to the query gene cluster.\tOnly in the "Gene Cluster Instance - Report" tab.
scaffold\tThe scaffold/contig identifier from which the gene cluster was extracted.\tOnly in the "Gene Cluster Instance - Report" tab.
start-coordinate\tThe start coordinate on the scaffold for the extracted gene cluster region.\tOnly in the "Gene Cluster Instance - Report" tab.
end-coordinate\tThe end coordinate on the scaffold for the extracted gene cluster region.\tOnly in the "Gene Cluster Instance - Report" tab.
aggregate-bitscore\tThe aggregate bitscore of hits to the query gene cluster genes.\tOnly the best hit for each query gene/ortholog-group is retained (based on bitscore).
aai-to-query\tThe average amino-acid identity of the proteins in the target genome to the query gene cluster genes.\tOnly the best hit for each query gene/ortholog-group is retained (based on bitscore).
mean-sequence-to-query-ratio\tThe average sequence-to-query ratio of the proteins.\tOnly the best hit for each query gene/ortholog-group is retained (based on bitscore).
proportion-query-genes-found\tThe proportion of query genes/ortholog-groups found in homologous gene clusters across the target genome ("Genome Wide - Report") or in a specific discrete neighborhood ("Gene Cluster Instance - Report")\t\t
avg-syntenic-correlation\tPearson product-moment correlation coefficient for global syntenic similarity of a specific discrete neighborhood to the query gene cluster or the average of these values across all discrete neighborhoods which meet user defined filters.\t\t
number-background-genes\tThe number of background genes in the delineated region within query gene hits which are not represented in the query.\t\t
number-gene-clusters\tThe number of discrete gene neighborhoods which individually meet reporting criteria.\tOnly in the "Genome Wide - Report" tab
copy-counts\tA string separated by commas listing the copy-count of individual query genes/ortholog-groups.\t\t
remaining-columns ...\tThe percent identity & sequence-to-query ratios for best hits of each query gene/ortholog-group. Key proteins defined by the user will have the prefix KEY-PROTEIN.\t\t
    """
    return(dd)
