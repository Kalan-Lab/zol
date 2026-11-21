import os
import sys

def zol_dd():
    dd = """
Column\tDescription\tNotes
Ortholog Group (OG) ID\tThe identifier of the ortholog/homolog group.\t\t
OG is Single Copy?\tWhether the ortholog/homolog group is single copy in the context of the gene-cluster.\tEvolutionary statistics should be evaluated carefully if False or multiple gene-clusters are from the same genome.
Proportion of Total Gene Clusters with OG\tThe proportion of input gene-clusters/GenBanks which feature the homolog group.\t\t
OG Median Length (bp)\tThe median length of the homolog group in basepairs.\t\t
Single-Linkage Full Protein Cluster\tIdentifier of the cluster of full-length proteins that the domain ortholog group belongs to, determined by single-linkage clustering of full proteins which share domain ortholog groups.\tOnly produced in domain mode.
OG Consensus Order\tThe consensus order of the homolog group across all gene clusters.\t\t
OG Consensus Direction\tThe consensus direction of the homolog group across all gene clusters.\t\t
Custom Annotation (E-value)\tCustom annotation based on user providing custom protein database.\t\t
KO Annotation (E-value)\tBest KEGG ortholog annotation(s) (the HMMER3 E-value associated with the best score)\tHMM-based: HMM profiles are queries, consensus OG sequences are targets.
Pfam Domains\tPfam domains meeting criteria of `-pfp` parameter, listed as "Pfam accession: description" and separated by semicolons. By default, the criteria are E-value < 10 and meeting the "gathering" score thresholds.\tHMM-based: HMM profiles are queries, consensus OG sequences are targets.
Proportion of Focal Gene Clusters with OG\t\tOnly produced if comparative analysis is requested by user.
Proportion of Comparator Gene Clusters with OG\t\tOnly produced if comparative analysis is requested by user.
Fixation Index\tFst estimate based on measuring pairwise differences in codon alignments and the statistic developed by Hudson, Slatkin, and Maddison 1992: https://academic.oup.com/genetics/article/132/2/583/6009124?login=true\tOnly produced if comparative analysis is requested by user.
Upstream Region Fixation Index\tFst estimate based on upstream 100 bp nucleotide alignments.\t\t
McDonald-Kreitman Dn\tNumber of non-synonymous fixed differences between focal and comparator gene clusters for the McDonald-Kreitman test (McDonald and Kreitman 1991: https://www.nature.com/articles/351652a0). A fixed difference is a nucleotide position within a codon that is monomorphic in both sets but differs between them.\tOnly produced if comparative analysis is requested by user. Computed on untrimmed, in-frame codon alignments (codon sites with gaps/ambiguity in >10% of sequences are skipped). Changes in codons differing at multiple positions are averaged over mutational pathways (Nei and Gojobori 1986), so counts can be fractional.
McDonald-Kreitman Ds\tNumber of synonymous fixed differences between focal and comparator gene clusters.\tOnly produced if comparative analysis is requested by user. Computed on untrimmed, in-frame codon alignments (codon sites with gaps/ambiguity in >10% of sequences are skipped).
McDonald-Kreitman Pn\tNumber of non-synonymous polymorphisms segregating within the focal and/or comparator gene clusters, counted along the minimum spanning network of distinct codons observed at each codon site.\tOnly produced if comparative analysis is requested by user. Computed on untrimmed, in-frame codon alignments (codon sites with gaps/ambiguity in >10% of sequences are skipped).
McDonald-Kreitman Ps\tNumber of synonymous polymorphisms segregating within the focal and/or comparator gene clusters.\tOnly produced if comparative analysis is requested by user. Computed on untrimmed, in-frame codon alignments (codon sites with gaps/ambiguity in >10% of sequences are skipped).
McDonald-Kreitman Neutrality Index\tNeutrality index, NI = (Pn/Ps)/(Dn/Ds) (Rand and Kann 1996). Values < 1 suggest an excess of non-synonymous divergence (adaptive evolution / positive selection) and values > 1 suggest an excess of non-synonymous polymorphism (e.g. slightly deleterious variants under purifying selection). The proportion of adaptive substitutions can be estimated as alpha = 1 - NI.\tNA if Ps or Dn is zero.
McDonald-Kreitman Direction of Selection\tDirection of selection, DoS = Dn/(Dn+Ds) - Pn/(Pn+Ps) (Stoletzki and Eyre-Walker 2011). Positive values suggest adaptive evolution and negative values suggest segregating slightly deleterious variation. Less biased than the NI when counts are small.\tNA if there are no fixed differences or no polymorphisms.
McDonald-Kreitman Fisher's Exact P-value\tTwo-sided Fisher's exact test P-value for the 2x2 McDonald-Kreitman contingency table (Dn, Ds vs. Pn, Ps). Counts are rounded to the nearest integer for the test.\tNot corrected for multiple testing. Interpret with care when focal and comparator sets are not well-separated lineages (e.g. recombination or horizontal transfer between sets).
Average Pairwise Nucleotide Identity\tThe average pairwise nucleotide identity computed over all unique sequence pairs in the trimmed codon alignment. Identity for each pair is calculated as the proportion of non-gap positions that are identical.\tOnly computed for single-copy ortholog groups by default. Maximum of 5000 pairwise comparisons performed.
Average Pairwise Protein Identity\tThe average pairwise protein identity computed over all unique sequence pairs in the protein alignment. Identity for each pair is calculated as the proportion of non-gap positions that are identical.\tOnly computed for single-copy ortholog groups by default. Maximum of 5000 pairwise comparisons performed.
Std Dev of Pairwise Protein Identity\tThe standard deviation of pairwise protein identity values across all unique sequence pairs in the protein alignment.\tOnly computed for single-copy ortholog groups by default.
Protein to Nucleotide Identity Ratio\tThe average pairwise protein identity divided by the average pairwise nucleotide identity for the ortholog group.\tOnly computed for single-copy ortholog groups by default.
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
Viral score (V-Score)\tA score from 0 (less viral) to 10 (more viral) - based on maximum V-Score (Zhou et al. 2026) for Pfam domains. A histogram of values can be found on the wiki documentation link above.\tV-Score paper: https://www.nature.com/articles/s41467-026-72028-0
Hydrophobicity Mean\tThe mean of Hydrophobicity values for proteins belonging to the ortholog group. More info can be found on the [peptides.py documentation](https://peptides.readthedocs.io/en/stable/api/peptide.html#peptides.Peptide.hydrophobicity)\tNote, this metric is computed via the extremely useful [peptides.py](https://github.com/althonos/peptides.py) library. Sequences detected as outliers based on prior distributions from Swiss-Prot are skipped.
Hydrophobicity Std Dev\tThe standard deviation of Hydrophobicity values.\tNote, this metric is computed via the extremely useful [peptides.py](https://github.com/althonos/peptides.py) library. Sequences detected as outliers based on prior distributions from Swiss-Prot are skipped.
Aliphatic Index Mean\tThe mean of Aliphatic Index values for proteins belonging to the ortholog group. Higher Aliphatic Index values coincide with greater thermostability for globular proteins. More info can be found on the [peptides.py documentation](https://peptides.readthedocs.io/en/stable/api/peptide.html#peptides.Peptide.aliphatic_index)\tNote, this metric is computed via the extremely useful [peptides.py](https://github.com/althonos/peptides.py) library. Sequences detected as outliers based on prior distributions from Swiss-Prot are skipped.
Aliphatic Index Std Dev\tThe standard deviation of Aliphatic Index values.\tNote, this metric is computed via the extremely useful [peptides.py](https://github.com/althonos/peptides.py) library. Sequences detected as outliers based on prior distributions from Swiss-Prot are skipped.
m/z Mean\tThe mean of mass over charge (m/z) values for proteins belonging to the ortholog group. More info can be found on the [peptides.py documentation](https://peptides.readthedocs.io/en/stable/api/peptide.html#peptides.Peptide.mz)\tNote, this metric is computed via the extremely useful [peptides.py](https://github.com/althonos/peptides.py) library. Sequences detected as outliers based on prior distributions from Swiss-Prot are skipped.
m/z Std Dev\tThe standard deviation of mean of mass over charge (m/z) values.\tNote, this metric is computed via the extremely useful [peptides.py](https://github.com/althonos/peptides.py) library. Sequences detected as outliers based on prior distributions from Swiss-Prot are skipped.
GARD Partitions Based on Recombination Breakpoints\tNumber of recombination segments detected by HyPhy's GARD method by Kosakovsky Pond et al. 2006: https://academic.oup.com/bioinformatics/article/22/24/3096/208339 .\tNot run by default due to time requirements.
Number of Sites Identified as Under Positive or Negative Selection by FUBAR\tThe number of sites inferred as under positive Prob[α<β] or negative selection Prob[α>β] based on Fubar method: Not run by default due to time requirements. Uses HyPhy's FUBAR method: Murrell et al. 2013: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3670733/\t\t
Average delta(Beta, Alpha) by FUBAR across sites\tThe average difference of β-α across sites in the codon alignment as calculated by FUBAR.\tMore negative values imply greater purifying selection whereas more positive values imply greater positive selection.
Gene-wide dN/dS by MG94\tAlignment-wide non-synonymous to synonymous substitution rate ratio (omega) estimated by fitting the MG94xREV codon model (Muse and Gaut 1994: https://academic.oup.com/mbe/article/11/5/715/1008710) with a single omega across sites and branches using HyPhy's fitmodel analysis and the ortholog group's gene tree (or GARD recombination-aware partitions if available).\tOnly computed if selection analysis is requested. Values < 1 suggest purifying selection and > 1 positive selection; gene-wide estimates are typically < 1 even when a few sites are under positive selection. Interpret with care for very closely related sequences (omega is biased towards 1) or highly divergent sequences (synonymous saturation).
Gene-wide dN/dS by FUBAR (Sum of Beta / Sum of Alpha)\tGene-wide dN/dS computed as the ratio of the summed site-level posterior mean non-synonymous (β) and synonymous (α) rates inferred by FUBAR.\tOnly computed if selection analysis is requested. Unlike the MG94 estimate, this allows synonymous rates to vary across sites. Site-level posterior means are shrunk towards FUBAR's prior when sites carry little information (e.g. few or highly similar sequences), which inflates this ratio towards 1 relative to the MG94 estimate - prefer the MG94 estimate for gene-wide comparisons and use this as a complementary, site-heterogeneity-aware view.
Proportion of Sites Under Selection which are Positive\tProportion of the number of sites identified as under either positive or negative selection by FUBAR analysis which are under positive selection.\t\t
P-value for gene-wide episodic selection by BUSTED\tP-value indicating if gene-wide episodic selection has occurred. For more information please checkout the manuscript by Murrell et al. 2015: https://pmc.ncbi.nlm.nih.gov/articles/PMC4408417/.\tWarning: BUSTED p-values can vary substantially between runs on identical input because of stochasticity in its likelihood optimization (e.g. 0.035 to 0.455 across five runs of the same ortholog group in zol's test case). Treat values near significance thresholds with caution and confirm with repeated runs before drawing conclusions. Not corrected for multiple testing.
PGAP Annotation (E-value)\tBest PGAP annotation(s) (the HMMER3 E-value associated with the best score)\tHMM-based: HMM profiles are queries, consensus OG sequences are targets.
PaperBLAST Annotation (E-value)\tBest PaperBLAST annotation(s) (the DIAMOND E-value associated with the best bitscore). For associated papers BLAST the consensus sequence or the ID here to on the PaperBLAST webpage: https://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi.\tFASTA-based: Consensus OG sequences are queries, database is target.
CARD Annotation (E-value)\tBest CARD annotation(s) of antimicrobial resistance genes (the DIAMOND E-value associated with the best bitscore)\tFASTA-based: Consensus OG sequences are queries, database is target.
IS Finder (E-value)\tBest ISFinder annotation(s) of IS elements / transposons (the DIAMOND E-value associated with the best bitscore)\tFASTA-based: Consensus OG sequences are queries, database is target.
MIBiG Annotation (E-value)\tBest MIBiG annotation(s) for genes in characterized BGCs (the DIAMOND E-value associated with the best bitscore)\tFASTA-based: Consensus OG sequences are queries, database is target.
VOG Annotation (E-value)\tBest VOG annotation(s) for viral/phage ortholog groups (the HMMER3 E-value associated with the best score)\tHMM-based: HMM profiles are queries, consensus OG sequences are targets.
VFDB Annotation (E-value)\tBest VFDB annotation(s) for virulence factor genes (the DIAMOND E-value associated with the best bitscore)\tFASTA-based: Consensus OG sequences are queries, database is target.
DefenseFinder Annotation (E-value)\tBest DefenseFinder annotation(s) of anti-phage defense system genes or anti-defense factors, listed as "gene (system)" (the HMMER3 E-value associated with the best score). Only hits meeting the gathering (GA) score thresholds of the DefenseFinder models are considered. Note, this is a gene-level annotation and does not assess whether a complete defense system is present - see DefenseFinder for system-level calls.\tHMM-based: HMM profiles are queries, consensus OG sequences are targets. Models from https://github.com/mdmparis/defense-finder-models (Tesson et al. 2022; Tesson et al. 2024).
CDS Locus Tags\tLocus tag identifiers of genes belonging to the ortholog group.\t\t
OG Consensus Sequence\tThe consensus sequence for the ortholog group.\t\t
    """
    return(dd)

def zol_scape_dd():
    """
    Data dictionary for the zol-scape results spreadsheet, which consolidates zol reports for each BiG-SCAPE GCF and
    thus shares zol's columns, plus two columns identifying the clustering file and GCF.
    """
    zol_lines = zol_dd().strip("\n").split("\n")
    header_line = zol_lines[0]
    zol_rows = zol_lines[1:]
    scape_rows = [
        "Clustering File Path\tPath to the BiG-SCAPE clustering file (*_clustering_*.tsv) in which the GCF was defined.\t\t",
        "GCF ID\tIdentifier of the gene cluster family (GCF) as defined by BiG-SCAPE in the clustering file.\tzol was run separately for each GCF, so ortholog group identifiers are only meaningful within a GCF.",
    ]
    return "\n" + "\n".join([header_line] + scape_rows + zol_rows) + "\n"

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
