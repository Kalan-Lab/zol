library(ggplot2)
library(gridExtra)
info.file <- "/Users/raufs/Coding/zol/develop_v1.6.0/commit/zol/test_case/fai_Results_2/Plotting_Files/Tiny_AAI_Plot_Data.txt"
pdf_file <- "/Users/raufs/Coding/zol/develop_v1.6.0/commit/zol/test_case/fai_Results_2/Final_Results/Candidate_Homologous_Gene_Clusters.Tiny_AAI_Plot.pdf"

info.dat <- read.table(info.file, header=T, sep="\t")

pdf(pdf_file, height=10, width=10)
ggplot(info.dat, aes(x=         AAI, y=Prop_Genes_Found, color=Mean_Syntenic_Correlation)) + 
geom_point(alpha=         0.7) + theme_bw() + scale_color_gradient(low="#e6ffbd", high="#1754b0") + 
guides(color=guide_legend("Syntenic Correlation to Query")) + 
xlab("Average amino acid identity") + ylab("Proportion of query proteins with match")
dev.off()
