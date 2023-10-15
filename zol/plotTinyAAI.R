library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

info.file <- args[1]
pdf_file <- args[2]

info.dat <- read.table(info.file, header=T, sep='\t')

pdf(pdf_file, height=10, width=10)
ggplot(info.dat, aes(x=AAI, y=Prop_Genes_Found, color=Mean_Syntenic_Correlation)) + geom_point(alpha=0.7) + theme_bw() + scale_color_gradient(low="#e6ffbd", high="#1754b0") + guides(color=guide_legend("Syntenic\nCorrelation\nto\nQuery")) + xlab("Average Amino-Acid Identity") + ylab("Proportion of Query Proteins with Match")
dev.off()
