library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

input.data_file <- args[1]
pdf_file <- args[2]

dat <- read.table(input.data_file, header=T, sep='\t')

pdf(pdf_file, height=5, width=5)

ggplot(dat, aes(x=ribo_aai, y=gc_aai)) + geom_point(alpha=0.7) + geom_smooth(method='lm', formula= y~x, linetype=2, color='red') +
                                        geom_abline(slope=1, yintercept=0, linetype=2, color="grey") + theme_bw() + 
                                        xlab("Ribosomal Protein AAI") + ylab("Gene Cluster AAI")

dev.off()
