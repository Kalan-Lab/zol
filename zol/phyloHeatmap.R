library(ggplot2)
library(ape)
library(dplyr)
library(aplot)

args = commandArgs(trailingOnly=TRUE)

phylo.tree_file <- args[1]
heatmap.data_file <- args[2]
pdf_file <- args[3]
pdf.height <- args[4]
pdf.width <- args[5]

phylo.tree <- read.tree(phylo.tree_file)
heatmap.data <- read.table(heatmap.data_file, header=T, sep='\t')

pdf(pdf_file, height=pdf.height, width=pdf.width)
gg_tr <- ggtree(phylo.tree)
gg_hm <- ggplot(heatmap.data, aes(x=query_prot_id, y=label, fill=bitscore)) +
         theme_classic() + scale_fill_gradient(low='grey', high='black', na.values="white") +
         xlab("Query Proteins/Homolog-groups") + ylab("") + geom_tile(color="white", show.legend=F) +
		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
         theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
gg_hm %>% insert_left(gg_tr, width=0.4)
dev.off()