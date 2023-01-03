library(ggplot2)
library(ggtree)
library(gggenes)
library(ape)
library(phytools)
library(aplot)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

phylo.tree_file <- args[1]
genes.data_file <- args[2]
heatmap.data_file <- args[3]
pdf_file <- args[4]

phylo.tree <- read.tree(phylo.tree_file)
phylo.tree <- midpoint.root(phylo.tree)
genes.data <- read.table(genes.data_file, header=T, sep='\t')
heatmap.data <- read.table(heatmap.data_file, header=T, sep='\t')

tree.labels <- phylo.tree$tip.label

genes.sub.data <- distinct(genes.data[c("og", "og_color")])
og_colors <- na.omit(c(genes.sub.data$og_color))
names(og_colors) <- c(na.omit(c(genes.sub.data$og)))

print(og_colors)

pdf(pdf_file, height=30, width=30)
gg_tr <- ggtree(phylo.tree)#  + ggplot2::xlim(NA, 1)
gg_gn <- ggplot(genes.data, aes(xmin = start, xmax = end, y = label, fill = og, forward = forward)) +
          geom_gene_arrow(show.legend=F) + theme_void() + scale_fill_manual(values=og_colors)
gg_hm <- ggplot(heatmap.data, aes(x = reorder(og, -og_count), y = label, fill=og_presence)) +
         theme_classic() + scale_fill_manual(values=og_colors) +
         xlab("Homolog Group IDs") + ylab("BGC IDs") + geom_tile(color='white', show.legend=F) +
		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gg_hm %>% insert_left(gg_tr, width=0.4) %>% insert_right(gg_gn, width=1.0)
dev.off()
