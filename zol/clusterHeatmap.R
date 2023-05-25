library(ggplot2)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)

medlen.data_file <- args[1]
heatmap.data_file <- args[2]
pdf_file <- args[3]
height <- as.numeric(args[4])
width <- as.numeric(args[5])

medlen.data <- read.table(medlen.data_file, header=T, sep='\t')
heatmap.data <- read.table(heatmap.data_file, header=T, sep='\t')

gg_ml <- ggplot(medlen.data, aes(x = reorder(og, og_order), y = med_length)) + theme_classic() + xlab("") +
	ylab("Median Length\n(kbp)") + geom_bar(stat='identity', fill='black') +
	theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
gg_hm <- ggplot(heatmap.data, aes(x = reorder(og, og_order), y = genbank, fill=as.factor(og_presence), label=copy_count)) +
         theme_classic() +
         xlab("ortholog group IDs in Consensus Order") + ylab("") + geom_tile(color='white', show.legend=F) + geom_text() +
		 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	     scale_fill_manual(values=c('#FFFFFF', '#889cbd'))

pdf(pdf_file, height=height, width=width)
print(plot_grid(gg_ml, gg_hm, ncol=1, axis="l", align='v', rel_heights=c(1, 4)))
dev.off()