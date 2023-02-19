library(ggplot2)
library(gggenes)

args = commandArgs(trailingOnly=TRUE)
input.file <- args[1]
height <- as.numeric(args[2])
width <- as.numeric(args[3])
pdf.file <- args[4]

input.data <- read.table(file=input.file, sep='\t', header=T)
# 	pif_handle.write('\t'.join(['HG', 'Start', 'End', 'Direction', 'SC', 'Metric']) + '\n')

pdf(pdf.file, height=height, width=width)
ggplot(input.data, aes(xmin=Start, xmax = End, y = "", forward = Direction, label=SC)) +
  geom_gene_arrow(aes(fill=Metric)) + theme_classic() +
  scale_fill_gradient2(low='#e05c74', mid="#f2f2f2", high='#2087b3', breaks =c(-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0),
                      labels=c("", "-2", "", "0", "", "2", ""), limits=c(-3,3), na.value='grey50', guide='colourbar',
                      aesthetics='fill') +
  geom_gene_label(align='centre', min.size=5) + theme(legend.position="bottom")
dev.off()