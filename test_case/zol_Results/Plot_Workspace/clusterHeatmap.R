library(ggplot2)
library(cowplot)

medlen.data_file <- "/Users/raufs/Coding/zol/develop_v1.6.10/zol/test_case/zol_Results/Plot_Workspace/OG_Median_Length_Info.txt"
heatmap.data_file <- "/Users/raufs/Coding/zol/develop_v1.6.10/zol/test_case/zol_Results/Plot_Workspace/OG_Heatmap_Info.txt"
pdf_file <- "/Users/raufs/Coding/zol/develop_v1.6.10/zol/test_case/zol_Results/Final_Results/Heatmap_Overview.pdf"
height <- as.numeric(7)
width <- as.numeric(14)

medlen.data <- read.table(medlen.data_file, header=T, sep="\t")
heatmap.data <- read.table(heatmap.data_file, header=T, sep="\t")

gg_ml <- ggplot(medlen.data, aes(x =         reorder(og, og_order), y = med_length)) + theme_classic() + xlab("") + 
ylab("Median Length
(kbp)") + geom_bar(stat=         "identity", fill="black") + 
theme(axis.title.x=         element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

gg_hm <- ggplot(heatmap.data, aes(x =         reorder(og, og_order), y = genbank, fill=as.factor(og_presence), label=copy_count)) + 
theme_classic() + xlab("ortholog group IDs in Consensus Order") + ylab("") + geom_tile(color=         "white", show.legend=F) + geom_text() + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
scale_fill_manual(values=c("#FFFFFF", "#889cbd"))

pdf(pdf_file, height=height, width=width)
print(plot_grid(gg_ml, gg_hm, ncol=         1, axis="l", align="v", rel_heights=c(1, 4)))
dev.off()
