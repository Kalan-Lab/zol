library(ggplot2)
library(gridExtra)

que.info.file <- "/Users/raufs/Coding/zol/develop_v1.6.10/zol/test_case/fai_Results_3/Plotting_Files/Name_Mapping.txt"
seg.data.file <- "/Users/raufs/Coding/zol/develop_v1.6.10/zol/test_case/fai_Results_3/Plotting_Files/Bar_Plot_Input.txt"
pdf_file <- "/Users/raufs/Coding/zol/develop_v1.6.10/zol/test_case/fai_Results_3/Final_Results/Homologous_Gene_Cluster_Segments.pdf"
pdf_file2 <- "/Users/raufs/Coding/zol/develop_v1.6.10/zol/test_case/fai_Results_3/Final_Results/Name_Mapping_of_Queries_to_Original_LocusTags.pdf"
height <- as.numeric(5)
width <- as.numeric(15)

naming.data <- read.table(que.info.file, header=T, sep="\t")
pdf(pdf_file2, height=30, width=10)
grid.table(naming.data)
dev.off()

segment.data <- read.table(seg.data.file, header=T, sep="\t")
colors <- c("#000000", "#FFFFFF")
names(colors) <- c("True", "False")

samples <- unique(segment.data$sample)
pdf(pdf_file, height=height, width=width)
for (s in samples) {
sample.segment.data <- segment.data[segment.data$sample==s,]
print(unique(sample.segment.data$segment_title))
g<-ggplot(sample.segment.data, aes(x=         reorder(gene, gene_order), y=sql_ratio, fill=identity, color=key)) + 
geom_bar(stat=         "identity", position="dodge", size=1.5) + geom_hline(yintercept=1.0, color="blue", linetype=2) + 
ylab("CDS to Query Length Ratio") + facet_grid(.~segment_title, space=         "free", scales="free",  labeller = label_wrap_gen(width = 50, multi_line = TRUE)) + 
xlab("CDS Classifications") + theme_bw() + theme(legend.position=         "bottom", axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) + 
scale_color_manual(values=         colors) + scale_fill_gradient(limits=c(0.0, 100.0), breaks=c(0.0, 50.0, 100.0), low="grey", high="red")
print(g)
}
dev.off()

