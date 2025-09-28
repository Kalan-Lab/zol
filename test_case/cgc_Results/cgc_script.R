library(ggplot2)
library(cowplot)
library(gggenes)

dat <- read.table("/Users/raufs/Coding/zol/develop_v1.6.10/zol/test_case/cgc_Results/Plot_Input.txt", header = T, sep = "\t")

# track for tajimas_d
g1 <- ggplot(dat, aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = tajimas_d, fill = tajimas_d)) + theme_classic() + 
geom_rect(color = "black", show.legend = F) + scale_fill_gradient2(low = "#bd3131", mid = "#ffffff", high = "#548ce8", na.value = "grey50") +
theme(axis.line.x = element_line(colour = "white"), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
ggtitle("Tajima's D")

# track for conservation
g2 <- ggplot(dat, aes(xmin = x_start, xmax = x_end, y = 1, label = label, forward = direction)) +
geom_gene_arrow(aes(fill = conservation), color = "black") + theme_void() + scale_fill_gradient2(low = "#ffffff", mid = "#878787", high = "#000000", na.value = "grey50", guide = "colourbar", aesthetics = "fill") + 
geom_text(aes(x = x_midpoint), angle = 45, y = 0.975, size = 2.0) +
theme(legend.position = "bottom", legend.title = element_text(size = 10), legend.box = "horizontal")


pdf("/Users/raufs/Coding/zol/develop_v1.6.10/zol/test_case/cgc_Results/cgc_plot.pdf", height = 7, width = 7)
print(plot_grid(g1, NULL, g2, NULL, NULL, ncol = 1, axis = "lr", align = "v", rel_heights = c(1, 0.0, 1, 0.1)))
dev.off()
