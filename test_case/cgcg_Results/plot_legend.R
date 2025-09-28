library(ggplot2)
library(grid)
library(cowplot)

png("/Users/raufs/Coding/zol/develop_v1.6.10/zol/test_case/cgcg_Results/Final_Results/plot_legend.png", height = 10, width = 5, units = "in", res = 600)
V1 <- c(0.0, 1.0)
V2 <- c("low_value", "high_value")
dat <- data.frame(V1, V2)
my_hist <- ggplot(dat, aes(x = V1, y = 1, fill = V1)) + geom_bar(stat = "identity") + theme(legend.title = element_blank()) +
scale_fill_gradient(low = "#a2b7e8", high = "#102f75", limits = c(0.0, 1.0))
legend <- cowplot::get_legend(my_hist)
grid.newpage()
grid.draw(legend)
dev.off()

