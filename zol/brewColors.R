# Title     : TODO
# Objective : TODO
# Created by: rauf
# Created on: 7/24/21

library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

num_colors <- as.numeric(args[1])
output_file <- args[2]

mycolors <- colorRampPalette(brewer.pal(12, "Spectral"))(num_colors)
write(mycolors, file=output_file)