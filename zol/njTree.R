library(ape)

args = commandArgs(trailingOnly=TRUE)

input.dist.file <- args[1]
output.nwk.file <- args[2]

dat <- read.table(input.dist.file, header=T, sep='\t', row.names=1)
d <- as.dist(as.matrix(dat))
print(d)
njt <- nj(d)
write.tree(njt, file=output.nwk.file)
