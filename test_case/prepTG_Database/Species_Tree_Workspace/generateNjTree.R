library(ape)

input.dist.file <- "/Users/raufs/Coding/zol/develop_v1.6.0/commit/zol/test_case/prepTG_Database/Species_Tree_Workspace/Skani_Based_Distance_Matrix.txt"
output.nwk.file <- "/Users/raufs/Coding/zol/develop_v1.6.0/commit/zol/test_case/prepTG_Database/Species_Tree_Workspace/Unrooted_Species_Tree.nwk"

dat <- read.table(input.dist.file, header=T, sep="\t", row.names=1)
d <- as.dist(as.matrix(dat))
njt <- nj(d)
write.tree(njt, file=output.nwk.file)
