library(ape)
library(ggtree)

# read input arguments
args = commandArgs(trailingOnly=T)
tree = read.tree(args[[1]])
boot = as.integer(args[[2]])
out = args[[3]]

# collapse nodes with low bootstrap support into polytomies
tree = as.polytomy(tree, feature='node.label', fun=function(x) as.numeric(x) < boot)

# write output tree
write.tree(tree, file=out)
