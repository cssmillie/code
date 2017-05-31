library(ape)
library(phangorn)
ifn = commandArgs(trailingOnly=T)[1]
ofn = commandArgs(trailingOnly=T)[2]
t = read.tree(ifn)
T = midpoint(t)
write.tree(T, file=ofn)
