seurat()
library(DoubletFinder)

# read input arguments
args = commandArgs(trailingOnly=T)
seur = readRDS(args[[1]])
sample = args[[2]]
num_pcs = as.integer(args[[3]])
out = args[[4]]

# subset seurat object
cells.use = grep(sample, colnames(seur@data), value=T)
temp = run_seurat(name='temp', dge=seur@raw.data[,cells.use], var_genes=seur@var.genes, num_pcs=num_pcs, max_iter=1, write_out=F)

# run doublet finder
sweep.res = paramSweep_v3(temp, PCs=1:num_pcs, num_pcs=num_pcs)
sweep.stats = summarizeSweep(sweep.res, GT=FALSE)
bcmvn = find.pK(sweep.stats, do.plot=FALSE)
nExp = round(0.1*ncol(temp@data))
pk.use = as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric),'pK']))
print(paste('Optimal pK:', pk.use))
temp = doubletFinder_v3(temp, PCs=1:num_pcs, pN=0.25, pK=pk.use, nExp=nExp)

# save output
saveRDS(temp@data.info, file=out)
