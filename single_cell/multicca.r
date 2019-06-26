library(Seurat, lib.loc='~/lib/Seurat2')

args = commandArgs(trailingOnly=T)
seur = readRDS(args[[1]])
nccs = as.integer(args[[2]])

seur = UpdateSeuratObject(seur)
seur = ScaleData(seur)
ident = seur@ident
cells1 = colnames(seur@data)

test = sapply(levels(seur@ident), function(a){
    cells.use = colnames(seur@data)[seur@ident == a]
    SubsetData(seur, cells.use=cells.use)
}, simplify=FALSE)
cells2 = sapply(seur, function(a) colnames(a@data))


seur = RunMultiCCA(test, num.ccs=nccs)



new_ident = setNames(ident[colnames(seur@data)], colnames(seur@data))
