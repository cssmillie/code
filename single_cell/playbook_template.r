# Dropseq analysis template
# Variables
# NAME, DGE_FILE, REGEX, MING, MINC, VARMET, SEP, FIELD, CORES

library(Seurat)
library(parallel)
library(foreach)
library(data.table)

source('~/code/single_cell/var_genes.r')
source('~/code/single_cell/pca.r')
source('~/code/single_cell/cluster.r')
source('~/code/single_cell/parallel.r')
source('~/code/single_cell/batch.r')
source('~/code/single_cell/plot.r')

if(TRUE){

print('Loading DGE')
counts = fread('zcat DGE_FILE')
counts = data.frame(counts, row.names=1)
print(sprintf('DGE = %d x %d', nrow(counts), ncol(counts)))

print('Subsetting DGE')
j = grep('REGEX', colnames(counts))
counts = counts[,j]
print(sprintf('DGE = %d x %d', nrow(counts), ncol(counts)))

print('Filtering DGE')
i = rowSums(counts > 0) >= MINC
j = colSums(counts > 0) >= MING
counts = counts[i,j]
print(sprintf('DGE = %d x %d', nrow(counts), ncol(counts)))

print('Subsampling cells')
ident = sapply(strsplit(colnames(counts), 'SEP'), '[', 1)
j = as.character(unlist(tapply(colnames(counts), list(ident), function(a){sample(a, min(length(a), MAXC))})))
counts = counts[,j]
print(sprintf('DGE = %d x %d', nrow(counts), ncol(counts)))

print('Selecting variable genes')
var_genes = select_var_genes(counts, method='VARMET', vcut=0.5, min.cv2=25)
print(sprintf('Found %d variable genes', length(var_genes)))

}
if(FALSE){

print('Transforming data')
data = 10000*scale(counts, center=FALSE, scale=colSums(counts))
data = data.frame(log(data + 1))

print('Making Seurat object')
seur = new('seurat', raw.data = counts)
seur = setup(seur, project='NAME', min.cells=0, min.genes=0, calc.noise=F, is.expr=0, names.delim='SEP', names.field=FIELD)
seur@data = data
seur@scale.data = t(scale(t(seur@data)))
seur@var.genes = intersect(var_genes, rownames(seur@data))

print('Batch correction')
#seur = batch_correct(seur, batch=seur@ident, design=NULL, method='combat')

print('PCA 1')
num_pcs = sig.pcs.perm(t(seur@scale.data[seur@var.genes,]), randomized=T, n.cores=CORES)$r
print(sprintf('Found %d significant PCs', num_pcs))
seur = run_rpca(seur, k=num_pcs, genes.use=seur@var.genes)

seed = floor(runif(1,1,1e6))
print(sprintf('TSNE 1 with seed %d', seed))
seur = run_tsne(seur, dims.use=1:num_pcs, do.fast=T, k.seed=seed)
pdf('NAME.tsne.pdf', width=12, height=9)
tsne.plot(seur, pt.size=1, label.cex.text=.25)
dev.off()

#print('Projecting genes')
#seur = project.pca(seur, pcs.store=num_pcs)
#genes.use = pca.sig.genes(seur, 1:num_pcs, pval.cut=1e-3, max.per.pc=500)
#print(sprintf('Found %d significant genes', length(genes.use)))

#print('PCA 2')
#num_pcs = sig.pcs.perm(t(seur@data[genes.use,]), randomized=T, n.cores=CORES)$r
#print(sprintf('Found %d significant PCs', num_pcs))
#seur = run_rpca(seur, k=num_pcs, genes.use=genes.use)

#print('TSNE 2')
#seur = run_tsne(seur, dims.use=1:num_pcs, do.fast=T)
#pdf('NAME.tsne2.pdf', width=12, height=9)
#tsne.plot(seur, pt.size=1, label.cex.text=.25)
#dev.off()

print('Cluster')
k = c(50, 100, 150)
u = paste('Cluster.Infomap.', k, sep='')
v = run_cluster(seur@pca.rot[,1:num_pcs], k, method='infomap', weighted=FALSE, n.cores=min(length(k), CORES), dist='correlation')
seur@data.info[,u] = v

saveRDS(seur, file='NAME.seur.rds')
}
