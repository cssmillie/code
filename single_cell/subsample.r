# This code subsamples N cells per cell type
# It then does PCA and TSNE on the cell subsets
# Then it projects all cells onto these axes

msg(name, 'Subsampling cells for PCA', verbose)
ident = seur@data.info$orig.ident
cells.use = as.character(unlist(tapply(colnames(seur@data), list(ident), function(a){sample(a, min(length(a), maxc))})))
msg(name, sprintf('Selected %d total cells', length(cells.use)))

msg(name, 'Subsampled PCA', verbose)
if(num_pcs == 0){num_pcs = sig.pcs.perm(t(seur@scale.data[seur@var.genes,cells.use]), randomized=T, n.cores=ncores)$r}
seur@data.info$num_pcs = num_pcs
msg(name, sprintf('Found %d significant PCs', num_pcs), verbose)
seur = run_rpca(seur, k=25, genes.use=seur@var.genes, cells.use=cells.use)

msg(name, sprintf('Subsampled TSNE', seed), verbose)
tsne.rot = Rtsne(seur@pca.rot[cells.use,1:num_pcs], do.fast=T, max_iter=max_iter, verbose=T, perplexity=perplexity)@Y[,1:2]

msg(name, 'Projecting cells', verbose)
seur@pca.rot = project_pca(seur@pca.obj, seur@scale.data[var_genes,])
new_cells = setdiff(colnames(seur@data), cells.use)
seur@tsne.rot = data.frame(matrix(NA, nrow=ncol(seur@data), 2), row.names=colnames(seur@data))
seur@tsne.rot[cells.use,] = tsne.rot$Y
seur@tsne.rot[new_cells,] = project_tsne(seur@pca.rot[new_cells,1:num_pcs], seur@pca.rot[cells.use,1:num_pcs], seur@tsne.rot[cells.use,], perplexity=perplexity, n.cores=ncores)
colnames(seur@tsne.rot) = c('tSNE_1', 'tSNE_2')
