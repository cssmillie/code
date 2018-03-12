require(Seurat)
require(Matrix)
require(methods)
require(data.table)
require(tibble)

source('~/code/single_cell/batch.r')
source('~/code/single_cell/cluster.r')
source('~/code/single_cell/contamination.r')
source('~/code/single_cell/dmap.r')
source('~/code/single_cell/downsample.r')
source('~/code/single_cell/gsea.r')
source('~/code/single_cell/map.r')
source('~/code/single_cell/markers.r')
source('~/code/single_cell/parallel.r')
source('~/code/single_cell/pca.r')
source('~/code/single_cell/plot.r')
source('~/code/single_cell/specific.r')
source('~/code/single_cell/tsne.r')
source('~/code/single_cell/var_genes.r')
source('~/code/util/mtx.r')

msg = function(name, text, verbose){
    if(verbose == TRUE){
        print(paste(name, text))
    }
}


make_seurat = function(name, dge=NULL, regex='', regexv='', minc=10, maxc=1e6, ming=500, maxg=1e6, genes.use=NULL, cells.use=NULL, ident_fxn=NULL, verbose=FALSE, x11=FALSE){
    
    # Load packages
    require(data.table)
    require(Matrix)
    source('~/code/util/mtx.r')
    
    # Set graphics device
    options(device=pdf)

    # Load counts from Seurat object, matrix, or file
    msg(name, 'Loading DGE', verbose)
    if(typeof(dge) == typeof('')){
        if(file.exists(dge)){
	    counts = fread(paste('zcat', dge))
	    counts = data.frame(counts, row.names=1)
	} else {
	    counts = read_mtx(prefix=dge)
	}
    } else {
	counts = dge
    }
    msg(name, sprintf('DGE = %d x %d', nrow(counts), ncol(counts)), verbose)
    
    # Subset DGE with NA, regex, genes.use, and cells.use
    msg(name, 'Subsetting DGE', verbose)
    if(regex != ''){
        j = grep(regex, colnames(counts))
	counts = counts[,j,drop=F]
    }
    if(regexv != ''){
        j = grep(regexv, colnames(counts), invert=T)
	counts = counts[,j,drop=F]
    }
    if(!is.null(genes.use)){
	genes.use = intersect(rownames(counts), genes.use)
	counts = counts[genes.use,,drop=F]
    }
    if(!is.null(cells.use)){
	cells.use = intersect(colnames(counts), cells.use)
	counts = counts[,cells.use,drop=F]
    }
    genes.use = rowSums(is.na(counts)) == 0
    counts = counts[genes.use,,drop=F]
    
    msg(name, sprintf('DGE = %d x %d', nrow(counts), ncol(counts)), verbose)
    
    # Convert counts to sparse matrix
    counts = as(as.matrix(counts), 'sparseMatrix')
        
    # Get cell identities
    if(is.null(ident_fxn)){
        ident = sapply(strsplit(colnames(counts), '\\.'), '[', 1)
    } else {
        ident = sapply(colnames(counts), ident_fxn)
    }
    print(table(ident))

    # Select maxc cells from each group
    if(any(maxc < table(ident))){
        msg(name, paste('Selecting', maxc, 'cells from each group'))
        j = as.character(unlist(tapply(colnames(counts), list(ident), function(a) sample(a, min(length(a), maxc)))))
        counts = counts[,j]
    }
    
    # Filter cells by minc, ming, maxg
    msg(name, 'Filtering DGE', verbose)
    j1 = colSums(counts > 0) >= ming
    j2 = colSums(counts > 0) <= maxg
    counts = counts[,(j1 & j2)]
    i = rowSums(counts > 0) >= minc
    counts = counts[i,]
    msg(name, sprintf('DGE = %d x %d', nrow(counts), ncol(counts)), verbose)
    
    # Make Seurat object
    msg(name, 'Making Seurat object', verbose)
    seur = new('seurat', raw.data=counts)
    seur = Setup(seur, project=name, min.cells=0, min.genes=0, total.expr=1e4, names.delim='\\.', names.field=1, do.center=F, do.scale=F, do.logNormalize=F)
    
    # Convert data to log2(TPM)
    seur@data = as(log2(1e4*scale(seur@raw.data, center=F, scale=colSums(seur@raw.data)) + 1), 'sparseMatrix')
    
    # Get cell identities
    if(!is.null(ident_fxn)){
        ident = sapply(colnames(seur@data), ident_fxn)
	seur = set.ident(seur, ident.use=ident)
	seur@data.info$orig.ident = seur@ident
    }
    
    if(length(unique(seur@ident)) > 100){
        seur@ident = '1'
	seur@orig.ident = '1'
    }
    
    print(table(seur@ident))
    return(seur)
}


run_seurat = function(name, seur=NULL, dge=NULL, regex='', regexv='', cells.use=NULL, genes.use=NULL, minc=5, maxc=1e6, ming=500, maxg=1e6, ident_fxn=NULL, varmet='loess', var_regexv=NULL,
             min_cv2=.25, var_genes=NULL,
	     num_genes=1500, do.batch='none', batch.use=NULL, design=NULL, pc.data=NULL, num_pcs=0, robust_pca=F, perplexity=25, max_iter=1000, dist.use='cosine', do.largevis=FALSE, largevis.k=50,
	     cluster='infomap', k=c(), verbose=T, write_out=T, do.backup=F, ncores=1, stop_cells=50, marker.test=''){

    # check input arguments
    if(! do.batch %in% c('none', 'combat', 'mnn', 'cca')){stop('do.batch must be none, combat, mnn, or cca')}
    
    # Make Seurat object
    if(is.null(seur)){
        seur = make_seurat(name=name, dge=dge, regex=regex, regexv=regexv, minc=minc, maxc=maxc, ming=ming, maxg=maxg, genes.use=genes.use, cells.use=cells.use, ident_fxn=ident_fxn, verbose=verbose)
    }
    if(ncol(seur@data) <= stop_cells){return(seur)}
    
    msg(name, 'Selecting variable genes', verbose)
    ident = seur@ident
    if(is.null(var_genes)){var_genes = get_var_genes(seur@raw.data, ident=ident, method=varmet, num_genes=num_genes, min_ident=25)}
    if(!is.null(var_regexv)){var_genes = grep(var_regexv, var_genes, invert=T, value=T)}
    msg(name, sprintf('Found %d variable genes', length(var_genes)), verbose)
    seur@var.genes = intersect(var_genes, rownames(seur@data))
    
    # Batch correction
    if(do.batch != 'none'){
	msg(name, 'Batch correction', verbose)
	if(is.null(batch.use)){
	    batch.use = seur@ident
	} else {
	    batch.use = batch.use[names(seur@ident),1]
	}
	if(!is.null(design)){
	    design = design[names(seur@ident),,drop=F]
	}
	bc.data = batch_correct(seur, batch.use, design=design, method=do.batch, genes.use='var')
	pc.data = t(scale(t(bc.data)))
    }
    if(is.null(pc.data)){pc.data = seur@data}
    
    # Number of significant PCs (stored in seur@data.info$num_pcs)
    if(num_pcs == 0){
	num_pcs = sig.pcs.perm(scale(t(seur@data)), randomized=T, n.cores=ncores)$r + 2
    }
    if(is.na(num_pcs)){num_pcs = 5}
    seur@data.info$num_pcs = num_pcs
    msg(name, sprintf('Found %d significant PCs', num_pcs), verbose)

    # Fast PCA on data
    seur = run_rpca(seur, data=pc.data, k=50, genes.use=seur@var.genes, robust=robust_pca, rescale=T)
    msg(name, 'PC loadings', verbose)
    loaded_genes = get.loaded.genes(seur@pca.obj[[1]], components=1:num_pcs, n_genes=20)
    print(loaded_genes)

    # TSNE
    knn = NULL
    if(do.largevis == FALSE){
        msg(name, 'TSNE', verbose)
	require(Rtsne)
        if(dist.use == 'euclidean'){
            q = Rtsne(seur@pca.rot[,1:num_pcs], do.fast=T, max_iter=max_iter, perplexity=perplexity, verbose=T)$Y
        } else if(dist.use == 'cosine'){
            d = cosine_dist(t(seur@pca.rot[,1:num_pcs]))
	    q = Rtsne(d, is_distance=T, do.fast=T, max_iter=max_iter, perplexity=perplexity, verbose=T)$Y
        } else {
            d = dist(seur@pca.rot[,1:num_pcs], method=dist.use)
	    q = Rtsne(d, is_distance=T, do.fast=T, max_iter=max_iter, perplexity=perplexity, verbose=T)$Y
        }
    } else {
        msg(name, 'largeVis', verbose)
	q = run_largevis(seur@pca.rot[,1:num_pcs], k=largevis.k, save_knn=TRUE, save_weights=FALSE, dist.use=dist.use, verbose=T)
	knn = q$knn
	q = q$coords
    }
    rownames(q) = colnames(seur@data)
    colnames(q) = c('tSNE_1', 'tSNE_2')
    seur@tsne.rot = as.data.frame(q)
    
    msg(name, 'Calculate signatures', verbose)
    seur = update_signatures(seur)

    # Cluster cells and run DE tests
    if(length(k) > 0){

        # Save backup Seurat object
	if(do.backup){saveRDS(seur, file=paste0(name, '.seur.rds'))}
	
	msg(name, 'Clustering cells', verbose)
	k = k[k < ncol(seur@data)]
	u = paste('Cluster.Infomap.', k, sep='')
	v = run_cluster(seur@pca.rot[,1:num_pcs], k, method=cluster, weighted=FALSE, n.cores=min(length(k), ncores), dist='cosine', do.fast=T, knn=knn)
	seur@data.info[,u] = v

	if(marker.test != ''){
    	    msg(name, 'Differential expression', verbose)
            covariates = subset(seur@data.info, select=c(nGene, Cell_Cycle))
	    seur = set.ident(seur, ident.use=v[,1])
	    print(table(seur@ident))
	    markers = lapply(k, function(ki){
	        seur = set.ident(seur, ident.use=seur@data.info[,paste0('Cluster.Infomap.', ki)])
	        markers = p.find_all_markers(seur, test.use=marker.test)
	    })
	    names(markers) = k
	}
    }

    if(write_out){

        # Batch corrected data
	if(do.batch != 'none'){fwrite(as.data.table(bc.data), file=paste0(name, '.bc.dge.txt'), sep='\t')}
	
	# PCA loaded genes
	write.table(loaded_genes, paste0(name, '.loaded_genes.txt'), sep='\t', quote=F)

	# Plot TSNE
	png(paste0(name, '.tsne.png'), width=800, height=650)
	tsne.plot(seur, pt.size=1)
	dev.off()

	# Plot clusters
	if(length(k) > 0){
	    pdf(paste0(name, '.clusters.pdf'), width=9, height=9)
	    plot_clusters(seur)
	    dev.off()
	}

	# Marker genes
	if(marker.test != ''){
	    for(ki in names(markers)){write.table(markers[[ki]], file=paste0(name, '.k', ki, '.', marker.test, '.txt'), sep='\t', quote=F)}
	}

	# Plot summary statistics
	png(paste0(name, '.summary_stats.png'), width=1200, height=800)
	plot_tsne(seur, subset(seur@data.info, select=c(orig.ident, G1S, G2M, nGene)), do.label=F)
	dev.off()

	# Save Seurat object
	saveRDS(seur, file=paste0(name, '.seur.rds'))
    }

    return(seur)
}


safeRDS = function(object, file){
    temp = paste0(file, '.temp')
    saveRDS(object, file=temp)
    system(paste('mv', temp, file))
}


merge_seurat = function(seur1, seur2, rm_old=FALSE){
    
    # Merge objects
    seur = MergeSeurat(seur1, seur2, do.scale=F, do.center=F, do.logNormalize=F)
    
    # Fix idents
    ident = c(seur1@ident, seur2@ident)[colnames(seur@data)]
    seur = SetIdent(seur, ident.use=ident)

    # Remove old seurat objects (save memory)
    if(rm_old == TRUE){rm(seur1); rm(seur2)}
    
    # Calculate log2(TPM + 1)
    seur@data = log2(calc_tpm(seur=seur) + 1)
    
    # Return merged object
    return(seur)
}


make_mini = function(seur, num_genes=100, num_cells=100, ident.k=NULL){

    # Select genes and cells
    genes.use = rownames(seur@data)[order(rowMeans(seur@data), decreasing=T)[1:num_genes]]
    cells.use = sample(colnames(seur@data), num_cells)

    # Construct mini Seurat object
    mini = make_seurat(seur=seur, name='mini', genes.use=genes.use, cells.use=cells.use, ming=0, minc=0)

    # Set random identities
    if(!is.null(ident.k)){
        ident = sample(1:ident.k, ncol(mini@data), replace=T)
	mini = set.ident(mini, ident.use=ident)
    }
    return(mini)    
}


map_ident = function(seur, old_ident){    
    old_ident = as.data.frame(old_ident)
    new_ident = data.frame(ident=rep(NA, ncol(seur@data)), row.names=colnames(seur@data))
    i = intersect(rownames(old_ident), rownames(new_ident))
    new_ident[i,1] = old_ident[i,1]
    new_ident = structure(as.factor(new_ident[,1]), names=rownames(new_ident))
    return(new_ident)
}


update_signatures = function(seur){
    
    # Predict host (human or mouse)
    genes = rownames(seur@data)
    
    # Calculate signatures
    seur@data.info = cbind(seur@data.info, get_scores(seur, file='cell_cycle'))
    seur@data.info = cbind(seur@data.info, get_scores(seur, file='early'))
    seur@data.info$nGene = colSums(seur@data > 0)
    
    return(seur)
}
