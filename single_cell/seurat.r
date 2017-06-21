library(Seurat)
library(Matrix)
library(methods)

source('~/code/single_cell/batch.r')
source('~/code/single_cell/cluster.r')
source('~/code/single_cell/map.r')
source('~/code/single_cell/markers.r')
source('~/code/single_cell/parallel.r')
source('~/code/single_cell/pca.r')
source('~/code/single_cell/plot.r')
source('~/code/single_cell/tsne.r')
source('~/code/single_cell/var_genes.r')


msg = function(name, text, verbose){
    if(verbose == TRUE){
        print(paste(name, text))
    }
}


tpm = function(seur, data=NULL, genes.use=NULL, cells.use=NULL, total=1e4){

    # Calculate TPM from @raw.data or data argument

    # Get counts data for calculation
    if(is.null(data)){data = seur@raw.data}
    if(is.null(genes.use)){genes.use = rownames(seur@raw.data)}
    if(is.null(cells.use)){cells.use = colnames(seur@raw.data)}
    
    # Intersect genes.use and cells.use
    genes.use = intersect(genes.use, rownames(data))
    cells.use = intersect(cells.use, colnames(data))
    
    # Get number of UMIs for normalization
    if('nUMI' %in% colnames(seur@data.info)){
        numi = seur@data.info[cells.use,'nUMI']
    } else {
        numi = colSums(data[,cells.use])
    }
    
    # Calculate TPM
    tpm = total*scale(data[genes.use, cells.use, drop=F], center=F, scale=numi)
    return(tpm)
}


make_seurat = function(name, seur=NULL, dge=NULL, regex='', minc=10, ming=500, maxg=4000, genes.use=NULL, cells.use=NULL, ident_fxn=NULL, verbose=FALSE, x11=FALSE){
    
    # Load packages
    library(data.table)
    library(Matrix)
    
    # Set graphics device
    options(device=pdf)
    
    # Load counts from Seurat object, matrix, or file
    msg(name, 'Loading DGE', verbose)
    if(!is.null(seur)){
	counts = seur@raw.data
    } else if(!is.null(dge)){
        if(typeof(dge) == typeof('')){
	    counts = fread(paste('zcat', dge))
	    counts = data.frame(counts, row.names=1)
	} else {
	    counts = dge
	}
    } else {stop('Must specify seurat object or dge location')}
    msg(name, sprintf('DGE = %d x %d', nrow(counts), ncol(counts)), verbose)

    
    # Subset DGE with regex, genes.use, and cells.use
    msg(name, 'Subsetting DGE', verbose)
    j = grep(regex, colnames(counts))
    counts = counts[,j]
    if(!is.null(genes.use)){
	genes.use = intersect(rownames(counts), genes.use)
	counts = counts[genes.use,]
    }
    if(!is.null(cells.use)){
	cells.use = intersect(colnames(counts), cells.use)
	counts = counts[,cells.use]
    }
    msg(name, sprintf('DGE = %d x %d', nrow(counts), ncol(counts)), verbose)

    # Convert counts to sparse matrix
    counts = as(as.matrix(counts), 'sparseMatrix')
    
    # Filter cells by minc, ming, maxg
    msg(name, 'Filtering DGE', verbose)
    j1 = colSums(counts > 0) >= ming
    j2 = colSums(counts > 0) <= maxg
    counts = counts[, (j1 & j2)]
    i = rowSums(counts > 0) >= minc
    counts = counts[i,]
    msg(name, sprintf('DGE = %d x %d', nrow(counts), ncol(counts)), verbose)
    
    # Make Seurat object
    msg(name, 'Making Seurat object', verbose)
    seur = new('seurat', raw.data=counts)
    seur = Setup(seur, project=name, min.cells=0, min.genes=0, total.expr=1e4, names.delim='\\.', names.field=1)
    
    # Convert data to log2(TPM)
    seur@data = as(log2(1e4*scale(seur@raw.data, center=F, scale=colSums(seur@raw.data)) + 1), 'sparseMatrix')
    seur@scale.data = t(scale(t(seur@data)))
    
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


run_seurat = function(name, seur=NULL, dge=NULL, regex='', cells.use=NULL, genes.use=NULL, minc=10, maxc=0, ming=500, maxg=4000, ident_fxn=NULL, varmet='karthik', min_cv2=.25, num_genes=1500,
	     do.batch=F, batch.use=NULL, batch.replace=FALSE, design=NULL, num_pcs=0, robust_pca=F, perplexity=25, max_iter=1000, tsne_cor=F, cluster='infomap', k=c(), verbose=T, write_out=T,
	     do.backup=F, ncores=1, stop_cells=50, marker.test=''){

    seur = make_seurat(name=name, seur=seur, dge=dge, regex=regex, minc=minc, ming=ming, maxg=maxg, genes.use=genes.use, cells.use=cells.use, ident_fxn=ident_fxn, verbose=verbose)
    if(ncol(seur@data) <= stop_cells){return(seur)}
    
    msg(name, 'Selecting variable genes', verbose)
    var_genes = select_var_genes(seur@raw.data, method=varmet, vcut=NULL, min.cv2=min_cv2, num_genes=num_genes)
    msg(name, sprintf('Found %d variable genes', length(var_genes)), verbose)
    seur@var.genes = intersect(var_genes, rownames(seur@data))
    
    # Batch correction
    if(do.batch){
	msg(name, 'Batch correction', verbose)
	if(is.null(batch.use)){
	    batch.use = seur@ident
	} else {
	    batch.use = batch.use[names(seur@ident),1]
	    design = design[names(seur@ident),1]
	}
	bc.data = batch_correct(seur, batch.use, design=design, replace=batch.replace)
	pc.data = t(scale(t(bc.data)))
	if(batch.replace){
	    seur@data = bc.data
	    seur@scale.data = pc.data
	}
    } else {
        pc.data = seur@scale.data
    }
    
    # Number of significant PCs (stored in seur@data.info$num_pcs)
    if(num_pcs == 0){
        num_pcs = sig.pcs.perm(t(pc.data[seur@var.genes,]), randomized=T, n.cores=ncores)$r
    }
    if(is.na(num_pcs) | num_pcs <= 1){num_pcs = 3}
    seur@data.info$num_pcs = num_pcs
    msg(name, sprintf('Found %d significant PCs', num_pcs), verbose)

    # Fast PCA on data
    seur = run_rpca(seur, data=pc.data, k=20, genes.use=seur@var.genes, robust=robust_pca)
    
    msg(name, 'TSNE', verbose)
    if(tsne_cor == F){
        seur = RunTSNE(seur, dims.use=1:num_pcs, do.fast=T, max_iter=max_iter, verbose=T, perplexity=perplexity)
    } else {
        library(Rtsne)
	#d = as.dist(1 - cor(t(seur@pca.rot[,1:num_pcs])))
	d = cosine_dist(t(seur@pca.rot[,1:num_pcs]))
	q = Rtsne(d, is_distance=T, do.fast=T, max_iter=max_iter, verbose=T, perplexity=perplexity)$Y
	rownames(q) = colnames(seur@data)
	colnames(q) = c('tSNE_1', 'tSNE_2')
	seur@tsne.rot = as.data.frame(q)
    }
    
    msg(name, 'PC loadings', verbose)
    loaded_genes = get.loaded.genes(seur@pca.obj[[1]], components=1:num_pcs, n_genes=20)
    print(loaded_genes)

    msg(name, 'Calculate signatures', verbose)
    seur = update_signatures(seur)

    # Cluster cells and run DE tests
    if(length(k) > 0){

        # Save backup Seurat object
	if(do.backup){saveRDS(seur, file=paste0(name, '.seur.rds'))}
	
	msg(name, 'Clustering cells', verbose)
	k = k[k < ncol(seur@data)]
	u = paste('Cluster.Infomap.', k, sep='')
	v = run_cluster(seur@pca.rot[,1:num_pcs], k, method=cluster, weighted=FALSE, n.cores=min(length(k), ncores), dist='cosine', do.fast=T)
	seur@data.info[,u] = v

	msg(name, 'Differential expression', verbose)
	if(marker.test != ''){
	    covariates = subset(seur@data.info, select=c(nGene, Cell_Cycle))
	    markers = lapply(k, function(ki){
	        seur = set.ident(seur, ident.use=seur@data.info[,paste0('Cluster.Infomap.', ki)])
	        markers = p.find_all_markers(seur, covariates=covariates, test.use=marker.test)
	    })
	    names(markers) = k
	}
    }

    if(write_out){

        # Batch corrected data
	if(do.batch){write.table(bc.data, file=paste0(name, '.bc.dge.txt'), sep='\t', quote=F)}

	# PCA loaded genes
	write.table(loaded_genes, paste0(name, '.loaded_genes.txt'), sep='\t', quote=F)

	# Plot TSNE
	png(paste0(name, '.tsne.png'), width=800, height=650)
	tsne.plot(seur, pt.size=1, label.cex.text=.25)
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
	    for(ki in names(markers)){write.table(marker_table(markers[[ki]], pval=.05, top=50), file=paste0(name, '.k', ki, '.', marker.test, '.table.txt'), sep='\t', quote=F)}
	}

	# Plot summary statistics
	png(paste0(name, '.summary_stats.png'), width=1200, height=800)
	plot_tsne(seur, subset(seur@data.info, select=c(orig.ident, G1S, G2M, nGene, nUMI)), do.label=F)
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
