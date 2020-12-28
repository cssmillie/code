require(Seurat)
require(Matrix)
require(methods)
require(data.table)
require(tibble)
require(naturalsort)

source('~/code/single_cell/batch.r')
source('~/code/single_cell/cluster.r')
source('~/code/single_cell/contamination.r')
source('~/code/single_cell/dmap.r')
source('~/code/single_cell/downsample.r')
source('~/code/single_cell/frequencies.r')
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
source('~/code/FIt-SNE/fast_tsne.R')

msg = function(name, text, verbose){
    if(verbose == TRUE){
        print(paste(name, text))
    }
}


make_seurat = function(name, dge=NULL, regex='', regexv='', minc=10, maxc=NULL,maxc_per_group=NULL, ming=500, maxg=1e6, genes.use=NULL, cells.use=NULL, ident_fxn=NULL, verbose=FALSE, x11=FALSE, qnorm=F){
    
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
    if(is.data.frame(counts)){counts = as.matrix(counts)}
    counts = as(counts, 'sparseMatrix')
            
    # Get cell identities
    if(is.null(ident_fxn)){
        ident = sapply(strsplit(colnames(counts), '\\.'), '[', 1)
    } else {
        ident = sapply(colnames(counts), ident_fxn)
    }
        
    # Downsample cells
    if(!is.null(maxc) | !is.null(maxc_per_group)){
        cells.use = simple_downsample(cells=colnames(counts), groups=ident, total_cells=maxc, cells_per_group=maxc_per_group)
	msg(name, paste('Downsampling to', length(cells.use), 'total cells'), verbose)
	counts = counts[,cells.use]
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
        
    # Normalize data (default = TPM)
    msg(name, 'Normalizing data', verbose)
    if(qnorm == FALSE){
        seur@data = calc_tpm(counts=seur@raw.data)
	seur@data@x = log2(seur@data@x + 1)
    } else {
        seur = quantile_normalize(seur=seur)
    }
        
    # Get cell identities
    msg(name, 'Setting cell identities', verbose)
    if(!is.null(ident_fxn)){
        ident = sapply(colnames(seur@data), ident_fxn)
	seur = set.ident(seur, ident.use=ident)
	seur@data.info$orig.ident = seur@ident
    }
    
    if(length(unique(seur@ident)) > 100){
        msg(name, 'WARNING: nlevels(seur@ident) > 100', verbose)
        #seur@ident = '1'
    	#seur@data.info$orig.ident = '1'
    }
    
    print(table(seur@ident))
    return(seur)
}

quantile_normalize = function(seur=NULL, data=NULL, zero_cut=1e-4){
    print(paste0('Quantile normalization with zero_cut = ', zero_cut))
    library(preprocessCore)
    if(!is.null(seur)){
        print(dim(seur@raw.data))
        seur@raw.data[] = as(normalize.quantiles.robust(as.matrix(seur@raw.data)), 'sparseMatrix')
	seur@data[] = as(normalize.quantiles.robust(as.matrix(seur@data)), 'sparseMatrix')
        seur
    } else {
        print(dim(data))
        data[] = normalize.quantiles.robust(as.matrix(data))
	print(paste('Flooring', sum(data < zero_cut), 'zeros'))
	data[data < zero_cut] = 0
	as(data, 'sparseMatrix')
    }
}


run_seurat = function(name, seur=NULL, dge=NULL, regex='', regexv='', cells.use=NULL, genes.use=NULL, minc=5, maxc=1e6, ming=200, maxg=1e6, ident_fxn=NULL, varmet='loess', var_regexv=NULL,
             var_remove=NULL, min_cv2=.25, var_genes=NULL, qnorm=F, num_genes=1500, do.batch='none', batch.use=NULL,
	     design=NULL, pc.data=NULL, num_pcs=0, pcs.use=NULL, pcs.rmv=NULL, robust_pca=F,
	     perplexity=25, max_iter=1000, dist.use='euclidean', do.largevis=FALSE, do.umap=FALSE, largevis.k=50, do.fitsne=FALSE, fitsne.K=-1,
	     cluster='infomap', k=c(), verbose=T, write_out=T, do.backup=F, ncores=1, stop_cells=50, marker.test=''){

    # check input arguments
    if(! do.batch %in% c('none', 'combat', 'mnn', 'cca', 'multicca', 'liger')){stop('do.batch must be none, combat, mnn, or cca')}
    if(! is.null(batch.use)){if(is.null(names(batch.use))){stop('batch.use needs names')}}
    
    # Make Seurat object
    if(is.null(seur)){
        seur = make_seurat(name=name, dge=dge, regex=regex, regexv=regexv, minc=minc, maxc=maxc, ming=ming, maxg=maxg, genes.use=genes.use, cells.use=cells.use, ident_fxn=ident_fxn, verbose=verbose, qnorm=qnorm)
    }
    if(ncol(seur@data) <= stop_cells){return(seur)}
    
    msg(name, 'Selecting variable genes', verbose)
    ident = seur@ident
    if(is.null(var_genes)){
        gi = rownames(seur@raw.data)
	print(paste('Starting with', length(gi), 'genes'))
	if(!is.null(var_regexv)){gi = grep(var_regexv, gi, invert=T, value=T)}
	print(paste('var_regexv:', length(gi), 'genes'))
	if(!is.null(var_remove)){gi = setdiff(gi, var_remove)}
	print(paste('var_remove:', length(gi), 'genes'))
	var_genes = get_var_genes(seur@raw.data, ident=ident, method=varmet, genes.use=genes.use, num_genes=num_genes, min_ident=25)
    }
    #if(is.null(var_genes)){var_genes = get_var_genes(seur@raw.data, ident=ident, method=varmet, num_genes=num_genes, min_ident=25)}
    #if(!is.null(var_regexv)){var_genes = grep(var_regexv, var_genes, invert=T, value=T)}
    msg(name, sprintf('Found %d variable genes', length(var_genes)), verbose)
    seur@var.genes = intersect(var_genes, rownames(seur@data))
    print(var_genes)
    
    # Batch correction with variable genes
    if(do.batch != 'none'){
	msg(name, 'Batch correction', verbose)
	if(is.null(batch.use)){
	    batch.use = seur@ident
	}
	batch.use = batch.use[names(seur@ident)]
	print(table(batch.use))
	if(!is.null(design)){
	    design = design[names(seur@ident),,drop=F]
	}
	bc.data = batch_correct(seur, batch.use, design=design, method=do.batch, genes.use=seur@var.genes, ndim=num_pcs)
	
	# write batch corrected data to file
	if(write_out == TRUE){fwrite(as.data.table(bc.data), file=paste0(name, '.bc.data.txt'), sep='\t')}
	
	pc.data = t(scale(t(bc.data), center=F))
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
    if(do.batch %in% c('liger', 'multicca')){
        seur@pca.rot = as.data.frame(pc.data) # liger and multi-cca output saved in pc.data
	print(dim(seur@pca.rot))
    } else {
        seur = run_rpca(seur, data=pc.data, k=50, genes.use=seur@var.genes, robust=robust_pca, rescale=T)
	if(write_out == TRUE){saveRDS(seur@pca.obj, file=paste0(name, '.pca.rds'))}    	
        msg(name, 'PC loadings', verbose)
        loaded_genes = get.loaded.genes(seur@pca.obj[[1]], components=1:num_pcs, n_genes=20)
        print(loaded_genes)
    }
    
    # Fix problem with duplicates
    seur@pca.rot[,num_pcs] = seur@pca.rot[,num_pcs] + runif(nrow(seur@pca.rot), min=-1e-8, max=1e-8)
    
    # Regress out PCs
    if(is.null(pcs.use)){
        pcs.use = 1:(min(ncol(seur@pca.rot), num_pcs))
    }
    
    # TSNE
    knn = NULL
    if(max_iter > 0){
    if(do.fitsne == TRUE){
        msg(name, 'FIt-SNE', verbose)
	q = fftRtsne(seur@pca.rot[,pcs.use], max_iter=max_iter, perplexity=perplexity, K=fitsne.K, fast_tsne_path='~/code/FIt-SNE/bin/fast_tsne')
    } else if(do.largevis == TRUE){
        msg(name, 'largeVis', verbose)
	q = run_largevis(seur@pca.rot[,pcs.use], k=largevis.k, save_knn=TRUE, save_weights=FALSE, dist.use=dist.use, verbose=T)
	knn = q$knn
	q = q$coords
    } else if(do.umap == TRUE){
        msg(name, 'umap', verbose)
	library(umap)
	q = umap(seur@pca.rot[,pcs.use], method='umap-learn')$layout
    } else {
        msg(name, 'TSNE', verbose)
	require(Rtsne)
        if(dist.use == 'euclidean'){
            q = Rtsne(seur@pca.rot[,pcs.use], do.fast=T, max_iter=max_iter, perplexity=perplexity, verbose=T)$Y
        } else if(dist.use == 'cosine'){
            d = cosine_dist(t(seur@pca.rot[,pcs.use]))
	    q = Rtsne(d, is_distance=T, do.fast=T, max_iter=max_iter, perplexity=perplexity, verbose=T)$Y
        } else {
            d = dist(seur@pca.rot[,pcs.use], method=dist.use)
	    q = Rtsne(d, is_distance=T, do.fast=T, max_iter=max_iter, perplexity=perplexity, verbose=T)$Y
        }
    }
    rownames(q) = colnames(seur@data)
    colnames(q) = c('tSNE_1', 'tSNE_2')
    seur@tsne.rot = as.data.frame(q)
    }
    
    # Cluster cells and run DE tests
    if(length(k) > 0){

        # Save backup Seurat object
	if(do.backup){saveRDS(seur, file=paste0(name, '.seur.rds'))}
	
	msg(name, 'Clustering cells', verbose)
	k = k[k < ncol(seur@data)]
	u = paste('Cluster.Infomap.', k, sep='')
	v = run_cluster(seur@pca.rot[,pcs.use], k, method=cluster, weighted=FALSE, n.cores=min(length(k), ncores), dist='cosine', do.fast=T, knn=knn)
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


merge_seurat = function(seur1, seur2, rm_old=FALSE, mem=FALSE, id1=NULL, id2=NULL){

    # Fix cell names
    if(!is.null(id1)){
	colnames(seur1@raw.data) = paste(id1, colnames(seur1@raw.data), sep='.')
	colnames(seur1@data) = paste(id1, colnames(seur1@data), sep='.')
	rownames(seur1@data.info) = paste(id1, rownames(seur1@data.info), sep='.')
	rownames(seur1@tsne.rot) = paste(id1, rownames(seur1@tsne.rot), sep='.')
	names(seur1@ident) = paste(id1, names(seur1@ident), sep='.')
    }
    if(!is.null(id2)){
	colnames(seur2@raw.data) = paste(id2, colnames(seur2@raw.data), sep='.')
	colnames(seur2@data) = paste(id2, colnames(seur2@data), sep='.')	
	rownames(seur2@data.info) = paste(id2, rownames(seur2@data.info), sep='.')
	rownames(seur2@tsne.rot) = paste(id2, rownames(seur2@tsne.rot), sep='.')
	names(seur2@ident) = paste(id1, names(seur2@ident), sep='.')	
    }
    
    # Calculate idents
    ident = setNames(c(as.character(seur1@ident), as.character(seur2@ident)), c(colnames(seur1@data), colnames(seur2@data)))
    ident = factor(ident, levels=c(levels(seur1@ident), levels(seur2@ident)))
            
    # Get metadata
    rows = c(rownames(seur1@data.info), rownames(seur2@data.info))
    cols = sort(unique(c(colnames(seur1@data.info), colnames(seur2@data.info))))
    meta = matrix(NA, nrow=length(rows), ncol=length(cols))
    rownames(meta) = rows
    colnames(meta) = cols
    meta[rownames(seur1@data.info), colnames(seur1@data.info)] = as.matrix(seur1@data.info)
    meta[rownames(seur2@data.info), colnames(seur2@data.info)] = as.matrix(seur2@data.info)
    meta = as.data.frame(meta)
    
    # Get tsne coordinates
    tsne = rbind(seur1@tsne.rot[,1:2], seur2@tsne.rot[,1:2])
    
    if(mem == FALSE){
        
        # Merge objects
    	seur = MergeSeurat(seur1, seur2, do.scale=F, do.center=F, do.logNormalize=F, add.cell.id1=id1, add.cell.id2=id2)
		
    	# Remove old seurat objects (save memory)
    	if(rm_old == TRUE){rm(seur1); rm(seur2)}
    	
    	# Calculate log2(TPM + 1)
    	seur@data = log2(calc_tpm(seur=seur) + 1)
    	seur@ident = ident[colnames(seur@data)]
    
    } else {

	# Make counts matrix
	print('Merge counts')
	counts = mem_cbind(list(seur1@raw.data, seur2@raw.data))
	
	# Remove old seurat objects (save memory)
	if(rm_old == TRUE){rm(seur1); rm(seur2)}
	
	# Make Seurat object
	seur = new('seurat', raw.data=counts)
	print('Calculate TPM')
        seur@data = calc_tpm(counts=seur@raw.data)
	print('Log transform')
        seur@data@x = log2(seur@data@x + 1)
	seur@ident = ident[colnames(seur@data)]
    }
    
    # Add metadata
    seur@data.info = as.data.frame(meta)
    seur@tsne.rot = as.data.frame(tsne)
        
    # Return merged object
    return(seur)
}


merge_seurat_hm = function(seurs, target='human', ident_fxn=NULL, ortholog_filter='hm'){

    # Merge seurat objects, using only 1:1 orthologs in human and mouse
    # -----------------------------------------------------------------
    # Input arguments:
    # - seurs = list of seurat objects ('human' or 'mouse')
    # - target = organism to map genes to ('human' or 'mouse')
    # - ident_fxn = identity function
    # - ortholog_filter:
    #   - 'all' = require orthologs to be found in all seurat objects
    #   - 'hm'  = require orthologs to be found in human and mouse
    #   - 'none' = no ortholog filter
    
    # Fix input arguments
    seurs = as.list(seurs)
    orgs = sapply(seurs, function(a) predict_organism(rownames(a@data)))
    print(names(seurs))
    print(orgs)
    
    # Calculate metadata
    cells.use = unname(unlist(sapply(seurs, function(a) colnames(a@data))))
    ident.use = unname(unlist(sapply(names(seurs), function(a) paste(a, as.character(seurs[[a]]@ident), sep='.'))))
    levels(ident.use) = unname(unlist(sapply(names(seurs), function(a) paste(a, as.character(levels(seurs[[a]]@ident)), sep='.'))))
    org.use = unlist(sapply(1:length(seurs), function(i) rep(orgs[[i]], ncol(seurs[[i]]@data))))
    dset.use = unlist(sapply(names(seurs), function(a) rep(a, ncol(seurs[[a]]@data))))
    names(ident.use) = names(org.use) = names(dset.use) = cells.use
        
    # Define human and mouse gene sets
    h_genes = readLines('~/aviv/db/map_gene/hg19_genes.txt')
    m_genes = readLines('~/aviv/db/map_gene/mm10_genes.txt')
        
    # Map genes
    ortho_fn = '~/aviv/db/map_gene/hm_orthologs.rds'
    if(file.exists(ortho_fn)){
        print('Reading 1:1 orthologs')
        res = readRDS(ortho_fn)
	h2m = res$h2m
	m2h = res$m2h
    } else {
        print('Calculating 1:1 orthologs')
        h2m = sapply(h_genes, function(gene) {
            mi = map_gene(gene, source='human', target='mouse', do.unlist=F)[[1]]
            hi = map_gene(gene, source='mouse', target='human', do.unlist=F)[[1]]
            if(length(mi) == 1 & length(hi) == 1){
                if(gene == hi){
                    mi
                }
            }
        })
        h2m = sapply(h2m[lengths(h2m) == 1], '[[', 1)
        m2h = setNames(names(h2m), h2m)
	saveRDS(list(h2m=h2m, m2h=m2h), file=ortho_fn)
    }
    
    # Merge datasets
    print('Merging counts')
    counts = sapply(1:length(seurs), function(i) {
        
        # Get genes to use
	if(orgs[[i]] == 'human'){
	    genes.use = intersect(names(h2m), rownames(seurs[[i]]@data))
	} else {
	    genes.use = intersect(names(m2h), rownames(seurs[[i]]@data))
	}

	# Subset data
	ci = seurs[[i]]@raw.data[genes.use,]
	
	# Fix rownames
	if(orgs[[i]] == 'human' & target == 'mouse'){
	    rownames(ci) = unname(h2m[rownames(ci)])
	}
	if(orgs[[i]] == 'mouse' & target == 'human'){
	    rownames(ci) = unname(m2h[rownames(ci)])
	}
	
	ci

    }, simplify=F)
    counts = sparse_cbind(as.list(counts))
    print(dim(counts))
    
    # Filter genes
    if(ortholog_filter != 'none'){
        print('Filtering genes')
        if(ortholog_filter == 'all'){
	    groups = dset.use
	} else {
	    groups = org.use
	}
	genes.use = sapply(unique(groups), function(a){
	    j = (groups == a)
	    apply(counts[,j], 1, any)
	})
	counts = counts[apply(genes.use, 1, all),]
	print(dim(counts))
    }

    # Make seurat object
    seur = make_seurat(name='merge', dge=counts, minc=0, maxc=1e9, ming=0, maxg=1e9, ident_fxn=ident_fxn, verbose=T)
    seur@ident = ident.use[colnames(seur@data)]
    seur@data.info$organism = org.use[colnames(seur@data)]
    seur@data.info$dataset = dset.use[colnames(seur@data)]
        
    seur
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
    new_ident[i,1] = as.character(old_ident[i,1])
    new_ident = structure(as.factor(new_ident[,1]), names=rownames(new_ident))
    return(new_ident)
}

fast_ident = function(seur, ident_map, partial=F){
    ident_map = data.frame(stack(ident_map), row.names=1)
    ident_map = ident_map[levels(seur@ident),,drop=F]
    u = as.character(ident_map[,1])
    ident_map[,1] = ave(u, u, FUN=function(x){if(length(x) > 1){paste(x, 1:length(x), sep='_')} else {x}})
    ident = seur@ident
    if(partial == FALSE){
        levels(ident) = ident_map[as.character(levels(ident)),1]
    } else {
        i = as.character(levels(ident)) %in% rownames(ident_map)
	levels(ident)[i] = ident_map[as.character(levels(ident))[i], 1]
    }
    return(ident)
}

update_signatures = function(seur){
    
    # Predict host (human or mouse)
    genes = rownames(seur@data)
    
    # Calculate signatures
    seur@data.info = cbind(seur@data.info, score_cells(seur, files='cell_cycle'))
    seur@data.info = cbind(seur@data.info, score_cells(seur, files='early'))
    seur@data.info$nGene = colSums(seur@data > 0)
    
    return(seur)
}

hclust_ident = function(seur=NULL, data.use=NULL, ident.use=NULL, genes.use=NULL, agg='after', dmethod='euclidean', hmethod='complete'){

    # Input: data.use (features x cells) and ident.use (1 x cells), or seurat object
    # Output: hclust object
    
    # Setup input data
    if(is.null(data.use)){
        if(is.null(genes.use)){
	    genes.use = seur@var.genes
	}
	genes.use = intersect(rownames(seur@data), genes.use)
	data.use = seur@data[genes.use,]
    }
    if(is.null(ident.use)){
        ident.use = seur@ident
    }
    
    # hclust on distance matrix
    if(agg == 'before'){print('Aggregating data')
        print(dim(data.use))
	print(length(ident.use))
        data.use = t(data.frame(aggregate(t(as.matrix(data.use)), list(ident.use), mean), row.names=1))
    }
    print('Calculating distance matrix')
    d = dist(as.data.frame(t(as.matrix(data.use))), method=dmethod)
    if(agg == 'after'){print('Aggregating distances')
        d = data.frame(aggregate(as.matrix(d), list(ident.use), mean), row.names=1)
        d = data.frame(aggregate(t(d), list(ident.use), mean), row.names=1)
    }
    print('Running hclust')
    hclust(as.dist(d), method=hmethod)
}


read_tome_vector = function(tome, name) {
    as.vector(unlist(rhdf5::h5read(tome, name)))
}


read_h5ad_dgCMatrix = function(h5ad, target = "/raw.X") {

    library(Matrix)
    library(rhdf5)
    
    root <- rhdf5::H5Fopen(h5ad)
    
    i_path <- paste0(target,"/indices")
    p_path <- paste0(target,"/indptr")
    x_path <- paste0(target,"/data")
    
    print("Reading indices")
    i <- read_tome_vector(root, i_path)
    print("Reading pointers")
    p <- read_tome_vector(root, p_path)
    print("Reading values")
    x <- read_tome_vector(root, x_path)
    print("Reading observations")
    o <- as.vector(rhdf5::h5read(root, "/obs/_index"))
    print("Reading variables")
    v <- as.vector(rhdf5::h5read(root, "/var/_index"))
    
    print("Reading dimensions")
    dims <- c(length(v), length(o))
    
    H5Fclose(root)
    
    print("Assembling dgCMatrix")
    m <- Matrix::sparseMatrix(i = i,
                              p = p,
                              x = x,
			      dims = dims,
                              index1 = FALSE)
    
    rownames(m) <- v
    colnames(m) <- o
    
    return(m)
}


read_h5ad = function(fn, seur=TRUE, counts_regex='layer.*count'){
    
    library(Matrix)
    library(rhdf5)
    
    # read metadata
    print('reading metadata')
    meta = h5read(fn, '/obs')
    cats = h5read(fn, '/obs/__categories')
    meta = sapply(names(meta), function(a){
        tryCatch({
        if(a %in% names(cats)){
	    b = factor(as.integer(meta[[a]]))
	    levels(b) = cats[[a]]
	} else {
	    b = meta[[a]]
	}
	as.character(b)
	}, error=function(e) {NULL})
    })
    meta = meta[names(meta) != '__categories']
    meta = as.data.frame(do.call(cbind, meta))
    rownames(meta) = h5read(fn, '/obs/_index')
    print(dim(meta))
    
    # read coordinates
    print('reading tsne coords from "/obsm/X_umap"')
    coords = tryCatch({t(h5read(fn, '/obsm/X_umap'))}, error=function(e){NULL})
    
    # read counts
    group = grep(counts_regex, sort(unique(h5ls(fn)$group)), value=T)
    print('reading counts data from h5 index:')
    print(group)
    if(length(group) == 1){
        counts = read_h5ad_dgCMatrix(fn, target=group)
    } else {
        stop('error: found multiple "count" layers')
    }
        
    if(seur == FALSE){
        return(list(counts=counts, meta=meta, coords=coords))
    } else {
        
        # make seurat object
        seur = make_seurat(name='test', dge=counts, ming=0, minc=0, maxg=1e10, ident_fxn=function(a){'all'})

	# add metadata
	seur@data.info = cbind(seur@data.info, meta)

	# add tsne coordinates
	if(!is.null(coords)){
	    seur@tsne.rot = as.data.frame(coords)
	    rownames(seur@tsne.rot) = colnames(seur@data)
	    colnames(seur@tsne.rot) = c('tSNE_1', 'tSNE_2')
	}
	return(seur)
    }
}


