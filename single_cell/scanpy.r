run_scanpy = function(data, n_pcs=0, method='dmap', dpt=FALSE, root_cell=NULL, n_neighbors=10, paga=FALSE, clusters=NULL, cleanup=TRUE, out=NULL, sparse=TRUE){

    # -----------------------------------------
    # Run scanpy on [cells x feats] data matrix
    # -----------------------------------------
    # For dmap, use t(seur@data) or seur@pca.rot
    # For dca, use t(seur@raw.data)
    # method should be 'dmap' or 'dca'
    
    # Filter data (remove zero genes)
    data = data[, colSums(data > 0) >= 3]
    
    # Write data
    if(is.null(out)){
        out = tempfile(pattern='scanpy.', tmpdir='~/tmp', fileext='.data.txt')
    } else {
        out = paste0(out, '.data.txt')
    }
    if(sparse == FALSE){
        write.table(as.matrix(data), file=out, sep='\t', quote=F, row.names=FALSE, col.names=FALSE)
	print(paste0('Writing dense matrix to ', out))
    } else {
        out = write_mtx(data, prefix=gsub('.data.txt', '', out))$data
	print(paste0('Writing sparse matrix to ', out))	
    }
    
    # Write clusters
    if(is.null(clusters)){clusters = 'louvain_groups'} else {
        clusters_fn = gsub('data', 'clusters', out)
	write.table(clusters, clusters_fn, quote=F, row.names=FALSE, col.names=FALSE)
    }
    
    # Root cell index
    if(is.null(root_cell)){iroot = 0} else {iroot = match(root_cell, rownames(data))-1}
    
    if(method == 'dmap'){
    
        print('Calculating diffusion map. Make sure data = t(seur@data) or seur@pca.rot')
        
    	# Diffusion map
    	command = paste('python ~/code/single_cell/run_scanpy.py --data', out, '--dmap', '--n_pcs', n_pcs, '--n_neighbors', n_neighbors, '--iroot', iroot, '--out', gsub('.data.*', '', out))
        
    	# Pseuodotime
    	if(dpt == TRUE){command = paste(command, '--dpt')}
	
    	# Approximate graph abstraction
        if(paga == TRUE){command = paste(command, '--paga', '--clusters', clusters_fn)}
	
    }
    if(method == 'dca'){

        print('Imputing data with DCA. Make sure data = t(seur@raw.data)')
    	
        # deep count autoencoder
	command = paste('python ~/code/single_cell/run_scanpy.py --data', out, '--dca', '--out', gsub('.data.*', '', out))
    }
    
    system(command)
    
    # Load & cleanup results
    res = list(dmap='X_diffmap', dpt_time='dpt_pseudotime', paga_full='connectivities', paga_tree='connectivities_tree', categories='categories', dca='dca')
    
    for(name in names(res)){
        fn = gsub('data', res[[name]], out)
	if(file.exists(fn)){
	    res[[name]] = as.data.frame(fread(fn, header=FALSE))
	    if(nrow(data) == nrow(res[[name]])){
	        rownames(res[[name]]) = rownames(data)
	    }
	    if(ncol(data) == ncol(res[[name]])){
	        colnames(res[[name]]) = colnames(data)
	    }
	    if(cleanup == TRUE){system(paste('rm', fn))}
	}
    }
    tryCatch({for(name in c('paga_tree', 'paga_full')){rownames(res[[name]]) = colnames(res[[name]]) = res[['categories']][,1]}}, error=function(e){})
    if(cleanup == TRUE){system(paste('rm', out))}
    
    return(res)
}
