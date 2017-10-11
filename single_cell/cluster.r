require(Seurat)
require(proxy)

source('~/code/single_cell/parallel.r')

cosine_dist = function(x){as.dist(1 - crossprod(x)/crossprod(t(matrix(sqrt(colSums(x**2))))))}

run_graph_cluster = function(data, k=100,  method='infomap', weighted=FALSE, dist='cosine', do.fast=FALSE, out=NULL, knn=NULL){

    # Graph cluster rows of data
    require(cccd)

    if(is.null(knn)){
        print('Building kNN graph')
        knn = nng(data, k=k, method=dist, use.fnn=do.fast)
        if(weighted == TRUE){    
        
	    print('Calculating Jaccard similarity')
            s = similarity(knn, method='jaccard')
            
	    print('Building weighted graph')
            knn = graph.adjacency(s, mode='undirected', weighted=T)
        
        }
    }
    if(method == 'louvain'){
        print('Louvain clustering')
        m = cluster_louvain(as.undirected(knn))
    }

    if(method == 'infomap'){
        print('Infomap clustering')
        m = cluster_infomap(knn)
    }
    
    # Write output
    print(sprintf('Clustering with k = %d finished', k))
    if(!is.null(out)){
        write.table(m$membership, file=out, sep='\t', quote=F)
    }
    
    # Return clusters
    clusters = data.frame(x=as.character(m$membership), row.names=rownames(data), stringsAsFactors=F)
    colnames(clusters) = c(paste0(method, '.', dist, '.k', k))
    return(clusters)
}

run_phenograph = function(data, k=50, dist='cosine', out=NULL){

    # Phenograph cluster rows of data
    
    # Write data
    if(is.null(out)){
	out = tempfile(pattern='phenograph.', tmpdir='~/tmp', fileext='.txt')
    }
    write.table(data, file=out, sep='\t', quote=F)
    
    # Run phenograph
    system(paste0('python ~/code/single_cell/run_phenograph.py --data ', out, ' -k ', k, ' --metric ', dist, ' --out ', out))
    
    # Cleanup files
    clusters = readLines(out)
    if(is.null(out)){
        system(paste0('rm ', out))
    }
    
    # Return clusters
    clusters = data.frame(x=as.character(clusters), row.names=rownames(data), stringsAsFactors=F)
    colnames(clusters) = c(paste0('phenograph.', dist, '.k', k))
    return(clusters)
}

run_mcl = function(data, method='spearman'){

    # MCL cluster rows of data
    require(MCL)
	
    # Calculate similarity
    s = cor(data, method=method)
    
    # Cluster data
    m = mcl(s, addLoops=T)

    # Return clusters
    return(m$Cluster)
}

run_density = function(data, dist='cosine'){

    # Density cluster rows of data
    require(densityClust)    
    if(dist == 'cosine'){
        d = cosine_dist(t(data))
    } else {
        d = dist(data, method=dist)
    }
    p = densityClust(d)
    q = findClusters(p)
    clusters = data.frame(x=as.character(q$clusters), row.names=rownames(data), stringsAsFactors=F)
    colnames(clusters) = c(paste('density.', dist, '.k', k))
    return(clusters)
}

run_cluster = function(data, k, method='infomap', weighted=FALSE, n.cores=1, dist='cosine', do.fast=TRUE, prefix=NULL, knn=NULL){

    # Cluster rows of data
    
    g = run_parallel(
        foreach(i=k, .packages=c('cccd', 'MCL'), .export=c('run_graph_cluster', 'run_phenograph', 'run_mcl'), .combine=cbind) %dopar% {
	    
	    # Get output file
	    if(!is.null(prefix)){
	        out = paste0(prefix, '.', i, '.', method, '.membership.txt')
	    } else{
	        out = NULL
	    }
	    
	    # Get clusters
	    if(method %in% c('infomap', 'louvain')){
            	gi = run_graph_cluster(data, k=i, method=method, weighted=weighted, dist=dist, do.fast=do.fast, out=out, knn=knn)
	    }
	    if(method == 'phenograph'){
	        gi = run_phenograph(data, k=i, dist=dist, out=out)
	    }
	    if(method == 'mcl'){
	        gi = run_mcl(data, method='spearman')
	    }
	    return(gi)
        },
	n.cores = n.cores
    )
    return(g)
}


