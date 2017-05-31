library(Seurat)
library(proxy)

source('~/code/single_cell/parallel.r')

cosine_dist = function(x){as.dist(1 - crossprod(x)/crossprod(t(matrix(sqrt(colSums(x**2))))))}

run_graph_cluster = function(d, k=100,  method='infomap', weighted=FALSE, dist='cosine', do.fast=FALSE, out=NULL){
    
    library(cccd)
    
    print('Building kNN graph')
    g = nng(d, k=k, method=dist, use.fnn=do.fast)
    if(weighted == TRUE){    
    
	print('Calculating Jaccard similarity')
        s = similarity(g, method='jaccard')
    
	print('Building weighted graph')
        g = graph.adjacency(s, mode='undirected', weighted=T)
    
    }
    if(method == 'louvain'){
        print('Louvain clustering')
        m = cluster_louvain(as.undirected(g))
    }

    if(method == 'infomap'){
        print('Infomap clustering')
        m = cluster_infomap(g)
    }
    
    # Write output
    print(sprintf('Clustering with k = %d finished', k))
    if(!is.null(out)){
        write.table(m$membership, file=out, sep='\t', quote=F)
    }
    
    # Return clusters
    clusters = data.frame(x=m$membership, row.names=rownames(data), stringsAsFactors=F)
    colnames(clusters) = c(paste0(method, '.k', k))
    return(clusters)
}

run_phenograph = function(data, k=50, dist='cosine', out=NULL){
    
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
    clusters = data.frame(x=clusters, row.names=rownames(data), stringsAsFactors=F)
    colnames(clusters) = c(paste0('phenograph.k', k))
    return(clusters)
}

run_mcl = function(data, method='spearman'){

    library(MCL)
	
    # Calculate similarity
    s = cor(data, method=method)
    
    # Cluster data
    m = mcl(s, addLoops=T)

    # Return clusters
    return(m$Cluster)
}

run_density = function(data, dist='correlation'){

    library(densityClust)    
    d = cosine_dist(t(data))
    p = densityClust(d)
    q = findClusters(p)
    return(q$clusters)
    
}

run_cluster = function(data, k, method='infomap', weighted=FALSE, n.cores=1, dist='correlation', do.fast=TRUE, prefix=NULL){
    
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
            	gi = run_graph_cluster(data, k=i, method=method, weighted=weighted, dist=dist, do.fast=do.fast, out=out)
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


