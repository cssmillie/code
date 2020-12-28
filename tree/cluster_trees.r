
seurat()

cluster_trees = function(d=NULL, pre=NULL){
    
    if(is.null(d)){
        d = readRDS(paste0('/home/unix/csmillie/Gut_Human/hyper/', pre, '/tree/trees.qdist.rds'))
    }
    
    # copy distances across diagonal
    if(all(d[lower.tri(d)] == 0)){
        print('adding lower.tri distances')
        d = d + t(d)
    }
    
    # get labels, etc
    lab = tryCatch({sapply(strsplit(rownames(d), '\\.'), '[[', 2)}, error=function(e){rownames(d)})
    
    # umap visualization
    library(umap)
    print('calculating umap')
    v = umap(d, method='umap-learn', n_neighbors=5, n_epochs=2500)$layout
    
    # pre-compute distance matrix (of distances)
    print('distance matrix')
    dx = as.matrix(proxy::dist(d, method=NULL))
    
    # pre-compute nng
    print('nearest neighbor graph')
    k.use = c(2, 3, 5, 10, 15, 20, 25, 50, 100, 150, 200, 250)
    G = sapply(k.use, function(ki){print(ki)
      nng(dx=dx, k=ki, use.fnn=FALSE)
    }, simplify=F)
    names(G) = as.character(k.use)
    
    # infomap clustering
    print('infomap clustering')
    ms = sapply(k.use, function(ki) as.factor(cluster_infomap(as.undirected(G[[as.character(ki)]]))$membership), simplify=F)
    mi = ms[['25']]
    
    return(list(d=d, v=v, lab=lab, dx=dx, G=G, ms=ms, mi=mi))
}

plot_clusters = function(res, mi, lab.n=0, lab.type='random', lab.regex='', do.label=T, lab.use='auto', rename=T, umap=NULL, col.dist=FALSE){
    
    if(rename == TRUE){
        i = ! grepl('group', res$lab)
        names = tapply(substr(res$lab, 1, 3)[i], mi[i], function(a) names(sort(table(a), dec=T))[[1]])
	names[] = paste(names(names), names)
	print(names)
    } else {
        names = setNames(unique(mi), unique(mi))
    }
    
    if(lab.use == 'auto'){
        i = ! grepl('group', res$lab)
	regex = tapply(substr(res$lab, 1, 3)[i], mi[i], function(a) names(sort(table(a), dec=T))[[1]])
	print(regex)
	lab.use = unname(unlist(sapply(regex, function(a) grep(a, res$lab, value=T)[[1]])))
	print(lab.use)
    }

    lab.use = c(lab.use, grep(lab.regex, res$lab, value=T))

    if(is.null(umap)){x = res$v[,1]; y = res$v[,2]} else {x = umap[,1]; y = umap[,2]}

    if(col.dist == TRUE){
        D = collapse_dist(res, mi)
	col = setNames(diag(as.matrix(D)), rownames(D))
        simple_scatter(x, y, lab=res$lab, lab.n=lab.n, lab.type=lab.type, xlab='UMAP 1', ylab='UMAP 2', col=col[as.character(mi)], do.label=do.label, lab.use=lab.use)	
    } else {
        simple_scatter(x, y, lab=res$lab, lab.n=lab.n, lab.type=lab.type, xlab='UMAP 1', ylab='UMAP 2', col=names[as.character(mi)], pal=tsne.colors, do.label=do.label, lab.use=lab.use)
    }
    
}

collapse_dist = function(res, mi){
    u = data.frame(aggregate(res$d, list(mi), mean), row.names=1)
    v = data.frame(aggregate(t(u), list(mi), mean), row.names=1)
    v
}



plot_order = function(cluster, n=1){
gorder = load_signature('~/Gut_Human/hyper/ecoli/ecoli.gene_order.txt')
genes.use = lab[m.use == cluster]
    ecoli.use = names(sort(lengths(sapply(gorder, intersect, genes.use)), decreasing=T))[[n]]
    order.use = gorder[[ecoli.use]]
    simple_scatter(1:length(order.use) %% 500, as.integer(1:length(order.use)/500), lab=order.use, lab.n=0, lab.use=genes.use, pal=c('grey', 'tomato'), max_size=1)
}



test_distance_metrics = function(genes.use, clusters){
    trees = sapply(genes.use, function(a) readRDS(paste0('~/Gut_Human/hyper/Escherichia_coli/tree/RAxML_bestTree.', a)), simplify=FALSE)
    dists = sapply(meths.use, function(a){
        tdist(trees, label_fxn=function(b) gsub('\\..*', '', b), type=a, normalize=TRUE)
    }, simplify=F)
    ratio = sapply(dists, function(di){
        di[is.na(di)] = quantile(di, .75, na.rm=T)
        ri = data.frame(aggregate(di, list(clusters), mean), row.names=1)
	ri = data.frame(aggregate(t(ri), list(clusters), mean), row.names=1)
	ri
    }, simplify=F)
}
