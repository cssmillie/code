
expand.grid.df = function(...) Reduce(function(...) merge(..., by=NULL), list(...))

smart_shuffle = function(data, cols.shuf, cols.test=NULL, max_tries=10, max_iter=100){
    
    # Shuffle columns while avoiding duplicate pairs
    # data = [m x n] matrix
    # cols.shuf = list of columns to shuffle
    # cols.test = list of columns to test

    # Fix input data
    data = as.data.frame(data)
    if(is.null(cols.test)){cols.test = colnames(data)}
    if(! all(cols.shuf %in% cols.test)){stop('error: cols.shuf not in cols.test')}

    # For each attempt, shuffle columns
    for(a in 1:max_tries){
	data[, cols.shuf] = data[sample(1:nrow(data)), cols.shuf]
	
	# Iteratively fix duplicates
	for(b in 1:max_iter){
	    new_data = data[, cols.shuf, drop=F]
	    i = duplicated(data[,cols.test])
	    if(sum(i) == 0){return(as.data.table(data))}
	    j = sample(which(!i), sum(i))
	    new_data[i,] = data[j, cols.shuf]
	    new_data[j,] = data[i, cols.shuf]
	    data[, cols.shuf] = new_data
	}
    }
    as.data.table(data)
}

make_ix_network = function(data, diff=NULL, ix=NULL, weights=NULL, permute=FALSE, perm.col='gene', symm=TRUE, method='sum', do.intersect=FALSE){
    require(data.table)
    
    # Build cell-cell interaction network from markers and gene pairs
    #
    # Input arguments:
    # - data = matrix(1=ident, 2=gene) = cell type markers
    # - diff = matrix(1=ident, 2=gene) = DE genes (optional)
    # - ix = matrix(1=gene, 2=gene) = interaction matrix
    # - weights = matrix(1=ident, 2=gene, 3=weight)    
    # - permute = permute gene list?
    # - method sum = add total unique edge weights
    # - symm = make graph symmetric by: x = x + t(x)
    # - do.intersect = intersect diff with data, i.e.
    #   require DE genes to also be cell type markers
    # - returns list(graph, edges) for data or diff
    
    # Fix inputs
    data = as.data.table(data[,c('ident', 'gene')])
    data[, ct := 1]
    
    # Merge markers with DE genes
    if(!is.null(diff)){
        diff = as.data.table(diff[,c('ident', 'gene')])
	diff[, de := 1]
	data = merge(data, diff, by=c('ident', 'gene'), all=TRUE)
    } else {
        data[, de := 0]
    }
    data = setkeyv(data, c('ident', 'gene'))
    
    # Remove DE genes that are not cell type markers
    if(do.intersect == TRUE){data = data[ct == 1]}
    
    # Add weights to data
    if(!is.null(weights)){
        weights = as.data.table(weights[,c('ident', 'gene', 'w')])
	data = merge(data, weights, by=c('ident', 'gene'), all.x=TRUE)
    } else {
        data[, w := 1]
    }
    print(head(data))
    
    # Read interactions
    if(is.null(ix)){
        ix = read.table('~/aviv/db/fantom/PairsLigRec.txt', sep='\t', header=T, stringsAsFactors=F, comment.char='', quote='')
	ix = ix[ix[,'Pair.Evidence'] == 'literature supported',]
	ix = ix[,c('Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol')]
	colnames(ix) = c('lig', 'rec')
    }
    ix = as.data.table(sapply(ix, as.character))
    colnames(ix) = c('lig', 'rec')
    
    # Intersect genes
    if(is.null(diff)){
        ix = ix[lig %in% data$gene & rec %in% data$gene]
    } else {
        ix = ix[(lig %in% data$gene & rec %in% diff$gene) | (rec %in% data$gene & lig %in% diff$gene)]
    }
    genes.use = sort(unique(c(ix$lig, ix$rec)))
    data = data[gene %in% genes.use]
    
    # Permute ligands and receptors
    if(permute == TRUE){
        i = data$gene %in% ix$lig
	data[i] = smart_shuffle(data[i], perm.col, c('ident', 'gene'))
	i = data$gene %in% ix$rec
	data[i] = smart_shuffle(data[i], perm.col, c('ident', 'gene'))
    }
    
    # Initialize network
    nodes = sort(unique(data$ident))
    edges = c()
    graph = matrix(0, nrow=length(nodes), ncol=length(nodes))
    rownames(graph) = colnames(graph) = nodes
    n = list(edges=edges, graph=graph)
    d = list(edges=edges, graph=graph)
    
    # Iterate over ligand-receptor pairs
    for(i in 1:nrow(ix)){
        
        # Ligands and receptors
	l = ix[i, lig]
	r = ix[i, rec]
	
	# Baseline edges
	n.edges = c()
	n.l = data[gene == l & ct == TRUE, ident]
	n.r = data[gene == r & ct == TRUE, ident]
	if(length(n.l) > 0 & length(n.r) > 0){
	    n.edges = as.matrix(expand.grid(n.l, n.r))
	}
	
	# DE edges
	d.edges = c()
	d.l = data[gene == l & de == TRUE, ident]
	D.l = unique(c(n.l, d.l))
	d.r = data[gene == r & de == TRUE, ident]
	D.r = unique(c(n.r, d.r))
	
	if(length(d.l) > 0 & length(D.r) > 0){
	    d.edges = rbind(d.edges, as.matrix(expand.grid(d.l, D.r)))
	}
	if(length(D.l) > 0 & length(d.l) > 0){
	    d.edges = rbind(d.edges, as.matrix(expand.grid(D.l, d.r)))
	}
	d.edges = unique(d.edges)
	
	# Fix names
	if(!is.null(nrow(n.edges))){if(nrow(n.edges) > 0){
	    n.edges = cbind(n.edges, l, r)
	    colnames(n.edges) = c('lcell', 'rcell', 'lig', 'rec')
	    n$edges = rbind(n$edges, n.edges)
	}}
	if(!is.null(nrow(d.edges))){if(nrow(d.edges) > 0){
	    d.edges = cbind(d.edges, l, r)
	    colnames(d.edges) = c('lcell', 'rcell', 'lig', 'rec')
	    d$edges = rbind(d$edges, d.edges)
	}}
    }
    
    # Get weights from data
    weights = data[,.(ident, gene, w)]
    
    # Build graph from edgelist
    if(is.null(diff)){
        edges = as.data.table(n$edges)
    } else {
        edges = as.data.table(d$edges)
    }
    graph = edgelist2graph(nodes=nodes, edges=edges, weights=weights, symm=symm, method=method)
    
    # Remove edges from permuted data
    if(permute == TRUE){edges = NULL}
    return(list(edges=edges, graph=graph))
}

edgelist2graph = function(nodes, edges, weights=NULL, symm=TRUE, method='sum'){
    require(data.table)
    
    # Converts an edgelist to an adjacency matrix
    # edgelist = list(1=lcell, 2=rcell, 3=lig, 4=rec)
    # weights = matrix(1=ident, 2=gene, 3=weight)
    # symm = make graph symmetric by: x = x + t(x)
    # method sum = add total unique edge weights
    
    # Initialize adjacency matrix
    edges = as.data.table(edges)
    graph = as.data.frame(matrix(0, nrow=length(nodes), ncol=length(nodes)))
    rownames(graph) = colnames(graph) = nodes
    
    # Get unique ligands and receptors
    lig = unique(edges[,.(lcell, rcell, lig)])
    rec = unique(edges[,.(lcell, rcell, rec)])
    
    # Add edge weights
    if(!is.null(weights)){
        setkeyv(weights, c('ident', 'gene'))
	lig[, w := weights[.(lcell, lig)]$w]
	rec[, w := weights[.(rcell, rec)]$w]
    }
    
    # Aggregate edge weights
    lig = lig[, .(w=sum(w)), .(lcell, rcell)]
    rec = rec[, .(w=sum(w)), .(lcell, rcell)]
        
    # Convert to matrix
    lig = data.frame(spread(lig, rcell, w), row.names=1)
    rec = data.frame(spread(rec, rcell, w), row.names=1)
    lig[is.na(lig)] = 0
    rec[is.na(rec)] = 0
    
    # Align data and combine
    x = y = graph
    x[rownames(lig), colnames(lig)] = lig
    y[rownames(rec), colnames(rec)] = rec   
    if(method == 'sum'){graph = as.matrix(x + y)}

    # Make symmetric
    if(symm == TRUE){
        graph = graph + t(graph)
    }
    
    return(graph)
}

matrix_pvals = function(graph, shuf, type='gt'){
    
    # Calculate empirical p-values
    # true = [m x n] matrix
    # shuf = list([m x n] matrices)
    # type 'gt' = test true data > shuffled data
    # returns [m x n] matrix of p-values
        
    if(type == 'lt'){
        Reduce('+', sapply(shuf, function(a) graph >= a, simplify=F))/length(shuf)
    }
    if(type == 'gt'){
        Reduce('+', sapply(shuf, function(a) graph <= a, simplify=F))/length(shuf)
    }
}

matrix_mean = function(graphs){Reduce('+', graphs)/length(graphs)}

matrix_qvals = function(graph, shuf, type='gt', n=100, pvals=c(0, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2)){

    true_pvals = matrix_pvals(graph, shuf)
    true_pvals = true_pvals[lower.tri(true_pvals, diag=T)]

    shuf_pvals = as.vector(sapply(1:min(n, length(shuf)), function(i){print(i)
        p = matrix_pvals(shuf[[i]], shuf[setdiff(1:length(shuf), i)])
        p = p[lower.tri(p, diag=T)]
        as.vector(p)
    }))
    num = sapply(pvals, function(p) sum(true_pvals <= p))
    fdr = sapply(pvals, function(p) mean(shuf_pvals <= p)/mean(true_pvals <= p))
    return(list(pvals=pvals, num=num, fdr=fdr))
}

