

smart_shuffle = function(data, max_tries=10, max_iter=100){

    # shuffle genes while avoiding duplicate ident/gene pairs
    # data = matrix(1=ident, 2=gene, 3=health, 4=de)
    # returns data with shuffled 'ident' column
    
    for(try_count in 1:max_tries){
        data$ident = sample(data$ident)
        for(iter_count in 1:max_iter){
            ident = data$ident
	    i = duplicated(data[,.(gene, ident)])
	    if(sum(i) == 0){return(data)}
	    j = sample(which(!i), sum(i))
	    ident[i] = data$ident[j]
	    ident[j] = data$ident[i]
	    data$ident = ident
        }
    }
    return(data)
}


build_ix_network = function(data, ix=NULL, max_n=Inf, max_ix=Inf, permute=FALSE, unique=TRUE, symm=TRUE, method='psum'){
    require(data.table)
    
    # build cell-cell interaction network from markers and gene pairs
    #
    # input arguments:
    # - data = matrix(1=ident, 2=gene, 3=health, 4=de)
    # - ix = matrix(1=gene, 2=gene)
    # - max_n = maximum number of cells per gene
    # - max_ix = maxiumum number of interactions per pair
    # - permute = permute gene list?
    # - unique = count unique ligands/receptors (see edgelist2graph)
    # - symm = make graph symmetric by: x = x + t(x)
    
    # process input data
    data = as.data.table(sapply(data, as.character))
    if(any(colnames(data) != c('gene', 'ident', 'health', 'de'))){Rstop('data: check colnames')}
    data$health = as.logical(data$health)
    data$de = as.logical(data$de)

    # interaction list
    if(is.null(ix)){
        ix = read.table('~/aviv/db/fantom/PairsLigRec.txt', sep='\t', header=T, stringsAsFactors=F, comment.char='', quote='')
	ix = ix[ix[,'Pair.Evidence'] == 'literature supported',]
	ix = ix[,c('Ligand.ApprovedSymbol', 'Receptor.ApprovedSymbol')]
	colnames(ix) = c('lig', 'rec')
    }
    ix = as.data.table(sapply(ix, as.character))    
    if(any(colnames(ix) != c('lig', 'rec'))){stop('ix: check colnames')}
    ix = ix[lig %in% data$gene & rec %in% data$gene]
    ix.unique = sort(unique(c(ix$lig, ix$rec)))
    
    # subset data
    data = data[gene %in% ix.unique]
        
    # permute genes
    if(permute == TRUE){
        # permute ligands and receptors
        i = data$gene %in% ix$lig
	data[i] = smart_shuffle(data[i])
	i = data$gene %in% ix$rec
	data[i] = smart_shuffle(data[i])
    }
    
    # initialize network
    nodes = sort(unique(data$ident))
    edges = c()
    graph = matrix(0, nrow=length(nodes), ncol=length(nodes))
    rownames(graph) = colnames(graph) = nodes
    h = list(edges=edges, graph=graph)
    d = list(edges=edges, graph=graph)
    
    # iterate over pairs
    for(i in 1:nrow(ix)){
        
        # ligands and receptors
	l = ix[i, lig]
	r = ix[i, rec]
	
	# healthy edges
	h.edges = c()
	h.l = data[gene == l & health == TRUE, ident]
	h.r = data[gene == r & health == TRUE, ident]
	if(length(h.l) > 0 & length(h.r) > 0){
	    h.edges = expand.grid(h.l, h.r)
	}
	
	# de edges
	d.edges = c()
	d.l = data[gene == l & de == TRUE, ident]
	D.l = unique(h.l, d.l)
	d.r = data[gene == r & de == TRUE, ident]
	D.r = unique(h.r, d.r)
	if(length(d.l) > 0 & length(D.r) > 0){
	    d.edges = rbind(d.edges, as.matrix(expand.grid(d.l, D.r)))
	}
	if(length(D.l) > 0 & length(d.l) > 0){
	    d.edges = rbind(d.edges, as.matrix(expand.grid(D.l, d.r)))
	}
	d.edges = unique(d.edges)
	
	# test & add to network
	if(length(h.edges) > 0){
	    if(length(h.l) <= max_n & length(h.r) <= max_n){
	        if(nrow(h.edges) <= max_ix){
	            h.edges = cbind(h.edges, l, r)
		    colnames(h.edges) = c('lcell', 'rcell', 'lig', 'rec')
	            h$edges = rbind(h$edges, h.edges)
		}
	    }
	}
	
	if(length(d.edges) > 0){
	    if(length(d.l) <= max_n & length(d.r) <= max_n){
	        if(nrow(d.edges) <= max_ix){
	            d.edges = cbind(d.edges, l, r)
		    colnames(d.edges) = c('lcell', 'rcell', 'lig', 'rec')
	            d$edges = rbind(d$edges, d.edges)
		}
	    }
	}
    }
    
    # build graph
    if(!is.null(h$edges)){
        h$edges = as.data.table(h$edges)
	h$graph = edgelist2graph(nodes=nodes, edges=h$edges, unique=unique, symm=symm, method=method)
    }
    if(!is.null(d$edges)){
        d$edges = as.data.table(d$edges)
        d$graph = edgelist2graph(nodes=nodes, edges=d$edges, unique=unique, symm=symm, method=method)
    }
    
    return(list(h=h,d=d))
}


edgelist2graph = function(nodes, edges, unique=TRUE, symm=TRUE, method='psum'){
    require(data.table)
    
    # converts an edgelist to an adjacency matrix
    # edgelist = (1=lcell, 2=rcell, 3=lig, 4=rec)
    # if unique == FALSE, count all ligand/receptor pairs
    # if unique == TRUE, count unique ligands/receptors and take minimum for each pair
    
    # initialize graph
    if(! all(c('lcell', 'rcell', 'lig', 'rec') %in% colnames(edges))){stop('edges: check colnames')}
    edges = as.data.table(edges)
    graph = matrix(0, nrow=length(nodes), ncol=length(nodes))
    rownames(graph) = colnames(graph) = nodes
    
    if(unique == FALSE){
        g = as.matrix(table(edges[,.(lcell, rcell)]))
    } else {
        u = table(unique(edges[,.(lcell, rcell, lig)])[,.(lcell, rcell)])
       	v = table(unique(edges[,.(lcell, rcell, rec)])[,.(lcell, rcell)])
	if(method == 'pmin'){g = as.matrix(pmin(u,v))}
	if(method == 'pmax'){g = as.matrix(pmax(u,v))}
	if(method == 'psum'){g = as.matrix(u + v)}
    }
    
    # make symmetric
    graph[rownames(g), colnames(g)] = g
    if(symm == TRUE){
        graph = graph + t(graph)
    }
    
    return(graph)
}
