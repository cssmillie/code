seurat()

make_edges = function(data, pairs, subset=NULL, binary=FALSE, symmetric=TRUE, max_cells_per_gene=2, verbose=TRUE){
    
    # this function creates a (cells x cells) edge matrix
    # weights wij are calculated as the product of the node weights, wij=wi*wj
    # input arguments:
    # data = (cells x genes) weights matrix or (cells, genes, weights) list
    # pairs = (gene, gene) list
    # subset = (gene -> cells) list

    library(Matrix)
    library(tidyr)
    library(data.table)
    
    # format data
    print(data)
    if(ncol(data) == 3){
        colnames(data) = c('cell', 'gene', 'value')
	data = data.frame(spread(data, gene, value, fill=0), row.names=1)
	data = as(as.matrix(data), 'sparseMatrix')
    }
    print(dim(data))
    
    # make binary
    if(binary == FALSE & max(data) > 1){
        print('[warning] max(data) > 1')
    }
    
    if(binary == TRUE){
        data = data > 0
    }

    # remove promiscuous genes
    data = data[,colSums(data > 0, na.rm=T) <= max_cells_per_gene]
    print(dim(data))
    
    # make adjacency matrix
    cells = rownames(data)
    A = matrix(0, nrow=length(cells), ncol=length(cells))
    rownames(A) = colnames(A) = cells
    
    # make labels
    L = sapply(cells, function(a) sapply(cells, function(b) c(), simplify=F), simplify=F)
    
    # add edges
    count = 0
    for(i in 1:nrow(pairs)){
        u = as.character(pairs[i,1]) # source gene
	v = as.character(pairs[i,2]) # target gene
	if(u %in% colnames(data) & v %in% colnames(data)){
	    
	    Ai = data[,u] %*% t(data[,v])
	    
	    if(!is.null(subset)){

	        # subset on gene u
		Au = Ai
		i = rownames(data) %in% subset[[u]]
		Au[!i,] = 0
		
		# subset on gene v
		Av = Ai
		i = rownames(data) %in% subset[[v]]
		Av[!i,] = 0
		
		# sum components
		Ai = Au + Av

		# skip on empty
		if(sum(Ai) == 0){
		    next
		}
	    }
	    
	    if(symmetric == TRUE){
	        Ai = Ai + t(Ai)
	    }
	    
	    A = A + Ai
	    
	    # update labels
	    cells1 = rownames(data)[data[,u] > 0]
	    cells2 = rownames(data)[data[,v] > 0]
	    for(ca in cells1){
	        for(cb in cells2){
		    L[[ca]][[cb]] = c(L[[ca]][[cb]], paste(u, v, sep='-'))
		    if(symmetric == TRUE){
		        L[[cb]][[ca]] = c(L[[cb]][[ca]], paste(u, v, sep='-'))
		    }
		}
	    }
	    
	    # print message
	    if(verbose == TRUE){
	        print(paste0('Edge: ', u, ', ', v, ' (', paste(cells1, collapse=', '), ') <--> (', paste(cells2, collapse=', '), ')'))
	    }
	    count = count + 1
        }
    }
    
    # fix labels
    for(ca in names(L)){
        for(cb in names(L[[ca]])){
	    L[[ca]][[cb]] = unique(L[[ca]][[cb]])
	}
    }
    
    print(paste(count, 'total interactions'))
    return(list(edges=A, labels=L))
}


rescale_vector = function(x, target=c(0, 1), f=identity){
    # map vector onto target interval, then apply f
    a = target[[1]]
    b = target[[2]]
    x = (x - min(x))/(max(x) - min(x))*(b - a) + a
    f(x)
}


plot_network = function(edges, labels=NULL, label_n=25, f_weights=identity, rescale_weights=c(0,1), f_coords=identity, rescale_coords=c(0,1), alpha=.5, color_edges=FALSE){
    
    # plots a network diagram from edges
    # data = (nodes x nodes) weights matrix or (node, node, weight) list
    # labels = list, labels[node][node] = label
    
    library(igraph)
    library(tidyverse)
    
    # get edgelist
    if(ncol(edges) != 3){
        edges = gather(as.data.frame(edges) %>% rownames_to_column('source'), target, weight, -source)
    }
    
    # de-duplicate
    i = apply(edges, 1, function(a) paste(sort(c(a[[1]], a[[2]])), collapse=' '))
    edges = edges[!duplicated(i),]
    
    # modify edges
    edges$weight = f_weights(edges$weight)
    
    # delete zero and self edges
    edges = edges[edges$weight > 0,]
    edges = edges[edges$source != edges$target,]
    
    # get nodes
    nodes = data.frame(id=sort(unique(c(edges[,1], edges[,2]))))
    
    # make graph
    g = graph.data.frame(d=edges, vertices=nodes, directed=FALSE)
    l = layout_with_fr(g, weight=edges$weight)
    
    # rescale coordinates
    l = apply(l, 2, function(a) rescale_vector(a, target=rescale_coords, f=f_coords))
    
    # rescale edge widths
    w = edges$weight
    w = rescale_vector(w, target=rescale_weights, f=f_weights)

    # add labels
    if(!is.null(labels)){
        edges$label = apply(edges, 1, function(a) sample(labels[[a[[1]]]][[a[[2]]]], 1))
	edges$label[edges$weight < sort(edges$weight, decreasing=T)[label_n]] = NA	
	edges$color = set.colors[1:nrow(edges)]
    }
    print(edges)
    
    # plot and save
    #plot(g, edge.width=w, edge.color=hsv(0,0,0,alpha=alpha), edge.label=edges$label, vertex.label.color='black', edge.label.color='black', layout=l, vertex.label.family='arial', edge.label.family='arial')
    plot(g, edge.width=w, edge.color=edges$color, edge.label=edges$label, edge.label.color=edges$color, vertex.label.color='black', layout=l, vertex.label.family='arial', edge.label.family='arial', edge.label.cex=1.5, vertex.label.cex=1.5)

}


ggplot_network = function(edges, labels=NULL, label_n=25, f_weights=identity, rescale_weights=NULL, f_coords=identity, rescale_coords=NULL, alpha=.5, color_edges=FALSE, symmetric=FALSE){

    library(ggplot2)
    library(tidyverse)
    library(cowplot)

    # get edgelist
    if(ncol(edges) != 3){
        edges = gather(as.data.frame(edges) %>% rownames_to_column('source'), target, weight, -source)
    }
    
    # de-duplicate
    if(symmetric == TRUE){
        i = apply(edges, 1, function(a) paste(sort(c(a[[1]], a[[2]])), collapse=' '))
	edges = edges[!duplicated(i),]
    }
    
    # modify edges
    edges$weight = f_weights(edges$weight)
    
    # delete zero and self edges
    edges = edges[edges$weight > 0,]
    edges = edges[edges$source != edges$target,]
    
    # rescale edge widths
    if(!is.null(rescale_weights)){
        edges$weight = rescale_vector(edges$weight, target=rescale_weights, f=f_weights)
    }
    
    # add labels
    if(!is.null(labels)){
    	edges$label = apply(edges, 1, function(a) sample(labels[[a[[1]]]][[a[[2]]]], 1))
	edges$label[edges$weight < sort(edges$weight, decreasing=T)[label_n]] = ''
    }
    print(edges)
    # make network
    g = graph.data.frame(edges)
    
    # make ggnetwork
    n = ggnetwork(g, weights='weight')
    print(head(n))
    
    # ggplot
    p = ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color='grey50', size=na.omit(n$weight), alpha=.5) +
    geom_nodes(size=10, color='lightskyblue2') + 
    geom_nodelabel_repel(aes(label=vertex.names), size=3, fill='white', label.r=unit(.05, 'lines'), box.padding=unit(.15, 'lines'), label.padding=unit(.15, 'lines')) +
    theme_blank()
    p
    #geom_nodetext_repel(aes(label=vertex.names), size=3, fill='gold', box.padding=unit(.15, 'lines'), label.padding=unit(.15,'lines'), label.r=unit(0,'lines'))
    #geom_edgetext_repel(aes(label=label), fill=NA, size=2.5, label.r=unit(0, 'lines'))
    #p$data$xmid = (p$data$x + p$data$xend)/2
    #p$data$ymid = (p$data$y + p$data$yend)/2
    #p + geom_label(aes(x=xmid, y=ymid, label=label))
}







