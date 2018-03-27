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


rescale_vector = function(x, target=c(0, 1), f=identity, abs=TRUE){
    # map vector onto target interval, then apply f
    a = target[[1]]
    b = target[[2]]
    if(min(x) == max(x)){
        rep(max(target), length(x))
    } else {
        if(abs == TRUE){
            x = (x - min(abs(x)))/(max(abs(x)) - min(abs(x)))*(b - a) + a
	} else {
            x = (x - min(x))/(max(x) - min(x))*(b - a) + a
	}
        f(x)
    }
}


plot_network = function(edges, labels=NULL, label_n=25, f_weights=identity, rescale_weights=c(0,1), f_coords=identity, rescale_coords=c(0,1), alpha=1, color_edges=FALSE){
    
    # plots a network diagram from edges
    # data = (nodes x nodes) weights matrix or (node, node, weight) list
    # labels = list, labels[node][node] = label
    
    library(igraph)
    library(tidyverse)
    
    # get edgelist
    if(ncol(edges) != 3){
        edges = gather(as.data.frame(edges) %>% rownames_to_column('source'), target, weight, -source)
    }
    
    # remove zeros and de-duplicate
    edges = edges[edges$weight != 0,]
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
    l = apply(l, 2, function(a) rescale_vector(a, target=rescale_coords, f=f_coords, abs=TRUE))
    
    # rescale edge widths
    w = edges$weight
    w = rescale_vector(w, target=rescale_weights, f=f_weights, abs=TRUE)
    edges$color = ifelse(edges$weight > 0, set.colors[[1]], set.colors[[2]])
    print(edges)
    
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


ggplot_network = function(graph=NULL, edges=NULL, node_sizes=NULL, node_colors=NULL, scale_nodes=c(0,10), scale_edges=c(0,1), alpha=1, symm=TRUE, do.legend=TRUE, legend_title=NULL, out=NULL, nrow=1, ncol=1){
    
    library(ggplot2)
    library(tidyverse)
    library(cowplot)
    library(ggnetwork)
    source('~/code/single_cell/colors.r')

    # plot network from graph or edgelist
    # -----------------------------------

    # graph = [m x n] adjacency graph
    # edges = data.frame(1=source, 2=target, 3=weight)
    # node_sizes = unnamed list (size1, size2, ...) or named list, list(node1=size1, node2=size2, ...)
    # node_colors = unnamed list (color1, color2, ...) or named list, list(node1=color1, node2=color2, ...)
    
    # convert data to long format
    if(!is.null(graph)){edges = gather(as.data.frame(graph) %>% rownames_to_column('source'), target, weight, -source)}
    if(!is.null(node_sizes) & is.null(names(node_sizes))){names(node_sizes) = rownames(graph)}
    if(!is.null(node_colors) & is.null(names(node_colors))){names(node_colors) = rownames(graph)}
    nodes = sort(unique(c(rownames(graph), colnames(graph))))
    if(is.null(node_sizes)){node_sizes = structure(rep(1, length(nodes)), names=nodes)}
    if(is.null(node_colors)){node_colors = structure(rep(1, length(nodes)), names=nodes)}
        
    # remove duplicate edges
    if(symm == TRUE){
        i = apply(edges[,1:2], 1, function(a) paste(sort(a), collapse=' '))
	edges = edges[!duplicated(i),]
    }
    
    # adjust edge weights
    edges = edges[edges$weight > 0,]
    edges = edges[edges$source != edges$target,]
    if(!is.null(scale_edges)){
        edges$weight = rescale_vector(edges$weight, target=scale_edges, abs=TRUE)
    }
    
    # convert edges to igraph object
    g = graph.data.frame(edges)
    
    # vertex attributes / igraph craziness
    nodes = V(g)$name
    
    # node sizes
    if(length(node_sizes) > 1){
        node_sizes = node_sizes[nodes]
        if(!is.null(scale_nodes)){
            node_sizes = rescale_vector(node_sizes, target=scale_nodes)
        }
    }
    V(g)$node_size = node_sizes
    
    # node colors
    if(length(node_colors) > 1){
        node_colors = node_colors[nodes]
    }
    V(g)$node_color = as.character(node_colors)
    
    # edge colors
    if(TRUE){
        E(g)$edge_color = '#cccccc'
    } else {
        E(g)$edge_color = '#cccccc'
    }

    # igraph -> data.frame for ggplot2
    n = ggnetwork(g, weights='weight')
    n$node_color = as.character(n$node_color)
    
    # plot with ggnetwork
    if(is.null(out)){alpha = 1}
    p = ggplot(n, aes(x=x, y=y, xend=xend, yend=yend)) +
	geom_edges(color='#cccccc', size=na.omit(n$weight), alpha=alpha) +
        geom_nodes(aes(color=node_color), size=6) + 
	geom_nodetext_repel(aes(label=vertex.names)) +
	scale_color_manual(legend_title, values=tsne.colors) +
	theme_blank()
    
    if(do.legend == FALSE){p = p + theme(legend.position='none')}

    if(!is.null(out)){save_plot(p, file=out, nrow=nrow, ncol=ncol)}
    
    p
}







