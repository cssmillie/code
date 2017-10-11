
network_diagram = function(data, pairs){

    # data = (cells x genes) expression data (binary or continuous)
    # pairs = (n x 2) interaction list
    
    # initialize edges
    edges = data.frame(NA, nrow=nrow(data), ncol=nrow(data))

    
    
}






network_diagram = function(edges, nodes, out){

    # edges = df(source, target, weight)
    # nodes = df(id, group)

    library(igraph)
    g = graph.data.frame(d=edges, vertices=nodes, directed=FALSE)
    l = layout_with_fr(g, weight=w)
    
    
    plot(g, edge.width=edges$weight, vertex.color=nodes$group, edge.color=hsv(0,0,0,alpha=.5), layout=l)
    
}
