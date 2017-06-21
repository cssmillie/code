library(ggplot2)
library(gtools)
source('~/code/single_cell/multiplot.r')
source('~/code/single_cell/map_gene.r')

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

material.heat <- function(n)
{
    mh = c(
        #"#607D8B", #blue grey
        "#283593", #indigo 800
        "#3F51B5",  #indigo
        "#2196F3", # blue
        #"#03A9F4", # light blue
        "#00BCD4", #cyan
        #"#009688", # teal
        "#4CAF50", #green
        "#8BC34A", #light green
        "#CDDC39", # lime
        "#FFEB3B", #yellow
        "#FFC107", # amber
        "#FF9800", # organe
        "#FF5722", # deep orange)
        "#f44336")
    colorRampPalette(mh)(n)
}

load_signature = function(file=NULL){

    if(!file.exists(file)){file = paste0('~/aviv/db/markers/', file, '.txt')}
    sig = read.table(file, stringsAsFactors=F, row.names=1)
    sig = structure(strsplit(sig[,1], ','), names=rownames(sig))
    return(sig)
}

get_scores = function(seur, genes=NULL, file=NULL, top=NULL, combine='mean'){
    
    # Gene signature
    if(is.null(genes)){genes = load_signature(file)}
    if(!is.list(genes)){genes = structure(genes, names=genes)}

    # Organism name
    org = predict_organism(rownames(seur@data)[1:100])
    genes = lapply(genes, map_gene, target=org)
    
    # Select top genes
    if(!is.null(top)){genes = lapply(genes, function(a){as.character(na.omit(a[1:top]))})}
    
    # Calculate TPM
    data = tpm(seur, genes.use=unique(na.omit(as.character(unlist(genes)))))
    
    # Score modules
    scores = sapply(genes, function(a){
        i = intersect(rownames(data), a)
	colMeans(data[i,,drop=F], na.rm=T)
    })

    # Log transform
    scores = log2(scores + 1)
    return(scores)
}

plot_genes_tsne = function(seur, genes=NULL, file=NULL, top=NULL, combine='mean', ...){
    scores = get_scores(seur, genes=genes, file=file, top=top, combine=combine)
    plot_tsne(seur, scores, ...)
}

plot_pcs = function(seur, pcs=NULL, ...){
    if(is.null(pcs)){pcs = 1:seur@data.info$num_pcs[[1]]}
    pcs = seur@pca.rot[,pcs]
    plot_tsne(seur, pcs, ...)
}

plot_clusters = function(seur, ...){
    clusters = sort(grep('Cluster', colnames(seur@data.info), value=T))
    clusters = sapply(seur@data.info[,clusters], function(a){drop.levels(as.factor(a))})
    plot_tsne(seur, clusters, ...)
}

plot_stats = function(seur, ...){
    for(col in c('nGene', 'nUMI', 'G1S', 'G2M')){
        if(!col %in% colnames(seur@data.info)){
	    seur@data.info[,col] = get_scores(seur, col)
	}
    }
    quality = subset(seur@data.info, select=c(nUMI, nGene, G1S, G2M))
    plot_tsne(seur, quality, ...)
}

plot_tsne = function(seur, scores=NULL, ident=TRUE, cells.use=NULL, ymin=0, ymax=1, num_col='auto', base_size=12, label.size=6, do.label=T, do.title=TRUE, ...){
    
    # Plot columns of scores on TSNE
    
    # Dataframe for ggplot
    d = structure(seur@tsne.rot[,1:2], names=c('x', 'y'))
    ps = list()
    
    # Cell identities
    if(!is.logical(ident)){
        d$Identity = ident
    } else if(ident){
        d$Identity = seur@ident
    }
    
    # Cell scores
    if(!is.null(scores)){d = cbind(d, scores)}
    
    # Subset cells
    if(!is.null(cells.use)){d = d[cells.use,]}
    
    for(col in setdiff(colnames(d), c('x', 'y'))){
	
	if(is.numeric(d[,col])){
	    
	    # Continuous plot
	    u = quantile(d[,col], ymin, na.rm=T)
	    v = quantile(d[,col], ymax, na.rm=T)
	    d[,col][d[,col] < u] = u
	    d[,col][d[,col] > v] = v
	    p = ggplot(d) +
	        geom_point(aes_string(x='x',y='y',colour=col), ...) +
		scale_colour_gradientn(colours=material.heat(50), guide=guide_colourbar(barwidth=.5, title='log2(TPM)')) + 
		theme_minimal(base_size=base_size) +
		xlab('TSNE 1') + ylab('TSNE 2')

	} else {
	    
	    # Discrete plot
	    colors = material.heat(nlevels(d[,col]))
	    
	    p = ggplot(d) +
	        geom_point(aes_string(x='x',y='y',colour=col), ...) +
		theme_minimal(base_size=base_size) +
		xlab('TSNE 1') + ylab('TSNE 2') +
		scale_colour_brewer(palette='Set2') + 
		theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	    
	    if(do.label == T){
	        t = aggregate(d[,c('x', 'y')], list(d[,col]), median)
		colnames(t) = c('l', 'x', 'y')
		p = p + geom_text(data=t, aes(x=x, y=y, label=l), size=label.size) + theme(legend.position='none')
	    }
	}
	if(do.title == TRUE){
	    p = p + ggtitle(col)
	}
	ps[[col]] = p
    }
    
    if(num_col == 'auto'){
        num_col = ceiling(sqrt(length(ps)))
    }
    multiplot(plotlist=ps, cols=num_col)
}


plot_heatmap = function(seur, data.use='raw', genes.use=NULL, cells.use=NULL, group_by=NULL, scale='row', ...){

    library(NMF)
    # Plot gene expression with aheatmap
    
    # Get data to plot
    if(data.use == 'raw'){
        data = 10000*scale(seur@raw.data, center=F, scale=colSums(seur@raw.data))
    } else if(data.use == 'data'){
        data = seur@data
    } else if(data.use == 'scale'){
        data = seur@scale.data
    }
    
    # Subset genes and cells
    if(is.null(genes.use)){genes.use = rownames(seur@data)}
    if(is.null(cells.use)){cells.use = colnames(seur@data)}
    data = data[genes.use, cells.use]
    
    # Average within each group
    if(!is.null(group_by)){
        data = t(data.frame(aggregate(t(data), list(as.character(group_by)), mean), row.names=1))
    }
    
    # Plot heatmap
    aheatmap(data, scale='row', border_color='black', ...)
}
