require(ggplot2)
require(gtools)
require(cowplot)

source('~/code/single_cell/colors.r')
source('~/code/single_cell/map_gene.r')
source('~/code/single_cell/tpm.r')


load_signature = function(file=NULL){

    if(!file.exists(file)){file = paste0('~/aviv/db/markers/', file, '.txt')}
    sig = read.table(file, stringsAsFactors=F, row.names=1)
    sig = structure(strsplit(sig[,1], ','), names=rownames(sig))
    return(sig)
}

get_scores = function(seur=NULL, data=NULL, meta=NULL, scores=NULL, names=NULL, regex=NULL, file=NULL, file.cols=NULL, top=NULL, combine='mean'){
    
    # Map variables
    genes = feats = c()
    names = as.list(names)
    regex = as.list(regex)
    if(is.null(data) & !is.null(seur)){data = seur@data}
    if(is.null(meta) & !is.null(seur)){meta = seur@data.info}
    
    # Get names
    if(length(names) > 0){
        genes = names[sapply(names, function(a) all(a %in% rownames(data)))]
	feats = names[sapply(names, function(a) all(a %in% colnames(meta)))]
    }

    # Map regex
    if(length(regex) > 0){
        genes = c(genes, structure(sapply(regex, function(a) grep(a, rownames(data), value=T), simplify=F), names=regex))
	feats = c(feats, structure(sapply(regex, function(a) grep(a, colnames(meta), value=T), simplify=F), names=regex))
    }
    
    # Convert to named lists
    if(length(genes) > 0){
        ni = sapply(genes, function(a) paste(a, collapse='.'))
        if(is.null(names(genes))){
	    genes = structure(genes, names=ni)
	} else {
	    genes = structure(genes, names=ifelse(names(genes) == '', ni, names(genes)))
	}
    }
    
    if(length(feats) > 0){
        ni = sapply(feats, function(a) paste(a, collapse='.'))
	if(is.null(names(feats))){
	    feats = structure(feats, names=ni)
	} else {
	    feats = structure(feats, names=ifelse(names(feats) == '', ni, names(feats)))
	}
    }
    
    # Load signatures from file
    if(!is.null(file)){
        sig = load_signature(file)
	if(!is.null(file.cols)){sig = sig[,file.cols]}
	genes = c(genes, sig)
    }
    
    # Fix lists
    genes = genes[sapply(genes, length) > 0]
    feats = feats[sapply(feats, length) > 0]
    
    # Get gene scores
    if(length(genes) > 0){
    
        # Organism name
	org = predict_organism(rownames(data)[1:100])
	genes = lapply(genes, map_gene, target=org)
	
	# Select top genes
	if(!is.null(top)){genes = lapply(genes, function(a){as.character(na.omit(a[1:top]))})}
	
	# Calculate TPM
	data = calc_tpm(data=data, genes.use=unique(na.omit(as.character(unlist(genes)))))
	
	# Score modules
	si = sapply(genes, function(a){
            i = intersect(rownames(data), a)
            colMeans(data[i,,drop=F], na.rm=T)
        })
	
    	# Log transform
	if(length(scores) > 0){
    	    scores = cbind.data.frame(scores, log2(si + 1))
	} else {
	    scores = log2(si + 1)
	}
    }
    if(length(feats) > 0){
        if(length(scores) > 0){
            scores = cbind.data.frame(scores, as.matrix(as.data.frame(lapply(feats, function(a) meta[,a,drop=F]))))
	} else {
	    scores = as.matrix(as.data.frame(lapply(feats, function(a) meta[,a,drop=F])))
	}
    }
    return(scores)
}


plot_feats = function(seur=NULL, names=NULL, scores=NULL, coords=NULL, data=NULL, meta=NULL, regex=NULL, file=NULL, file.cols=NULL, top=NULL, combine='mean', ...){
    scores = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, top=top, file.cols=file.cols, combine=combine)
    plot_tsne(seur=seur, coords=coords, scores, ...)
}


plot_pcs = function(seur, pcs=NULL, ...){
    if(is.null(pcs)){pcs = 1:seur@data.info$num_pcs[[1]]}
    pcs = seur@pca.rot[,pcs]
    plot_tsne(seur, pcs, ...)
}


plot_clusters = function(seur, ...){
    clusters = sort(grep('Cluster', colnames(seur@data.info), value=T))
    clusters = sapply(seur@data.info[,clusters], function(a){droplevels(as.factor(a))})
    plot_tsne(seur, clusters, ...)
}


plot_tsne = function(seur=NULL, scores=NULL, coords=NULL, ident=TRUE, cells.use=NULL, ymin=0, ymax=1, num_col='auto', pt.size=.75, font.size=11, label.size=5,
	             do.label=T, do.title=TRUE, do.legend=TRUE, na.value='transparent', legend.title='log2(TPM)', out=NULL, nrow=1.5, ncol=1.5, ...){
    
    # Plot columns of scores on TSNE
    # Example: plot_tsne(seur, seur@data.info$nGene)
    # Example: plot_tsne(NULL, seur@data.info$nGene, coords=coords)
    
    # Get xy coordinates
    if(is.null(coords)){
        d = structure(seur@tsne.rot[,1:2], names=c('x', 'y'))
    } else {
        d = structure(as.data.frame(coords[,1:2]), names=c('x', 'y'))
    }
    ps = list()
    
    # Cell identities
    if(!is.logical(ident)){
        d$Identity = ident
    } else if(ident & !is.null(seur)){
        d$Identity = seur@ident
    }
    
    # Cell scores
    if(!is.null(scores)){d = cbind.data.frame(d, scores)}
    
    # Subset cells
    if(!is.null(cells.use)){d = d[cells.use,]}
    d = data.frame(d)
    
    for(col in setdiff(colnames(d), c('x', 'y'))){
	
	if(is.numeric(d[,col])){
	    
	    # Continuous plot
	    u = quantile(d[,col], ymin, na.rm=T)
	    v = quantile(d[,col], ymax, na.rm=T)
	    d[,col][d[,col] < u] = u
	    d[,col][d[,col] > v] = v
	    p = ggplot(d) +
	        geom_point(aes_string(x='x',y='y',colour=col), size=pt.size, ...) +
		scale_colour_gradientn(colours=material.heat(50), guide=guide_colourbar(barwidth=.5, title=legend.title), na.value=na.value) + 
		theme_cowplot(font_size=font.size) +
		xlab('TSNE 1') + ylab('TSNE 2')

	} else {
	    
	    # Discrete plot
	    
	    p = ggplot(d) +
	        geom_point(aes_string(x='x',y='y',colour=col), size=pt.size, ...) +
		theme_cowplot(font_size=font.size) +
		xlab('TSNE 1') + ylab('TSNE 2') +
		scale_colour_manual(values=set.colors, na.value=na.value) + 
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
	if(do.legend == FALSE){
	    p = p + theme(legend.position='none')
	}
	ps[[col]] = p
    }
    
    if(num_col == 'auto'){
        num_col = ceiling(sqrt(length(ps)))
    }
    p = plot_grid(plotlist=ps, ncol=num_col)
    if(is.null(out)){
        p
    } else {
        save_plot(file=out, p, nrow=nrow, ncol=ncol)
    }
}


plot_bars = function(data, pvals=NULL, xlab='Cell_Type', ylab='Frequency', group.by='Samples', palette='Set1'){
    require(tidyverse)
    data = gather_(as.data.frame(data) %>% rownames_to_column(group.by), xlab, ylab, setdiff(colnames(data), group.by))
    p = ggplot(data, aes_string(x=xlab, y=ylab)) +
        geom_bar(aes_string(fill=group.by), stat='identity', position='dodge') +
	scale_fill_brewer(palette=palette)
    
}


plot_heatmap = function(seur, data.use='raw', genes.use=NULL, cells.use=NULL, group_by=NULL, scale='row', ...){

    require(NMF)
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
