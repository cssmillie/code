require(ggplot2)
require(ggrepel)
require(gtools)
require(cowplot)
require(tidyr)

source('~/code/single_cell/colors.r')
source('~/code/single_cell/map_gene.r')
source('~/code/single_cell/tpm.r')


qtrim = function(x, qmin=0, qmax=1, vmin=-Inf, vmax=Inf){
    
    # Trim by value
    x[x < vmin] = vmin
    x[x > vmax] = vmax
    
    # Trim by quantile
    u = quantile(x, qmin, na.rm=T)
    v = quantile(x, qmax, na.rm=T)
    x[x < u] = u
    x[x > v] = v
    
    return(x)
}


load_signature = function(file=NULL){

    if(!file.exists(file)){file = paste0('~/aviv/db/markers/', file, '.txt')}
    sig = read.table(file, stringsAsFactors=F, row.names=1)
    sig = structure(strsplit(sig[,1], ','), names=rownames(sig))
    return(sig)
}


get_scores = function(seur=NULL, data=NULL, meta=NULL, scores=NULL, names=NULL, regex=NULL, file=NULL, file.cols=NULL, top=NULL, combine='mean', cells.use=NULL, group_by=NULL){
    
    # Map variables
    genes = feats = c()
    names = as.list(names)
    regex = as.list(regex)
    if(is.null(data) & !is.null(seur)){data = seur@data}
    if(is.null(meta) & !is.null(seur)){meta = seur@data.info}
    
    # Subset data
    if(!is.null(cells.use)){
        data = data[,cells.use]
	meta = meta[cells.use,]
	if(!is.null(group_by)){
	    names(group_by) = colnames(seur@data)
	    group_by = group_by[cells.use]
	}
    }
    
    # Get names
    if(length(names) > 0){
        genes = sapply(names, function(a) unique(map_gene(a, target=predict_organism(rownames(data)[1:100]))))
        genes = genes[sapply(genes, function(a) all(a %in% rownames(data)))]
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
	if(!is.null(file.cols)){sig = sig[file.cols]}
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
            colMeans(as.matrix(data[i,,drop=F]), na.rm=T)
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
	    scores = cbind.data.frame(scores, lapply(feats, function(a) meta[,a,drop=F]))
	} else {
	    scores = as.data.frame(lapply(feats, function(a) meta[,a,drop=F]))
	}
    }
    
    # Group scores
    if(!is.null(group_by)){
        scores = t(data.frame(aggregate(scores, list(as.character(group_by)), mean), row.names=1))
    }
    
    return(scores)
}


plot_tsne = function(seur=NULL, names=NULL, scores=NULL, coords=NULL, data=NULL, meta=NULL, regex=NULL, file=NULL, file.cols=NULL, top=NULL, combine='mean',
                     ident=TRUE, cells.use=NULL, ymin=0, ymax=1, num_col='auto', pt.size=.75, font.size=11, label.size=5,
	             do.label=T, do.title=TRUE, do.legend=TRUE, na.value='transparent', legend.title='log2(TPM)', out=NULL, nrow=1.5, ncol=1.5, ...){
    
    # Cell scores
    scores = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, top=top, file.cols=file.cols, combine=combine)
    
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
	    d[,col] = qtrim(d[,col], qmin=ymin, qmax=ymax)
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
		scale_colour_manual(values=tsne.colors, na.value=na.value) + 
		theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	    
	    if(do.label == T){
	        t = aggregate(d[,c('x', 'y')], list(d[,col]), median)
		colnames(t) = c('l', 'x', 'y')
		p = p + geom_text_repel(data=t, aes(x=x, y=y, label=l, lineheight=.8), point.padding=NA, size=label.size, family='Helvetica') + theme(legend.position='none')
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


plot_heatmap = function(seur=NULL, names=NULL, scores=NULL, data=NULL, meta=NULL, regex=NULL, file=NULL, file.cols=NULL, top=NULL, combine='mean', group_by=NULL,
                        cells.use=NULL, do.scale=TRUE, border_color='black', Rowv=TRUE, Colv=TRUE, qmin=0, qmax=1, vmin=-Inf, vmax=Inf, out=NA, ...){
    
    require(NMF)
    
    # Cell scores
    if(is.null(group_by)){group_by = seur@ident}
    data = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, top=top, file.cols=file.cols, combine=combine, cells.use=cells.use, group_by=group_by)
    
    # Scale data
    if(do.scale == TRUE){
        data = t(scale(t(data)))
    }
    
    # Adjust scale
    data = qtrim(data, qmin=qmin, qmax=qmax, vmin=vmin, vmax=vmax)
    
    # Cluster columns
    d = hclust(dist(t(data), method='euclidean'), method='average')
    data = data[,rev(d$order)]
    
    # Re-order rows
    i = order(apply(data, 1, which.max))
    data = data[i,]
    
    # Plot heatmap
    aheatmap(data, scale='none', border_color=border_color, Rowv=Rowv, Colv=Colv, hclustfun='average', filename=out, gp=gpar(lineheight=.8))
}


plot_violin = function(seur=NULL, names=NULL, scores=NULL, data=NULL, meta=NULL, regex=NULL, file=NULL, file.cols=NULL, top=NULL, combine='mean', group_by=NULL, color_by=NULL, pt.size=.25,
	               do.facet=FALSE, facet_by=NULL, ymin=0, ymax=1, do.scale=FALSE,
                       ident=TRUE, cells.use=NULL, do.title=TRUE, do.legend=FALSE, xlab='Cell Type', ylab='log2(TPM)', out=NULL, nrow=1.5, ncol=1.5, legend.title='Group',
		       coord_flip=FALSE, alpha=.5){
    
    # Initialize plots
    ps = list()

    # Get facet formula
    if(do.facet == FALSE){
        facet_formula = ifelse(is.null(facet_by), '~ .', 'Facet ~ .')
    } else {
        do.title = FALSE
        facet_formula = ifelse(is.null(facet_by), 'Feature ~ .', 'Feature ~ Facet + .')
    }
    
    # Fix input arguments
    if(is.null(group_by)){group_by = seur@ident}
    if(is.null(color_by)){color_by = group_by}
    if(is.null(facet_by)){facet_by = rep('', ncol(seur@data))}
    if(is.null(cells.use)){cells.use = colnames(seur@data)}
    if(is.null(out)){alpha=1} else {alpha=alpha}
    
    # Plot data
    d = data.frame(Group=as.factor(group_by), Color=as.factor(color_by), Facet=as.factor(facet_by), row.names=colnames(seur@data))
    d = d[cells.use,,drop=F]
    if(coord_flip == TRUE){d$Group = factor(d$Group, levels=rev(levels(d$Group)))}
    
    # Cell scores
    scores = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, top=top, file.cols=file.cols, combine=combine, cells.use=cells.use)
    scores = scores[cells.use,,drop=F]
    d = cbind.data.frame(d, scores)
    
    # Fix NAs
    d = d[!is.na(d$Facet),,drop=F]
    
    # Scale data
    j = sapply(d, function(a) !is.factor(a))
    d[,j] = apply(d[,j,drop=F], 2, function(a) qtrim(a, qmin=ymin, qmax=ymax))
    
    # Facet data
    if(do.facet == TRUE){
        d = gather(d, Feature, Value, -Group, -Color, -Facet)
    }
    
    # Violin plots
    for(col in setdiff(colnames(d), c('Group', 'Color', 'Facet', 'Value'))){

        # Convert to numeric?
	d[,col] = as.numeric(d[,col])
        
        # Facet if necessary
        if(do.facet == TRUE){
	    p = ggplot(data=d, aes_string(x='Group', y='Value', fill='Color'))
	} else {
	    p = ggplot(data=d, aes_string(x='Group', y=col, fill='Color'))
	}
	
	# Make violin plot
	p = p +
	    geom_point(position=position_jitterdodge(dodge.width=0.6, jitter.width=2.5), size=pt.size, show.legend=F) +
	    geom_violin(scale='width', alpha=alpha) +
	    scale_fill_manual(values=set.colors) + theme_cowplot() +
	    xlab(xlab) + ylab(ylab) + labs(fill=legend.title) +
	    stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, geom='crossbar', width=.5, show.legend=F) +
	    theme(axis.text.x = element_text(angle = 45, hjust = 1))
	
	if(do.title == TRUE){
	    p = p + ggtitle(col)
	}

	if(do.legend == FALSE){
	    p = p + guides(fill=FALSE)
	}

	if(coord_flip == TRUE){
	    p = p + coord_flip()
	}

	if(facet_formula != '~ .'){
	    p = p + facet_grid(as.formula(facet_formula))
	}
	
	ps[[col]] = p
    }
    
    p = plot_grid(plotlist=ps, ncol=ceiling(length(ps)/3))
    if(is.null(out)){
        p
    } else {
        save_plot(file=out, p, nrow=nrow, ncol=ncol)
    }
}




plot_feats = function(...){plot_tsne(...)}

plot_pcs = function(seur, pcs=NULL, ...){
    if(is.null(pcs)){pcs = 1:seur@data.info$num_pcs[[1]]}
    pcs = seur@pca.rot[,pcs]
    plot_tsne(seur, scores=pcs, ...)
}

plot_clusters = function(seur, ...){
    clusters = sort(grep('Cluster', colnames(seur@data.info), value=T))
    clusters = sapply(seur@data.info[,clusters], function(a){droplevels(as.factor(a))})
    plot_tsne(seur, scores=clusters, ...)
}


matrix_barplot = function(data, group_by=NULL, pvals=NULL, xlab='', ylab='Frequency', error='se', legend.title='Groups', palette='Paired', out=NULL, nrow=1.5, ncol=1.5, coord_flip=FALSE){
    
    # Plot barplot of [M x N] matrix
    # x-axis = matrix columns (e.g. cell types)
    # y-axis = matrix values (e.g. frequencies)
    # fill = matrix rows (e.g. samples) or groups (e.g. conditions)
    
    # Arguments:
    # group.by = the group of each row
    # pvals = named list of p-values for each column
    # error = sd, se, or none
    
    require(tidyverse)
    require(data.table)
    require(ggsignif)

    # Groups (default = rows)
    if(is.null(group_by)){group_by = rownames(data)}
    
    # Construct input data
    names = colnames(data)
    data = data.frame(groups=group_by, data)
    colnames(data)[2:ncol(data)] = names
    data = as.data.table(gather_(data, 'x', 'y', setdiff(colnames(data), 'groups')))
    
    # Error function
    se = function(x){sd(x)/sqrt(length(x))}    
    if(error == 'sd'){ef = sd} else if(error == 'se'){ef = se} else {ef = function(x){0}}

    # Estimate error bars
    data = data[,.(u=mean(y), s=ef(y)),.(groups, x)]
    data$x = factor(data$x, levels=names)
    
    # Plot data
    p = ggplot(data) +
    geom_bar(aes(x=x, y=u, fill=groups), stat='identity', position=position_dodge(.9)) +
    geom_errorbar(aes(x=x, ymin=u-s, ymax=u+s, fill=groups), stat='identity', position=position_dodge(.9), width=.25) +
    scale_fill_manual(values=set.colors, name=legend.title) + xlab(xlab) + ylab(ylab)
    
    # P-values
    if(!is.null(pvals)){
        pvals = as.data.frame(pvals) %>% rownames_to_column('x')
        pvals$label = ifelse(pvals$pvals <= .001, '***', ifelse(pvals$pvals <= .01, '**', ifelse(pvals$pvals <= .05, '*', '')))
        pvals$pos = -.025*max(data$u)
        p = p + geom_text(data=pvals, aes(x=x, y=pos, label=label), hjust='center', vjust='center', size=4, angle=0, nudge_x=-.1)
    }
    
    if(coord_flip == TRUE){p = p + coord_flip()}
    
    # Save plot
    if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}
    p
}


plot_volcano = function(fcs, pvals, color_by=NULL, facet_by=NULL, labels=NULL, lab.x=c(-2.5, 2.5), lab.pval=0, lab.n=0, font.size=11, legend.title='Groups', out=NULL, nrow=1.5, ncol=1.5,
                        xlab='Fold change', ylab='-log10(pval)'){

    # Create a volcano plot using supplied p-values and fold-changes
    # Label points with fold-changes defined by lab.x and pval < lab.pval
    # Label points with best p-values along entire fold-change axis

    # Fix input arguments
    if(is.null(color_by)){color_by = rep('Group', length(fcs))}
    if(is.null(facet_by)){facet_by = rep('Group', length(fcs))}
    
    data = data.frame(x=fcs, y=-log10(pvals), Color=color_by, Facet=facet_by, stringsAsFactors=F)
    data$Facet = as.factor(data$Facet)
    
    # Select labels to keep
    if(!is.null(labels)){

        # Hard cutoffs
        data$labels = FALSE
        data$labels[data$x < lab.x[[1]] & data$y >= -log10(lab.pval)] = TRUE
        data$labels[data$x > lab.x[[2]] & data$y >= -log10(lab.pval)] = TRUE
	
	# Select top N per facet
	for(facet in unique(data$Facet)){
	    i = which(data$Facet == facet)
	    u = sort(unname(unlist(tapply(1:nrow(data[i,]), cut(data[i,]$x, lab.n), function(j) j[which.max(data[i,]$y[j])]))))
	    data[i[u], 'labels'] = TRUE
	}
    } else {
        data$labels = FALSE
    }	
    
    # Set labels
    data$labels = ifelse(data$labels, labels, '')
    
    # Make volcano plot
    p = ggplot(data, aes(x=x, y=y)) + geom_point(aes(colour=Color)) + geom_text_repel(aes(label=labels)) + theme_cowplot(font_size=font.size) + xlab(xlab) + ylab(ylab) +
        scale_colour_manual(values=tsne.colors)
    
    if(nlevels(data$Facet) > 1){
        print('Faceting')
        p = p + facet_wrap(~ Facet, scales='free') +
	    theme(axis.line = element_line(colour = "black"),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    panel.background = element_blank(),
	    strip.background = element_blank())
	    
    }
    
    if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}
    p
}
