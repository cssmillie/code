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


get_scores = function(seur=NULL, data=NULL, meta=NULL, scores=NULL, names=NULL, regex=NULL, file=NULL, file.cols=NULL, top=NULL, cells.use=NULL, group_by=NULL, type='mean', do.log2=TRUE){
    
    # Score gene expression for every cell (type = alpha, mu, or mean)
    
    # Map variables
    genes = feats = gene_scores = feat_scores = c()
    names = as.list(names)
    regex = as.list(regex)
    if(is.null(data) & !is.null(seur)){data = seur@data}
    if(is.null(meta) & !is.null(seur)){meta = seur@data.info}
    
    # Fix input data
    if(!is.null(scores)){scores = as.data.frame(scores)}
    if(!is.null(group_by)){names(group_by) = colnames(seur@data)}
    
    # Subset cells
    if(!is.null(cells.use)){
        data = data[,cells.use]
	meta = meta[cells.use,]
	if(!is.null(group_by)){group_by = group_by[cells.use]}
	if(!is.null(scores)){scores = scores[cells.use,,drop=F]}
    }
    
    # Get genes and feats
    if(length(names) > 0){
        org = predict_organism(rownames(data)[1:100])
        genes = sapply(names, function(a) unique(map_gene(a, target=org)), simplify=F)
	genes = sapply(genes, function(a) a[a %in% rownames(data)])
	feats = sapply(names, function(a) a[a %in% colnames(meta)])
    }
    
    # Add regex terms
    if(length(regex) > 0){
        genes = c(genes, structure(sapply(regex, function(a) grep(a, rownames(data), value=T), simplify=F), names=regex))
	feats = c(feats, structure(sapply(regex, function(a) grep(a, colnames(meta), value=T), simplify=F), names=regex))
    }
    
    # Convert to named lists
    if(length(genes) > 0){
        if(is.null(names(genes))){names(genes) = rep('', length(genes))}
	names(genes) = ifelse(names(genes) == '', sapply(genes, function(a) paste(a, collapse='.')), names(genes))
    }    
    if(length(feats) > 0){
        if(is.null(names(feats))){names(feats) = rep('', length(feats))}
	names(feats) = ifelse(names(feats) == '', sapply(feats, function(a) paste(a, collapse='.')), names(feats))
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
	
	# Score cells
	si = sapply(genes, function(a){
            i = intersect(rownames(data), a)
            colMeans(as.matrix(data[i,,drop=F]), na.rm=T)
        })

	# Convert to alpha, mu, or mean
	if(type == 'alpha'){si = si > 0}
	if(type == 'mu'){si[si == 0] = NA}
	if(do.log2 == TRUE & type %in% c('mu', 'mean')){si = log2(si + 1)}
	
	scores = cbind(scores, si)
    }
    
    if(length(feats) > 0){
        si = as.data.frame(lapply(feats, function(a) meta[,a,drop=F]))
	if(is.null(scores)){
	    scores = si
	} else {
	    scores = cbind.data.frame(scores, si)
	}
    }
    
    # Group scores
    if(!is.null(group_by)){
        scores = data.frame(aggregate(scores, list(as.character(group_by)), mean, na.rm=T), row.names=1)
	if(nrow(scores) == nlevels(group_by)){
	    scores = scores[levels(group_by),,drop=F]
	}
    }
    
    return(scores)
}


plot_tsne = function(seur=NULL, names=NULL, scores=NULL, coords=NULL, data=NULL, meta=NULL, regex=NULL, file=NULL, file.cols=NULL, top=NULL, type='mean', ident=TRUE,
                     cells.use=NULL, ymin=0, ymax=1, num_col='auto', pt.size=.75, font.size=11, do.label=T, label.size=5, do.title=TRUE, title.use=NULL,
	             do.legend=TRUE, legend.title='log2(TPM)', share_legend=FALSE, legend_width=.05, vmin=NA, vmax=NA, na.value='transparent',
		     out=NULL, nrow=1.5, ncol=1.5, ...){
    
    
    # TSNE coordinates
    if(is.null(coords)){
        d = structure(seur@tsne.rot[,1:2], names=c('x', 'y'))
    } else {
        d = structure(as.data.frame(coords[,1:2]), names=c('x', 'y'))
    }
    
    # Cell identities
    if(!is.logical(ident)){
        d$Identity = ident
    } else if(ident & !is.null(seur)){
        d$Identity = seur@ident
    }

    # Subset cells
    if(is.null(cells.use)){cells.use = rownames(d)}
    d = data.frame(d[cells.use,])
    
    # Cell scores
    scores = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, top=top, file.cols=file.cols, type=type, cells.use=cells.use)
    if(!is.null(scores)){d = cbind.data.frame(d, scores)}
    
    # Initialize plotlist
    cat('\nPlotting:', paste(colnames(subset(d, select=-c(x,y))), collapse=', '), '\n')
    ps = list()
    
    # Get limits for shared legend
    if(share_legend == TRUE){
        j = names(which(sapply(subset(d, select=-c(x,y)), is.numeric)))
	cat('\nShared limits:', paste(j, collapse=', '), '\n')
        if(is.na(vmin)){vmin = na.omit(min(d[,j]))}
	if(is.na(vmax)){vmax = na.omit(max(d[,j]))}
	cat('> vmin =', vmin, '\n> vmax =', vmax, '\n')
    }
    
    for(col in setdiff(colnames(d), c('x', 'y'))){
    	
	if(is.numeric(d[,col])){
		    
	    # Continuous plot
	    d[,col] = qtrim(d[,col], qmin=ymin, qmax=ymax)
	    p = ggplot(d) +
	        geom_point(aes_string(x='x',y='y',colour=col), size=pt.size, ...) +
		scale_colour_gradientn(colours=material.heat(50), guide=guide_colourbar(barwidth=.5, title=legend.title), na.value=na.value, limits=c(vmin, vmax)) + 
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
	    if(is.null(title.use)){title = col} else {title = title.use}
	    p = p + ggtitle(title)
	}
	if(do.legend == FALSE){
	    p = p + theme(legend.position='none')
	}
	ps[[col]] = p
    }
    
    if(num_col == 'auto'){
        num_col = ceiling(sqrt(length(ps)))
    }
    
    ps = make_compact(plotlist=ps, num_col=num_col)
    if(length(ps) > 1){
        if(share_legend == TRUE){
            p = share_legend(ps, num_col=num_col, width=legend_width)
        } else {
            p = plot_grid(plotlist=ps, ncol=num_col, align='h')
        }
    }
    
    if(is.null(out)){
        p
    } else {
        save_plot(file=out, p, nrow=nrow, ncol=ncol)
    }
}


make_compact = function(plotlist, num_col, labels=TRUE, ticks=TRUE){
    
    # Make a "plot_grid" plotlist compact by removing axes from interior plots

    # x-axis
    if(length(plotlist) > num_col){
        i = setdiff(1:length(plotlist), rev(1:length(plotlist))[1:min(num_col, length(plotlist))])
	if(labels == TRUE){plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.title.x=element_blank()))}
	if(ticks == TRUE){plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.text.x=element_blank()))}
    }
    
    # y-axis
    if(num_col > 1){
        i = setdiff(1:length(plotlist), which(1:length(plotlist) %% num_col == 1))
	if(labels == TRUE){plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.title.y=element_blank()))}
	if(ticks == TRUE){plotlist[i] = lapply(plotlist[i], function(p) p + theme(axis.text.y=element_blank()))}
    }
    
    return(plotlist)
}


share_legend = function(plotlist, num_col, width=0.05){
    
    # Get first legend in plotlist
    i = min(which(sapply(plotlist, function(p) 'guide-box' %in% ggplotGrob(p)$layout$name)))
    cat(paste('\nUsing shared legend:', names(plotlist)[i], '\n'))
    legend = get_legend(plotlist[[i]])
    
    # Remove all legends
    plotlist = lapply(plotlist, function(p) p + theme(legend.position='none'))
    
    # Make combined plot
    p = plot_grid(plotlist=plotlist, ncol=num_col, align='h')
    p = plot_grid(p, legend, ncol=2, rel_widths=c(1-width, width))
    
    return(p)
}


plot_heatmap = function(seur=NULL, names=NULL, scores=NULL, data=NULL, meta=NULL, regex=NULL, file=NULL, file.cols=NULL, top=NULL, type='mean', group_by=NULL,
                        cells.use=NULL, do.scale=FALSE, scale_method='max', border='black', Rowv=TRUE, Colv=TRUE, labRow=NULL, labCol=NULL, qmin=0, qmax=1, vmin=-Inf, vmax=Inf, out=NA, ...){
    
    require(NMF)
    
    # Cell scores
    if(is.null(group_by)){group_by = seur@ident}
    data = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, top=top, file.cols=file.cols, type=type, cells.use=cells.use, group_by=group_by)
    
    # Scale data
    if(do.scale == TRUE){
        if(scale_method == 'max'){
	    data = as.data.frame(t(t(data)/apply(data, 2, max)))
	} else {
	    data = scale(data)
	}
    }
    
    # Adjust scale
    data = qtrim(data, qmin=qmin, qmax=qmax, vmin=vmin, vmax=vmax)
    
    # Re-order rows
    i = order(apply(data, 1, which.max))
    data = data[i,,drop=F]
    
    # Plot heatmap
    if(ncol(data) == 1){Colv=NA}
    aheatmap(data, scale='none', border_color=border, Rowv=Rowv, Colv=Colv, labRow=labRow, labCol=labCol, hclustfun='average', filename=out, gp=gpar(lineheight=.8), ...)
}


heatmap2 = function(x, col='nmf', lmat=NULL, lwid=NULL, lhei=NULL, margins=c(5,5), rowSize=1, colSize=1, show.tree=TRUE, ...){
    library(gplots)

    # Make gplots heatmap.2 look like aheatmap
    
    # Adjust layout
    if(is.null(lmat)){lmat=matrix(c(0,2,3,1,4,0), nrow=2)}
    if(is.null(lwid)){lwid=c(.05,.9,.05)}
    if(is.null(lhei)){lhei=c(.05,.95)}

    # Label sizes
    cexRow = rowSize*(0.2 + 1/log10(nrow(x)))
    cexCol = colSize*(0.2 + 1/log10(ncol(x)))
    
    # NMF colors
    if(col == 'nmf'){col = colorRampPalette(nmf.colors)(100)} else {col = colorRampPalette(brewer.pal(9, col))(100)}

    # Hide dendrogram
    if(show.tree == FALSE){dendrogram='none'} else {dendrogram='both'}
    
    # Plot data
    heatmap.2(x, trace='none', lmat=lmat, lhei=lhei, lwid=lwid, col=col, cexRow=cexRow, cexCol=cexCol, dendrogram=dendrogram, ...)
}


plot_dots = function(seur=NULL, names=NULL, scores=NULL, data=NULL, meta=NULL, regex=NULL, file=NULL, file.cols=NULL, top=NULL, group_by=NULL, cells.use=NULL, dot_size='alpha', dot_color='mu',
	             rescale=FALSE, reorder=NULL, do.title=TRUE, do.legend=TRUE, xlab='Cell Type', ylab='Gene', out=NULL, nrow=1.5, ncol=1.5, coord_flip=FALSE, max_size=5){

    require(tibble)
    
    # Fix input arguments
    if(is.null(group_by)){group_by = seur@ident}
    if(is.null(cells.use)){cells.use = colnames(seur@data)}
    
    # Scores
    x = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, file.cols=file.cols, top=top, cells.use=cells.use, group_by=group_by, type=dot_size)
    y = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, file.cols=file.cols, top=top, cells.use=cells.use, group_by=group_by, type=dot_color)
    
    # Rescale
    if(rescale == TRUE){
        x.max = apply(x, 2, max, na.rm=T)
	y.max = apply(y, 2, max, na.rm=T)
	x = as.data.frame(t(t(x)/x.max))
	y = as.data.frame(t(t(y)/y.max))
    }
    
    # Gather
    x = x %>% rownames_to_column('Group') %>% gather(Feature, Size, -Group)
    y = y %>% rownames_to_column('Group') %>% gather(Feature, Color, -Group)
    d = cbind(x, y)
    
    # Legend titles
    uppercase = function(x){paste0(toupper(substr(x,1,1)), substr(x,2,nchar(x)))}
    
    # Re-order data
    if(!is.null(reorder)){if(coord_flip == TRUE){d$Group = factor(d$Group, levels=rev(reorder))} else {d$Group = factor(d$Group, levels=reorder)}}
    
    # Dot plot
    p = ggplot(d, aes(x=Group, y=Feature, size=Size, colour=Color)) +
        geom_point() +
	scale_size_area(uppercase(dot_size), max_size=max_size) +
	scale_colour_gradientn(uppercase(dot_color), colours=brewer.pal(9,'YlOrRd')) +
	theme(axis.line = element_blank(), axis.title=element_blank(), panel.grid.major=element_line(colour='black')) +
	background_grid(major='y')
    
    if(coord_flip == TRUE){p = p + coord_flip()}
    
    return(p)
}


plot_violin = function(seur=NULL, names=NULL, scores=NULL, data=NULL, meta=NULL, regex=NULL, file=NULL, file.cols=NULL, top=NULL, type='mean', group_by=NULL, color_by=NULL, pt.size=.25,
	               do.facet=FALSE, facet_by=NULL, facet_scales='free_y', ymin=0, ymax=1, do.scale=FALSE,
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
    scores = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, top=top, file.cols=file.cols, type=type, cells.use=cells.use)
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
	if(!is.character(d[,col])){d[,col] = as.numeric(d[,col])}
        
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
	    scale_fill_manual(values=rev(set.colors)) + theme_cowplot() +
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
	    p = p + facet_grid(as.formula(facet_formula), scales=facet_scales)
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


matrix_barplot = function(data, group_by=NULL, pvals=NULL, xlab='', ylab='Frequency', value='mean', error='se', legend.title='Groups', colors='Paired',
                          out=NULL, nrow=1.5, ncol=1.5, coord_flip=FALSE, sig_only=F, do.facet=F){
    
    # Plot barplot of [M x N] matrix
    # x-axis = matrix columns (e.g. cell types)
    # y-axis = matrix values (e.g. frequencies)
    # fill = matrix rows (e.g. samples) or groups (e.g. conditions)
    
    # Arguments:
    # group.by = the group of each row
    # pvals = [G x N] matrix of p-values for each group and column
    # error = sd, se, or none
    
    require(tidyverse)
    require(data.table)
    require(ggsignif)

    # Groups (default = rows)
    if(is.null(group_by)){group_by = rownames(data)}
    group_by = as.factor(group_by)

    # Select significant comparisons
    if(!is.null(pvals)){
        if(sum(! rownames(pvals) %in% group_by) > 0){stop('rownames(pvals) != group_by')}
	if(sum(! colnames(pvals) %in% colnames(data)) > 0){stop('colnames(pvals) != colnames(data)')}
	if(sig_only == TRUE){
	    j = apply(pvals, 2, min) <= .05
	    data = data[,j]
	    pvals = pvals[,j]
	}
    }
    
    # Construct input data
    names = colnames(data)
    data = data.frame(group=group_by, data)
    group_levels = levels(group_by)
    colnames(data)[2:ncol(data)] = names
    data = as.data.table(gather_(data, 'x', 'y', setdiff(colnames(data), 'group')))
    
    # Value function
    if(value == 'mean'){vf = mean} else if(value == 'median'){vf = median} else {stop()}
    
    # Error function
    se = function(x, na.rm=T){sd(x, na.rm=na.rm)/sqrt(length(x))}    
    if(error == 'sd'){ef = sd} else if(error == 'se'){ef = se} else {ef = function(x){0}}

    # Estimate error bars
    data = data[,.(u=vf(y, na.rm=T), s=ef(y, na.rm=T)),.(group, x)]
    data$x = factor(data$x, levels=names)
    
    # Add p-values
    if(!is.null(pvals)){
        pvals = as.data.frame(pvals) %>% rownames_to_column('group') %>% gather(x, pval, -group) %>% as.data.table()
	setkeyv(data, c('x', 'group'))
	setkeyv(pvals, c('x', 'group'))
	data = merge(data, pvals, all=T)
	data$label = ifelse(data$pval <= .001, '***', ifelse(data$pval <= .01, '**', ifelse(data$pval <= .05, '*', '')))
    }
    data$group = factor(data$group, levels=group_levels)

    # Get colors
    if(length(colors) == 1){colors = brewer.pal(9, colors)}
    
    # Plot data
    p = ggplot(data) + 
    geom_bar(aes(x=x, y=u, fill=group), stat='identity', position=position_dodge(.9)) +
    geom_errorbar(aes(x=x, ymin=u-s, ymax=u+s, fill=group), stat='identity', position=position_dodge(.9), width=.25) +
    scale_fill_manual(values=colors, name=legend.title) + xlab(xlab) + ylab(ylab)
    
    # Facet wrap
    if(do.facet == TRUE){
        p = p + facet_grid(group ~ ., scales='free')
    }

    dy = max(data$u + data$s)*.01
    if(coord_flip == FALSE){
        p = p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
	if(!is.null(pvals)){p = p + geom_text(aes(x=x, y=u+s+dy, label=label, group=group), hjust='center', vjust=0, size=4, angle=0, position=position_dodge(.9))}
    } else {
        p = p + coord_flip()
	if(!is.null(pvals)){p = p + geom_text(aes(x=x, y=u+s+dy, label=label, group=group), hjust='center', vjust=1, size=4, angle=90, position=position_dodge(.9))}
    }
        
    # Save plot
    if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}
    p
}


plot_volcano = function(fcs, pvals, color_by=NULL, facet_by=NULL, labels=NULL, lab.x=c(-2.5, 2.5), lab.nox=c(-1,1), lab.pval=0, lab.n=0, font.size=11, legend.title='Groups', out=NULL, nrow=1.5, ncol=1.5,
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

	# Remove labels
	data$labels[lab.nox[[1]] < data$x & data$x < lab.nox[[2]]] = FALSE
		
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
