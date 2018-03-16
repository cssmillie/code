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


load_signature = function(file=NULL, file.regex=NULL, file.cols=NULL){
    if(!file.exists(file)){file = paste0('~/aviv/db/markers/', file, '.txt')}
    sig = read.table(file, stringsAsFactors=F, row.names=1)
    sig = structure(strsplit(sig[,1], ','), names=rownames(sig))
    if(!is.null(file.regex)){file.cols = grep(file.regex, names(sig), value=T)}
    if(!is.null(file.cols)){sig = sig[file.cols]}
    return(sig)
}


get_scores = function(seur=NULL, data=NULL, meta=NULL, scores=NULL, names=NULL, regex=NULL, file=NULL, file.cols=NULL, file.regex=NULL, top=NULL, cells.use=NULL, group_by=NULL, order_by=NULL,
	              type='mean', do.log=TRUE, log_type='pre', log_zero=1){
    
    # Score gene expression for every cell (type = alpha, mu, or mean)
    
    # Map variables
    genes = feats = c()
    names = as.list(names)
    regex = as.list(regex)
    if(is.null(data) & !is.null(seur)){data = seur@data}
    if(is.null(meta) & !is.null(seur)){meta = seur@data.info}
    
    # Fix input data
    if(!is.null(scores)){scores = as.data.frame(scores, row.names=colnames(seur@data))}
    if(!is.null(group_by)){names(group_by) = colnames(seur@data)}
    name_map = structure(make.names(colnames(scores)), names=colnames(scores))
    
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
	genes = sapply(genes, function(a) a[a %in% rownames(data)], simplify=F)
	feats = sapply(names, function(a) a[a %in% colnames(meta)], simplify=F)
    }
    
    # Add regex terms
    if(length(regex) > 0){
        genes = c(genes, structure(sapply(regex, function(a) grep(a, rownames(data), value=T, perl=T), simplify=F), names=regex))
	feats = c(feats, structure(sapply(regex, function(a) grep(a, colnames(meta), value=T, perl=T), simplify=F), names=regex))
    }
    
    # Convert to named lists
    if(length(genes) > 0){
        if(is.null(names(genes))){names(genes) = rep('', length(genes))}
	names(genes) = ifelse(names(genes) == '', sapply(genes, function(a) paste(a, collapse='.')), names(genes))
	name_map[names(genes)[names(genes) != '']] = make.names(names(genes)[names(genes) != ''])
    }
    
    if(length(feats) > 0){
        if(is.null(names(feats))){names(feats) = rep('', length(feats))}
	names(feats) = ifelse(names(feats) == '', sapply(feats, function(a) paste(a, collapse='.')), names(feats))
	name_map[names(feats)[names(feats) != '']] = make.names(names(feats)[names(feats) != ''])
    }
    
    # Load signatures from file
    if(!is.null(file)){
        sig = load_signature(file, file.regex=file.regex, file.cols=file.cols)
	genes = c(genes, sig)
	name_map[names(sig)] = make.names(names(sig))
    }
    
    # Fix lists
    genes = genes[sapply(genes, length) > 0]
    feats = feats[sapply(feats, length) > 0]
    if(length(genes) > 0){names(genes) = name_map[names(genes)]}
    if(length(feats) > 0){names(feats) = name_map[names(feats)]}
    
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
            colSums(as.matrix(data[i,,drop=F]), na.rm=T)
        })
	
	# Convert to alpha, mu, or mean
	if(type == 'alpha'){si = si > 0}
	if(type == 'mu'){si[si == 0] = NA}
	
	# Combine scores
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
    
    # Calculate log2 before mean
    if((do.log == TRUE & log_type == 'pre') & length(genes) > 0){
        gene_names = names(genes)
	scores[,gene_names] = log2(scores[,gene_names] + log_zero)
    }
    
    # Group scores
    if(!is.null(group_by)){
        scores = data.frame(aggregate(scores, list(as.character(group_by)), mean, na.rm=T), row.names=1)
	if(nrow(scores) == nlevels(group_by)){
	    scores = scores[levels(group_by),,drop=F]
	}
    }
    
    # Calculate log2 after mean
    if((do.log == TRUE & log_type == 'post') & length(genes) > 0){
        gene_names = make.names(names(genes))
	scores[,gene_names] = log2(scores[,gene_names] + log_zero)
    }
    
    # Reverse map
    for(a in names(name_map)){if(a != ''){b = name_map[[a]]; name_map[[b]] = a}}
    if(!is.null(scores)){colnames(scores) = name_map[colnames(scores)]}
    
    return(scores)
}


plot_tsne = function(seur=NULL, names=NULL, scores=NULL, coords=NULL, data=NULL, meta=NULL, regex=NULL, file=NULL, file.cols=NULL, file.regex=NULL, top=NULL, type='mean', ident=TRUE,
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
    scores = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, top=top, file.cols=file.cols, file.regex=file.regex, type=type, cells.use=cells.use)
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

	    # Get colors
	    if(do.label == TRUE){tsne.colors = tsne.colors} else {tsne.colors = set.colors}
	    
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


plot_heatmap = function(seur=NULL, names=NULL, scores=NULL, data=NULL, meta=NULL, regex=NULL, file=NULL, file.regex=NULL, file.cols=NULL, top=NULL, type='mean', group_by=NULL,
                        cells.use=NULL, do.scale=FALSE, scale_method='max', border='black', Rowv=TRUE, Colv=TRUE, labRow=NULL, labCol=NULL, qmin=0, qmax=1, vmin=-Inf, vmax=Inf, out=NA, ...){
    
    require(NMF)
    
    # Cell scores
    if(is.null(group_by)){group_by = seur@ident}
    data = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, top=top, file.regex=file.regex, file.cols=file.cols, type=type, cells.use=cells.use, group_by=group_by)
    
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


ggheatmap = function(data, Rowv='hclust', Colv='hclust', xlab='', ylab='', xsec=FALSE, ysec=FALSE, xstag=FALSE, ystag=FALSE, title='', legend.title='', pal='nmf', do.legend=TRUE, font_size=14,
                     out=NULL, nrow=1.25, ncol=1.25, qmin=0, qmax=1, vmin=-Inf, vmax=Inf, symm=FALSE, xstrip=NULL, ystrip=NULL, hclust_met='average', outline='#cccccc', replace_na=NA,
		     labRow=TRUE, labCol=TRUE){
    require(ggplot2)
    require(tidyr)
    require(cowplot)
    require(tibble)
    
    # Scale values
    data = qtrim(data, qmin=qmin, qmax=qmax, vmin=vmin, vmax=vmax)
    data[is.na(data)] = replace_na
    if(!is.logical(labRow)){rownames(data) = labRow}
    if(!is.logical(labCol)){colnames(data) = labCol}
    
    # Convert to long format
    data = as.data.frame(data)
    x = data %>% rownames_to_column('row') %>% gather(col, value, -row)
    x$value = as.numeric(x$value)

    # Order rows
    if(length(Rowv) > 1){rowv = Rowv; Rowv = 'none'} else {rowv = rev(rownames(data))}
    if(length(Colv) > 1){colv = Colv; Colv = 'none'} else {colv = colnames(data)}
    if(nrow(data) <= 2){Rowv = 'none'}
    if(ncol(data) <= 2){Colv = 'none'}
    if(Rowv == 'hclust'){
        rowv = rev(rownames(data)[hclust(dist(data), method=hclust_met)$order])
    }
    if(Colv == 'hclust'){
        colv = colnames(data)[hclust(dist(t(data)), method=hclust_met)$order]
    }
    if(Rowv == 'none'){
        rowv = rev(rownames(data))
    }
    if(Colv == 'none'){
        colv = colnames(data)
    }
    if(Rowv == 'min'){
        rowv = rev(rownames(data)[order(apply(data, 1, which.min))])
    }
    if(Rowv == 'max'){
        rowv = rev(rownames(data)[order(apply(data, 1, which.max))])
    }
    if(Colv == 'min'){
        i = match(rownames(data)[apply(data, 2, which.min)], rowv)
	colv = rev(colnames(data)[order(i)])        
    }
    if(Colv == 'max'){        
        i = match(rownames(data)[apply(data, 2, which.max)], rowv)
	colv = rev(colnames(data)[order(i)])
    }
    Rowv = rowv
    Colv = colv

    # Set order of row and column labels
    x$row = factor(x$row, levels=Rowv)
    x$col = factor(x$col, levels=Colv)
    
    # Get odd/even indices
    r1 = seq(1, length(Rowv), by=2)
    r2 = seq(2, length(Rowv), by=2)
    c1 = seq(1, length(Colv), by=2)
    c2 = seq(2, length(Colv), by=2)
    
    # Stagger row and column labels
    if(ystag == TRUE){ystag = 4*max(nchar(as.character(x$row[r1])))} else {ystag = 0}
    if(xstag == TRUE){xstag = 4*max(nchar(as.character(x$col[c1])))} else {xstag = 0}
    levels(x$row)[r2] = paste0(levels(x$row)[r2], strrep(' ', ystag))
    Rowv = levels(x$row)
    levels(x$col)[c2] = paste0(levels(x$col)[c2], strrep(' ', xstag))
    Colv = levels(x$col)
    
    # Get plot data
    if(length(pal)==1){if(pal == 'nmf'){pal = rev(colorRampPalette(nmf.colors)(100))[10:100]} else {pal = colorRampPalette(brewer.pal(9, pal))(100)}}
    
    # Plot with geom_tile
    p = ggplot(x) +
        geom_tile(aes(x=as.numeric(col), y=as.numeric(row), fill=value), colour=outline) +
	labs(x=xlab, y=ylab, title=title, fill=legend.title) +
	theme_cowplot(font_size=font_size) +
	theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), axis.line=element_blank())
    
    # Scale
    if(symm == TRUE){
        max_value = max(abs(x$value))
	values = seq(-max_value, max_value, length=3)
        p = p + scale_fill_gradientn(colours=pal, values=values, rescaler=function(x, ...) x, oob=identity)
    } else {
        p = p + scale_fill_gradientn(colours=pal)
    }
	
    # Secondary x-axis
    if(xsec == FALSE){
        p = p + scale_x_continuous(breaks=1:length(Colv), labels=Colv, expand=c(0,0))
    } else {
	p = p + scale_x_continuous(breaks=c1, labels=Colv[c1], sec.axis=dup_axis(breaks=c2, labels=Colv[c2]), expand=c(0,0)) + theme(axis.text.x.top= element_text(angle=90, hjust=0, vjust=.5))
    }
    
    # Secondary y-axis
    if(ysec == FALSE){
        p = p + scale_y_continuous(breaks=1:length(Rowv), labels=Rowv, expand=c(0,0))
    } else {
	p = p + scale_y_continuous(breaks=r1, labels=Rowv[r1], sec.axis=dup_axis(breaks=r2, labels=Rowv[r2]), expand=c(0,0))
    }
    
    if(labRow[[1]] == FALSE){p = p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())}
    if(labCol[[1]] == FALSE){p = p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())}
    if(do.legend == FALSE){p = p + theme(legend.position='none')}

    # Save plot
    if(!is.null(out)){
        save_plot(p, file=out, nrow=nrow, ncol=ncol)
    }
    p
}


plot_dots = function(seur=NULL, names=NULL, scores=NULL, data=NULL, meta=NULL, regex=NULL, file=NULL, file.regex=NULL, file.cols=NULL, top=NULL, group_by=NULL, cells.use=NULL, dot_size='alpha', dot_color='mu',
	             rescale=FALSE, reorder=NULL, do.title=TRUE, do.legend=TRUE, xlab='Cell Type', ylab='Gene', out=NULL, nrow=1.5, ncol=1.5, coord_flip=FALSE, max_size=5){

    require(tibble)
    
    # Fix input arguments
    if(is.null(group_by)){group_by = seur@ident}
    if(is.null(cells.use)){cells.use = colnames(seur@data)}
    
    # Scores
    x = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, file.regex=file.regex, file.cols=file.cols, top=top, cells.use=cells.use, group_by=group_by, type=dot_size)
    y = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, file.regex=file.regex, file.cols=file.cols, top=top, cells.use=cells.use, group_by=group_by, type=dot_color)
    
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


plot_violin = function(seur=NULL, names=NULL, scores=NULL, data=NULL, meta=NULL, regex=NULL, file=NULL, file.regex=NULL, file.cols=NULL, top=NULL, type='mean', group_by=NULL, color_by=NULL, pt.size=.25,
	               do.facet=FALSE, facet_by=NULL, facet_formula=NULL, facet_scales='free_y', ymin=0, ymax=1, do.scale=FALSE, num_col='auto',
                       ident=TRUE, cells.use=NULL, do.title=TRUE, do.legend=FALSE, xlab='Cell Type', ylab='log2(TPM)', out=NULL, nrow=1.5, ncol=1.5, legend.title='Group',
		       coord_flip=FALSE, alpha=.5, order=NULL){
    
    # Initialize plots
    ps = list()

    # Get facet formula
    if(is.null(facet_formula)){
    if(do.facet == FALSE){
        facet_formula = ifelse(is.null(facet_by), '~ .', 'Facet ~ .')
    } else {
        do.title = FALSE
        facet_formula = ifelse(is.null(facet_by), 'Feature ~ .', 'Feature ~ Facet + .')
    }
    } else {facet_formula = as.formula(facet_formula)}
    
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
    scores = get_scores(seur=seur, data=data, meta=meta, scores=scores, names=names, regex=regex, file=file, top=top, file.regex=file.regex, file.cols=file.cols, type=type, cells.use=cells.use)
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
	    p = p + facet_grid(as.formula(facet_formula), scales=facet_scales)
	}
	
	ps[[col]] = p
    }

    if(num_col == 'auto'){
        num_col = ceiling(sqrt(length(ps)))
    }
    ps = make_compact(ps, num_col=num_col)
    p = plot_grid(plotlist=ps, ncol=num_col)
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
    if(nlevels(group_by) == 0){group_by = as.factor(group_by)}
    
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


plot_volcano = function(fcs, pvals, labels=NULL, color_by=NULL, facet_by=NULL, lab.x=c(-1, 1), max_pval=.05, lab.n=0, lab.p=0, lab.size=5, do.repel=TRUE, pval_floor=-Inf,
                        font.size=11, legend.title='Groups', out=NULL, nrow=1.5, ncol=1.5, xlab='Fold change', ylab='-log10(pval)', palette='RdBu', ret.labs=FALSE){
    
    # Volcano plot with automatic labeling:
    # best p-values in region lab.x (n = lab.n)
    # best p-values globally (n = lab.p)
    
    # Make plot data
    if(is.null(color_by)){color_by = rep('Group', length(fcs))}
    if(is.null(facet_by)){facet_by = rep('Group', length(fcs))}
    pvals[pvals < pval_floor] = pval_floor
    data = data.frame(x=fcs, y=-log10(pvals), Color=color_by, Facet=facet_by, stringsAsFactors=F)
    data$Facet = as.factor(data$Facet)
    
    if(!is.null(labels)){        
        for(facet in unique(data$Facet)){
		    
	    # Select facet data
	    i = which(data$Facet == facet & 10**(-data$y) <= max_pval)
	    di = data[i,]
	    
	    # Label points < lab.x
	    breaks = seq(from=min(di$x, na.rm=T), to=lab.x[[1]], length.out=10)
	    groups = cut(di$x, breaks=breaks, include.lowest=TRUE)
	    j1 = as.numeric(simple_downsample(cells=1:nrow(di), groups=groups, ngene=di$y, total_cells=as.integer(lab.n/2)))
	    
	    # Label points > lab.x
	    breaks = seq(from=lab.x[[2]], to=max(di$x, na.rm=T), length.out=10)
	    groups = cut(di$x, breaks=breaks, include.lowest=TRUE)
	    j2 = as.numeric(simple_downsample(cells=1:nrow(di), groups=groups, ngene=di$y, total_cells=as.integer(lab.n/2)))
	    
	    # Label best global p-values
	    j3 = which(order(-1*di$y) <= lab.p)

	    # Set labels
	    j = sort(unique(c(j1, j2, j3)))
	    data[i[j], 'labels'] = TRUE
	}
    } else {
        data$labels = FALSE
    }	
    
    # Set labels
    data$labels = ifelse(data$labels, labels, '')
    
    # Make volcano plot
    p = ggplot(data, aes(x=x, y=y)) +
        geom_point(aes(colour=Color)) +
	theme_cowplot(font_size=font.size) +
	xlab(xlab) +
	ylab(ylab)
    
    # Add labels
    if(do.repel == TRUE){
        p = p + geom_text_repel(aes(label=labels), size=lab.size, segment.color='#cccccc')
    } else {
        p = p + geom_text(aes(label=labels), size=lab.size)
    }
    
    # Add color scale
    if(is.numeric(data$Color)){
        p = p + scale_colour_distiller(palette=palette, trans='reverse')
    } else {
        p = p + scale_colour_manual(values=desat(set.colors, .5))
    }
    if(all(color_by == 'Group')){p = p + theme(legend.position='none')}
        
    # Facet wrap
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

    # Save to file
    if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}
    if(ret.labs){as.character(na.omit(data$labels))} else {p}
}


plot_volcanos = function(markers, color_by=NULL, outdir='volcano'){
    
    # Create output directory
    if(!file.exists(outdir)){dir.create(outdir)}
    if(!is.null(color_by) & (nrow(markers) != length(color_by))){stop()}
    
    # Split markers by contrast
    m = split(markers, markers$contrast)
    q = split(color_by, markers$contrast)
    
    # Make each volcano plot
    for(name in names(m)){
        out = gsub(':', '.', paste0(outdir, '/', name, '.volcano.png'))
	print(c(name, out))
        plot_volcano(m[[name]]$coefD, m[[name]]$pvalD, labels=m[[name]]$gene, color_by=q[[name]], lab.x=c(-.5, .5), max_pval=.05, lab.n=100, lab.p=20, lab.size=3, out=out, nrow=2, ncol=2)	
    }
}


simple_scatter = function(x, y, lab=NULL, col=NULL, lab.use=NULL, lab.near=0,  lab.n=0, lab.g=0, groups=NULL, lab.size=4, lab.type='up', palette='RdBu', xlab=NULL, ylab=NULL, out=NULL, nrow=1, ncol=1){

    require(ggplot2)
    require(ggrepel)
    require(cowplot)
    require(wordspace)

    if(is.null(col)){col = rep('', length(x))}
    d = data.frame(x=x, y=y, lab=lab, col=col, stringsAsFactors=FALSE)
    
    if(is.null(xlab)){xlab = deparse(substitute(x))}
    if(is.null(ylab)){ylab = deparse(substitute(y))}
    
    if(lab.n > 0){
        if(is.null(groups)){
            breaks = seq(from=min(d$x, na.rm=T), to=max(d$x, na.rm=T), length.out=20)
	    groups = cut(d$x, breaks=breaks)
	}
	i = c()
	if(lab.type %in% c('both', 'up')){
	    i = c(i, as.numeric(simple_downsample(cells=1:nrow(d), groups=groups, ngene=d$y, total_cells=lab.n)))
	    i = c(i, which(order(-1*d$y) <= lab.g))
	}
	if(lab.type %in% c('both', 'down')){
	    i = c(i, as.numeric(simple_downsample(cells=1:nrow(d), groups=groups, ngene=-1*d$y, total_cells=lab.n)))
	    i = c(i, which(order(d$y) <= lab.g))
	}
	d[unique(i), 'col'] = 'lab.n'
    }
    
    if(!is.null(lab.use)){

        # label neighbors
	u = as.matrix(d[, c('x','y')])
	v = as.matrix(d[lab %in% lab.use, c('x','y')])
	i = unique(sort(unlist(apply(dist.matrix(u,v,skip.missing=TRUE,method='euclidean'), 2, order)[1:(lab.near + 1),])))
	d$col[i] = 'lab.near'
	
        # label points
        d$col[d$lab %in% lab.use] = 'lab.use'
    }
    d$lab[d$col == ''] = ''
    d = d[order(d$col),]
    levels(d$col) = unique(c('', 'lab.n', 'lab.use', 'lab.near', sort(unique(col))))
    
    p = ggplot(d, aes(x=x, y=y)) +
        geom_point(aes(colour=col)) +
	geom_text_repel(aes(label=lab), size=lab.size, segment.color='grey') +
	theme_cowplot() + xlab(xlab) + ylab(ylab)
    
    if(all(col == '')){
        p = p + scale_colour_manual(values=c('lightgray', 'black', 'red', 'pink')) + theme(legend.position='none')	
    } else if(is.numeric(d$col)){
        p = p + scale_colour_distiller(palette=palette, trans='reverse')
    } else {
        p = p + scale_colour_manual()
    }
    
    if(!is.null(out)){save_plot(plot=p, filename=out, nrow=nrow, ncol=ncol)}

    return(p)
}
