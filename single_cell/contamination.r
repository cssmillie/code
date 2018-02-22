plot_contamination = function(u1, u2, coefs, residuals, lab.use=NULL, lab.fit=NULL, fit.cutoff=2, out=NULL){
    
    require(ggplot2)
    require(ggrepel)
    require(cowplot)
    
    # log-transform
    l1 = log2(u1 + .5*min(u1[u1 > 0]))
    l2 = log2(u2 + .5*min(u2[u2 > 0]))
    
    # contamination
    lab.con = names(which(residuals < fit.cutoff))
    
    # scatterplot data
    d = data.frame(x=l2, y=l1, lab=ifelse(names(l1) %in% i, names(l1), ''), Type=rep('Other', length(l1)), stringsAsFactors=F)
    d[lab.con, 'Type'] = 'Contamination'
    d[lab.fit, 'Type'] = 'Fit'
    lab.use = intersect(rownames(d), lab.use)
    d[lab.use, 'Type'] = 'Label'
    d[lab.use, 'lab'] = lab.use
    
    # rug and line data
    d.rug = data.frame(x=l2[u1 == 0], y=l2[u1 == 0])
    x0 = (min(l1, na.rm=T) - coefs[[1]] - fit.cutoff)/coefs[[2]]
    x1 = max(l2, na.rm=T)
    d.line = data.frame(x=c(x0, x1))
    d.line$y = coefs[[2]]*d.line$x + coefs[[1]] + fit.cutoff
    
    # plot data
    if(!is.null(out)){alpha=.25} else {alpha=1}
    p = ggplot(d, aes(x=x, y=y)) +
        geom_point(aes(colour=Type)) +
   	geom_text_repel(aes(label=lab), size=2, segment.color='grey') +
	geom_rug(data=d.rug, aes(x=x), sides='t', col='black', alpha=alpha) +
	geom_line(data=d.line, aes(x=x, y=y), lty=2) +
	xlab(paste0('Mean TPM (non-', group, ')')) +
    	ylab(paste0('Mean TPM (', group, ')')) +
    	scale_colour_manual(values=c('lightcoral', 'black', 'steelblue3', 'lightgray')) +	
    	theme_cowplot()
    
    # save or display plot
    if(!is.null(out)){
        save_plot(p, file=out, nrow=2.25, ncol=2.5)
    } else {
        p
    }
}


detect_contamination = function(tpm, groups, samples, global_coefs=NULL, fit.n=50, fit.cutoff=2, do.plot=TRUE, lab.use=NULL, prefix='test'){
    
    require(MASS)
    source('~/code/single_cell/regression.r')
    
    # summarize
    cat('\n\nDetecting ambient contamination\n\n')
    print(table(groups))
    print(table(samples))
    if(!is.null(global_coefs)){coefs = global_coefs}
    
    # initialize variables
    res = list()
    
    # iterate over groups
    for(group in unique(groups)){cat(paste0('\n\nGroup = ', group))
        
        # output file
        out = paste(prefix, group, 'fit.png', sep='.')
    	
        # subset data
	i = groups == group
	j = groups != group

	# sample frequencies
	f = table(as.factor(samples)[i])
	f = as.matrix(f/sum(f))

	# group mean
	u1 = rowMeans(tpm[,i])
	
	# other mean
	u2 = sapply(unique(samples), function(a){
	    rowSums(tpm[,j & (samples == a)])/sum(samples == a)
	})
	
	stopifnot(colnames(u2) == colnames(f))
	u2 = (u2 %*% f)[,1]
	
	# log-transform
	nice_log2 = function(x){y = log2(x); y[is.infinite(y)] = NA; y}
	l1 = nice_log2(u1)
	l2 = nice_log2(u2)
	
	if(is.null(global_coefs)){
	    
	    # fit boundaries
	    lo = quantile(l2[u1 == 0], .99, na.rm=T)
	    hi = sort(l2, decreasing=T)[100]
	    cat(paste0('\n\tLo Cutoff = ', lo, '\n\tHi Cutoff = ', hi))
	    exclude = list(c(-Inf, lo), c(hi, Inf))
	    
	    # select points for regression
	    lab.fit = names(select_points(l2, l1, n=fit.n, dir='down', nbins=10, loess=T, exclude=exclude))
	    cat(paste0('\n\tGenes for regression: ', paste(lab.fit, collapse=', ')))
	    
	    # robust linear model
	    cat('\n\tFitting rlm')
	    fit = rlm(l1[lab.fit] ~ l2[lab.fit])
	    coefs = coef(fit)
	    print(coefs)

	} else {fit = lab.fit = NULL}
	
	# calculate residuals
	residuals = l1 - (coefs[[2]]*l2 + coefs[[1]])
	lab.con = names(which(residuals < fit.cutoff))
	cat(paste0('\n\tLikely contaminants: ', paste(lab.con, collapse=', ')))
	
	# plot data
	if(do.plot == TRUE){
	    plot_contamination(u1, u2, coefs, residuals, lab.use=lab.use, lab.fit=lab.fit, fit.cutoff=fit.cutoff, out=out)
	}
	
	# update results
	res[[group]] = list(u1=u1, u2=u2, fit=fit, coefs=coefs, residuals=residuals, lab.use=lab.use, lab.fit=lab.fit, lab.con=lab.con)
    }
    res
}


full_detect_contamination = function(tpm, idents, groups, samples, fit.n=50, fit.cutoff=2, do.plot=TRUE, lab.use=NULL, prefix='test'){

    # fit models to cell groups
    cat('\n\nDetect contamination\n\n')
    cat('\nFitting group models\n')
    res.groups = detect_contamination(tpm, groups, samples, fit.n=fit.n, fit.cutoff=fit.cutoff, do.plot=do.plot, lab.use=lab.use, prefix=prefix)
    
    # calculate average model
    cat('\nModel coefficients\n')
    coefs = sapply(res.groups, function(a) a$coefs)
    print(coefs)
    cat('\nMedian coefficients\n')
    global_coefs = apply(coefs, 1, median)
    print(coefs)
    
    # run average model on cell idents
    cat('\nFitting ident models\n')
    res.idents = detect_contamination(tpm, idents, samples, global_coefs=global_coefs, fit.n=fit.n, fit.cutoff=fit.cutoff, do.plot=do.plot, lab.use=lab.use, prefix=prefix)
    
    # return all data
    return(res.groups=res.groups, res.idents=res.idents)
}
