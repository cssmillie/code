require(data.table)
source('~/code/single_cell/tpm.r')


load_mast = function(){
    library(BiocGenerics, pos=length(search()))
    library(S4Vectors, pos=length(search()))
    library(DelayedArray, pos=length(search()))
    library(MAST)
}


test_log_base = function(seur, base=2, total=1e4){
    # Test log base of seur@data
    j = sample(1:ncol(seur@data), 1)
    u = sum(base**seur@data[,j] - 1)
    (1e4 - 1e-2) <= u & u <= (1e4 + 1e-2)
}


unorder_factors = function(x){
    j = sapply(x, is.ordered)
    x[,j] = lapply(x[,j,drop=F], function(a) factor(as.character(a), levels=levels(a)))
    x
}


get_data = function(seur, data.use='tpm', tpm=NULL, cells.use=NULL){
    
    # Retrieve data from a Seurat object
    # data.use can be: counts, tpm, log2, data, any matrix
    # can pre-calculate tpm for speed
    
    if(!is.character(data.use)){
        return(data.use)
    }
    
    if(data.use == 'counts'){
        data = seur@raw.data
    }
    
    if(data.use == 'tpm'){
        if(is.null(tpm)){
	    data = calc_tpm(seur, cells.use=cells.use)
	} else {
	    data = tpm
	}
    }
    
    if(data.use == 'log2'){
        if(!test_log_base(seur, base=2, total=1e4)){
	    stop('Error: seur@data log base != 2')
        }
        data = seur@data
    }
    
    return(data)
}


select_cells = function(seur, covariates, batch.use=NULL, cells.use=NULL, max_cells=NULL){
    
    # Select cells from groups defined by covariates matrix
    # -----------------------------------------------------
    # 1. select cells.use
    # 2. build groups from covariates matrix
    # 3. select max_cells from each group
    # 4. sample evenly across batch.use
        
    cat('\nSelecting cells\n')
    
    # Subset data
    i = apply(covariates, 1, function(a) !any(is.na(a)))
    covariates = covariates[i,,drop=F]
    batch.use = batch.use[i]
    
    # Get cells to use
    cells = rownames(covariates)
    if(!is.null(cells.use)){
        cells = intersect(cells, cells.use)
    }
    
    # Construct cell groups
    j = sapply(covariates, function(a) !is.numeric(a))
    groups = as.factor(apply(covariates[cells,j,drop=F], 1, function(a) paste(a, collapse='.')))
    cat('\n')
    print(table(groups))
    
    # Select max_cells from each group
    if(!is.null(max_cells)){
    
        # Combine batch.use and groups to sample evenly across batches
        batch.use = as.factor(paste0(groups, batch.use))
	
	# Sample [max_cells] separately from each group
	cells = unname(unlist(lapply(levels(groups), function(group){
	    simple_downsample(cells=cells[groups == group], groups=batch.use[groups == group], total_cells=max_cells)
	})))
			
	# Print subsampled group sizes
	groups = apply(covariates[cells,j,drop=F], 1, function(a) paste(a, collapse='.'))
	cat('\n')
	print(table(groups))
    }
    return(cells)
}


select_genes = function(seur, stats, genes.use=NULL, min_cells=3, min_alpha=.025, min_fc=1.2, dir='both'){

    # Select genes from Seurat object by cells, genes.use, min_cells, min_alpha, and min_fc
    
    # Select genes by genes.use
    if(is.null(genes.use)){genes.use = rownames(seur@data)}
    
    # Select genes by min_cells
    g1 = as.character(stats[, (max(n) >= min_cells) | (max(ref_n) >= min_cells), .(gene)][V1 == TRUE, gene])

    # Select genes by min_alpha
    g2 = as.character(stats[, (max(alpha) >= min_alpha) | (max(ref_alpha) >= min_alpha), .(gene)][V1 == TRUE, gene])
    
    # Select genes by min_fc
    if(dir == 'pos'){
        g3 = as.character(stats[, max(log2fc) >= log2(min_fc), .(gene)][V1 == TRUE, gene])
    } else {
        g3 = as.character(stats[, max(abs(log2fc)) >= log2(min_fc), .(gene)][V1 == TRUE, gene])
    }
    
    # Intersect and return genes.use
    genes.use = Reduce(intersect, list(genes.use, g1, g2, g3))
    return(genes.use)
}


expression_stats = function(tpm, covariates, formula, lrt_regex, genes.use=NULL, cells.use=NULL, invert_method='auto', invert_logic='last', do.ratio=FALSE, samples=NULL, bg_tpm=NULL){
    
    # Calculate expression statistics for groups specified by formula
    # each column of the model matrix is either a single term A or an interaction term A:B
    # if single term, then select cells A and ~A
    # if interaction, then select cells A:B and A:~B
    # for multi-level factors, ~A is the next highest level of A
    source('~/code/util/mm_utils.r')
    
    # Select cells
    if(!is.null(cells.use)){
        tpm = tpm[,cells.use]
	covariates = covariates[cells.use,,drop=F]
    }
    
    # Select genes
    if(!is.null(genes.use)){tpm = tpm[genes.use,,drop=F]}
    total = rowSums(tpm)
    
    # Model matrix
    print('Model matrix')
    mm = as.matrix(model.matrix(as.formula(formula), data=unorder_factors(covariates)))
    
    # Invert matrix
    print('Invert matrix')
    u = mm_logical_not(mm, formula, covariates, method=invert_method, invert=invert_logic)
    
    MM = u$x
    refs = structure(u$names, names=colnames(MM))
    
    # For every column that matches lrt_regex
    print('Expression stats')
    stats = lapply(grep(lrt_regex, colnames(mm), value=T), function(a){print(a)
        
	# cell indices
	i = as.logical(mm[,a])
	j = as.logical(MM[,a])
	ref = refs[[a]]
	
	# number of expressing cells
	n1 = rowSums(tpm[,i,drop=F] > 0)
	n2 = rowSums(tpm[,j,drop=F] > 0)
	
	# total expression over all cells
	s1 = rowSums(tpm[,i,drop=F])
	s2 = rowSums(tpm[,j,drop=F])
	
	# fraction of expressing cells (alpha)
	a1 = n1/sum(i)
	a2 = n2/sum(j)
	
	# mean over expressing cells (mu)
	m1 = ifelse(n1 > 0, s1/n1, 0)
	m2 = ifelse(n2 > 0, s2/n2, 0)
	
	# mean over all cells
	u1 = s1/sum(i)
	u2 = s2/sum(j)
	
	# fraction of total expression
	t1 = s1/total
	t2 = s2/total
	
	# fix zeros for logs
	zero = .5*min(tpm[tpm > 0])/(sum(i) + sum(j))
	m1 = m1 + .5*zero
	m2 = m2 + .5*zero
	u1 = u1 + .5*zero
	u2 = u2 + .5*zero
	
	# log fold change
	log2fc = log2(u1) - log2(u2)

	if(do.ratio == TRUE){

	    # calculate ratio for group of interest
	    f1 = table(samples[colnames(tpm)[i]])
	    f1 = f1/sum(f1)
	    b1 = colSums(t(bg_tpm[,names(f1),drop=F])*as.vector(f1))	    
	    r1 = u1/b1[names(u1)]

	    # calculate ratio for reference group
	    f2 = table(samples[colnames(tpm)[j]])
	    f2 = f2/sum(f2)
	    b2 = colSums(t(bg_tpm[,names(f2),drop=F])*as.vector(f2))
	    r2 = u2/b2[names(u2)]
	    
	} else {r1 = r2 = NA}
	
	# combine in data frame
	res = data.frame(gene=rownames(tpm), contrast=a, ref=ref, n=n1, ref_n=n2, alpha=a1, ref_alpha=a2, mu=log2(m1), ref_mu=log2(m2), mean=log2(u1), ref_mean=log2(u2), total=t1, ref_total=t2, log2fc=log2fc, ratio=r1, ref_ratio=r2)
	return(res)
    })
    stats = as.data.table(do.call(rbind, stats))
    return(stats)
}


fdr_stats = function(data, covariates, formula, lrt_regex, genes.use=NULL, cells.use=NULL, invert_method='multi', invert_logic='last'){

    # Calculate FDR statistics for groups specified by formula
    # See mm_utils for information on mm_logical_not arguments
    source('~/code/util/mm_utils.r')

    # Select genes
    if(is.null(genes.use)){genes.use = rownames(data)}
    if(is.null(cells.use)){cells.use = colnames(data)}
    data = data[genes.use, cells.use]
    covariates = covariates[cells.use,,drop=F]
    
    # Model matrix
    mm = as.matrix(model.matrix(as.formula(formula), data=unorder_factors(covariates)))
    
    # Invert matrix
    u = mm_logical_not(mm, formula, covariates, method=invert_method, invert=invert_logic)
    MM = u$x
    refs = structure(u$names, names=colnames(MM))
    
    # For every column that matches lrt_regex
    stats = lapply(grep(lrt_regex, colnames(mm), value=T), function(a){
    
	# Select cells in each group
	i = as.logical(mm[,a])
	j = as.logical(MM[,a])
	ref = refs[[a]]
	
	# False discovery rate
	fdr = t(apply(data, 1, function(b){
	    predictions = b[i|j]
	    labels = ifelse(i, TRUE, NA)
	    labels[j] = FALSE
	    labels = na.omit(labels)
	    calc_fdr(predictions, factor(labels, levels=c(FALSE, TRUE)))
 	}))
	colnames(fdr) = c('cutoff', 'accuracy', 'sensitivity', 'specificity', 'fdr', 'f1')
	fdr = data.frame(gene=rownames(data), contrast=a, ref=ref, fdr)
	return(fdr)
    })
    stats = as.data.table(do.call(rbind, stats))
    return(stats)
}


p.find_markers = function(seur, ident.1=NULL, ident.2=NULL, data.use='log2', genes.use=NULL, cells.use=NULL, test.use='roc', min_cells=3, min_alpha=.05, min_fc=1.25, max_cells=1000, batch.use=NULL,
	       	          dir='pos', tpm=NULL, covariates=NULL, formula='~ ident', lrt_regex='ident', gsea.boot=100, invert_method='auto', invert_logic='last', do.stats=FALSE, n.cores=1,
			  filter_genes=TRUE, do.ratio=FALSE, samples=NULL){
    
    # Build covariates
    print(c(ident.1, ident.2))
    if(!is.null(ident.1)){
        ident.use = seur@ident
	if(is.null(ident.2)){
	    ident.use = as.factor(ifelse(ident.use == ident.1, ident.1, 'Other'))
	    ident.use = relevel(ident.use, 'Other')
	} else {
	    ident.use = factor(ifelse(ident.use %in% c(ident.1, ident.2), as.character(ident.use), NA), levels=c(ident.2, ident.1))
	}
	if(is.null(covariates)){
	    covariates = data.frame(ident=ident.use)
	} else {
	    covariates$ident = ident.use
	}
    }
    rownames(covariates) = colnames(seur@data)
    
    # Check covariates
    q = sapply(covariates, typeof)
    if('character' %in% q){print(q); stop('error: invalid covariates type')}
    
    # Select cells
    cells.use = select_cells(seur, covariates, cells.use=cells.use, max_cells=max_cells, batch.use=batch.use)
    
    # Mean TPM for contamination ratio [genes x samples]
    if(do.ratio == TRUE){
        
        # Fix samples
	if(is.null(samples)){
	    samples = as.character(seur@data.info$orig.ident)
	}
	if(is.null(names(samples))){
	    names(samples) = colnames(seur@data)
	}
	print(table(samples))
	
	# Background TPM
	bg_tpm = t(get_data(seur, data.use='tpm', tpm=tpm, cells.use=colnames(seur@data)))
	bg_tpm = sapply(unique(samples), function(a) colMeans(bg_tpm[samples == a,]))
	if(!is.null(genes.use)){bg_tpm = bg_tpm[genes.use,,drop=F]}
    
    } else {bg_tpm = NULL}
    
    # TPM for log fold changes [genes x cells]
    tpm = get_data(seur, data.use='tpm', tpm=tpm, cells.use=cells.use)
    
    # Data for DE test [genes x cells]
    data = get_data(seur, data.use=data.use, tpm=tpm, cells.use=cells.use)
    
    if(filter_genes == TRUE){
    
    # Expression stats
    stats = expression_stats(tpm, covariates, formula, lrt_regex, genes.use=genes.use, cells.use=cells.use, invert_method=invert_method, invert_logic=invert_logic,
                             do.ratio=do.ratio, samples=samples, bg_tpm=bg_tpm)
    
    # Select genes
    genes.use = select_genes(seur, stats, genes.use=genes.use, min_cells=min_cells, min_alpha=min_alpha, min_fc=min_fc, dir=dir)
    if(length(genes.use) == 0){return(c())}
    if(length(genes.use) <  5){genes.use = unique(c(genes.use, names(sort(rowSums(tpm), decreasing=T)[1:5])))}
    
    } else {
    
	stats = NULL
	genes.use = rownames(data)
    }
    
    # FDR stats
    if(do.stats == TRUE){
        fdr = fdr_stats(data, covariates, formula, lrt_regex, genes.use=genes.use, cells.use=cells.use, invert_method=invert_method, invert_logic=invert_logic)
	setkey(stats, gene, contrast, ref)
	setkey(fdr, gene, contrast, ref)
	stats = stats[fdr,]
    }
    
    # Subset data
    print(paste('Testing', length(genes.use), 'genes in', length(cells.use), 'cells'))
    data = data[genes.use, cells.use, drop=F]
    covariates = unorder_factors(covariates[cells.use, , drop=F])
    
    # Run marker tests
    labels = covariates[,1]
    if(test.use == 'f'){markers = de.rocr(data, labels, measures='f')}
    if(test.use == 'fdr'){markers = de.fdr(data, labels)}
    if(test.use == 'pr'){markers = de.rocr(data, labels, measures=c('prec', 'rec'))}
    if(test.use == 'mast'){markers = de.mast(data, covariates=covariates, formula=formula, lrt_regex=lrt_regex, n.cores=n.cores)}
    if(test.use == 'gsea'){markers = gsea.mast(data, covariates=covariates, formula=formula, lrt_regex=lrt_regex, gsea.boot=gsea.boot, n.cores=n.cores)}
    if(test.use == 'roc'){markers = de.rocr(data, labels, measures='auc')}
    if(test.use == 'mpath'){markers = de.mpath(seur@raw.data[genes.use,cells.use], covariates=covariates, formula=formula)}

    # Add cluster information
    if(! 'contrast' %in% colnames(markers)){
        markers$contrast = paste0('ident', levels(labels)[nlevels(labels)])
    }
    
    # Merge results
    markers = as.data.table(markers)
    setkey(stats, gene, contrast)
    setkey(markers, gene, contrast)
    markers = markers[stats,]
    
    # Sort markers
    if('auc' %in% colnames(markers)){
        markers = markers[order(contrast, -1*auc),]
    } else if('f1' %in% colnames(markers)){
        markers = markers[order(contrast, -1*f1),]
    } else if('pval' %in% colnames(markers)){
        markers = markers[order(contrast, pval),]
    } else if('pvalH' %in% colnames(markers)){
        markers = markers[order(contrast, pvalH),]
    } else {
        markers = markers[order(contrast, -1*alpha),]
    }

    # Add ident
    markers$ident = paste(ident.1, ident.2, sep=';')
    
    # Return marker genes
    return(markers)
}


p.find_all_markers = function(seur, data.use='log2', tpm=NULL, do.precalc=T, n.cores=1, ...){
    
    # Get cell identities
    idents = as.character(levels(seur@ident))
    
    # Pre-calculate TPM and data
    if(do.precalc == TRUE){
        tpm = get_data(seur, data.use='tpm', tpm=tpm)
        data.use = get_data(seur, data.use=data.use, tpm=tpm)
    }
    
    # Find marker genes
    run_parallel(
	foreach(i=idents, .combine=rbind) %dopar% {
	    print(i)
	    p.find_markers(seur, ident.1=i, tpm=tpm, data.use=data.use, n.cores=1, ...)
	},
	n.cores = n.cores
    )
}


p.pairwise_markers = function(seur, ident.1, data.use='log2', tpm=NULL, n.cores=1, ...){
    
    # Get cell identities
    idents = as.character(levels(seur@ident))
    
    # Pre-calculate TPM
    tpm = get_data(seur, data.use='tpm', tpm=tpm)
    
    # Pre-calculate data
    data = get_data(seur, data.use=data.use, tpm=tpm)
    
    # Find marker genes
    m = run_parallel(
        foreach(i=setdiff(idents, ident.1), .combine=rbind) %dopar% {
	    print(i)
	    p.find_markers(seur, ident.1=ident.1, ident.2=i, tpm=tpm, data.use=data, n.cores=1, ...)
	},
	n.cores = n.cores
    )
}


p.pairwise_all_markers = function(seur, data.use='log2', tpm=NULL, n.cores=1, ...){
    
    # Get cell identities
    idents = as.character(levels(seur@ident))
    
    # Pre-calculate TPM
    tpm = get_data(seur, data.use='tpm', tpm=tpm)
    
    # Pre-calculate data
    data = get_data(seur, data.use=data.use, tpm=tpm)
    
    # Find marker genes
    m = run_parallel(
        foreach(i=idents, .combine=rbind) %:% foreach(j=setdiff(idents, i), .combine=rbind) %dopar% {
	    print(c(i,j))
	    p.find_markers(seur, ident.1=i, ident.2=j, tpm=tpm, data.use=data, n.cores=1, dir='pos', ...)
	},
	n.cores = n.cores
    )
}


specific_markers = function(seur, markers, ident.use=NULL, test.use='mast', data.use='log2', tpm=NULL, n.cores=1, ...){
    
    # Check arguments
    markers[, ident := gsub(':.*', '', gsub('ident', '', contrast))]    
    markers[, keep := TRUE]
    if(! 'pvalH' %in% colnames(markers) | test.use != 'mast'){stop('test.use != mast')}
    
    # Get cell identities
    idents = as.character(levels(seur@ident))
    if(is.null(ident.use)){
        ident.use = idents
    }
    remove = setdiff(idents, unique(markers$ident))
    if(length(remove) > 0){print(paste('Removing', paste(remove, sep=', ')))}
    idents = intersect(idents, unique(markers$ident))
    
    # Pre-calculate TPM
    tpm = get_data(seur, data.use='tpm', tpm=tpm)
    
    # Pre-calculate data
    data = get_data(seur, data.use=data.use, tpm=tpm)
    
    # Find marker genes
    for(i in ident.use){

        # Sort comparisons by overlapping genes
        js = setdiff(idents, i)
	js = rev(names(sort(sapply(js, function(j) length(intersect(markers[ident == i, gene], markers[ident == j, gene]))))))
		
        for(j in js){
	    
	    # get genes to test
	    genes.1 = markers[ident == i & keep == TRUE, gene]
	    genes.2 = markers[ident == j & keep == TRUE, gene]
	    genes.use = intersect(genes.1, genes.2)
	    if(length(genes.use) == 0){next}
	    print(paste(i, j, 'testing', length(genes.use), 'genes'))	    
	    
	    # get lrt_regex
	    lrt_regex = paste(paste0(unique(markers[ident == i, contrast]), '$'), collapse='|')
	    	    
	    # run marker test
	    mi = p.find_markers(seur, ident.1=i, ident.2=j, tpm=tpm, data.use=data, genes.use=genes.use, dir='pos', test.use='mast', do.stats=F, min_fc=0, lrt_regex=lrt_regex, ...)
	    
	    # update genes to keep
	    if(is.null(mi)){
	        genes.remove = genes.use
	    } else {
	        genes.keep = mi[((pvalD <= .05 & coefD > 0) | (pvalC <= .05 & coefC > 0)) & log2fc > 0, gene]
	        genes.remove = setdiff(genes.use, genes.keep)
	    }
	    
	    print(paste('Removing', paste(genes.remove, collapse=', ')))
	    markers[ident == i & gene %in% genes.remove]$keep = FALSE
	    print(paste('Keeping', paste(markers[ident == i & keep == TRUE, gene], collapse=', ')))
	    print(paste('ident', i, 'found', sum(markers[ident == i, keep]), 'genes'))
	    print(dim(markers))
	}
    }
    markers = markers[ident %in% ident.use,]    
    return(markers)
}


calc_rocr = function(predictions, labels, measures='auc', retx=FALSE){
    require(ROCR)
    
    if(length(measures) == 1){
        q = performance(prediction(predictions, labels, label.ordering=levels(labels)), measures)
    } else {
        q = performance(prediction(predictions, labels, label.ordering=levels(labels)), measures[[1]], measures[[2]])
    }
    
    x = q@x.values
    if(length(x) > 0){x = ifelse(is.na(x[[1]]), 0, x[[1]])}
    y = ifelse(is.na(q@y.values[[1]]), 0, q@y.values[[1]])
    
    if(length(y) > 1){
        if(q@x.name == 'Cutoff'){
	    i = which.max(y[x != min(x)])
	    x = x[i]
	    y = y[i]
	} else {
            y = sum(diff(x) * (head(y,-1) + tail(y,-1)))/2
	}
    }
    
    if(retx) {c(x,y)} else {round(y, 3)}
}


de.rocr = function(data, labels, measures='auc'){
    
    # Calculate AUC of ROC curve and average difference
    scores = sapply(rownames(data), function(a){calc_rocr(as.numeric(data[a,]), labels, measures=measures)})
    
    # Return marker genes
    markers = data.frame(gene=rownames(data), stringsAsFactors=FALSE)
    cname = paste(measures, collapse='_')
    markers[,cname] = scores
    markers = markers[order(markers[,cname], decreasing=T),]
    return(markers)
}


calc_fdr = function(predictions, labels){

    # Get classification cutoff with F measure
    f = calc_rocr(as.numeric(predictions), labels, 'f', retx=T)
    
    # Get true/false negative/positives
    u = factor(as.numeric(predictions >= f[[1]]), levels=c(0,1), ordered=T)
    q = as.vector(table(u, labels))
    tn = q[[1]]; fp = q[[2]]; fn = q[[3]]; tp = q[[4]];
    
    # Calculate statistics
    fdr = fp/(tp+fp)
    tpr = tp/(tp+fn)
    tnr = tn/(tn+fp)
    acc = (tp + tn)/sum(q)
    return(c(f[[1]], acc, tpr, tnr, fdr, f[[2]]))
}


de.fdr = function(data, labels, sens_cut=.1){

    # Calculate classification statistics
    markers = t(sapply(rownames(data), function(a){calc_fdr(data[a,], labels)}))
    colnames(markers) = c('cutoff', 'accuracy', 'sensitivity', 'specificity', 'fdr', 'f1')
    
    # Calculate average difference
    avg_diff = sapply(rownames(data), function(a){as.numeric(diff(tapply(as.numeric(data[a,]), labels, mean)))})
    
    # Return marker genes
    markers = cbind(markers, data.frame(avg_diff=avg_diff, gene=rownames(markers), stringsAsFactors=F))
    markers = markers[markers$sens >= sens_cut,]
    return(markers)
    
}


de.mast = function(data, covariates, formula=NULL, lrt_regex=TRUE, n.cores=1){
    
    load_mast()
    options(mc.cores=n.cores)

    # Make single cell assay (SCA) object
    fdata = data.frame(matrix(rep(1, nrow(data))))
    covariates = as.data.frame(covariates)
    sca = MAST::FromMatrix(as.matrix(data), covariates, fdata)    
    
    # Fit MAST hurdle model
    if(is.null(formula)){
        formula = paste('~' , paste(colnames(covariates), collapse=' + '))
    }
    formula = as.formula(formula)
    zlm.obj = zlm(formula, sca)
    
    # Likelihood ratio test
    if(is.logical(lrt_regex)){lrt_regex = colnames(covariates)}
    contrasts = grep(paste(lrt_regex, collapse='|'), colnames(zlm.obj@coefC), value=T, perl=T)
    print(contrasts)
    res = summary(zlm.obj, doLRT=contrasts)$datatable
    
    # Get component information
    res.f = res[res$component == 'logFC', .(primerid, contrast, coef)]
    res.d = res[res$component == 'D', .(primerid, contrast, coef, `Pr(>Chisq)`)]
    res.c = res[res$component == 'C', .(primerid, contrast, coef, `Pr(>Chisq)`)]
    res.h = res[res$component == 'H', .(primerid, contrast, `Pr(>Chisq)`)]
    
    # Combine results
    res = merge(res.d, res.c, by=c('primerid', 'contrast'), all=T, suffixes=c('D', 'C'))
    res = Reduce(function(...) merge(..., by=c('primerid', 'contrast'), all=T), list(res, res.f, res.h))
    res = data.frame(subset(res, !is.na(`Pr(>Chisq)`)), stringsAsFactors=F)
    
    # Cleanup results
    colnames(res) = c('gene', 'contrast', 'coefD', 'pvalD', 'coefC', 'pvalC', 'mastfc', 'pvalH')
    res = res[order(res$contrast, res$pvalH),]

    # Replace NAs in mastfc with maximum
    res = as.data.table(res)
    res[, mastfc := ifelse(is.na(mastfc), max(mastfc, na.rm=T), mastfc), .(contrast)]
    
    # Adjust p-values
    res[, padjD := p.adjust(pvalD, 'fdr'), .(contrast)]
    res[, padjC := p.adjust(pvalC, 'fdr'), .(contrast)]
    res[, padjH := p.adjust(pvalH, 'fdr'), .(contrast)]
    
    options(mc.cores=1)
    return(res)
}


filter_mast = function(markers, regex=NULL, pval=NULL, padj=NULL, min_log2fc=NULL, type='all', seur=NULL, covariates=NULL, formula=NULL, top=NULL, sort.by=pval, split=FALSE){

    # Select by contrast
    if(!is.null(regex)){markers = markers[grep(regex, contrast)]}
    
    # Adjust p-values
    if(! 'padj' %in% colnames(markers)){markers[, padj := p.adjust(pval, 'fdr'), .(contrast)]}
    
    # Filter by p-value
    if(!is.null(pval)){markers = markers[pval <= pval]}
    if(!is.null(padj)){markers = markers[padj <= padj]}
    
    # Filter by log2fc
    if(!is.null(min_log2fc)){markers = markers[log2fc >= min_log2fc]}
    
    # Filter by type
    if(type == 'all'){
        markers = markers
    } else if(type == 'unique'){
        markers = markers[, if(.N == 1){.SD}, .(gene)]
    } else {
        markers = specific_markers(seur, markers, covariates=covariates, formula=formula)[keep == TRUE]
    }
    
    # Sort and select top markers
    i = order(markers[,contrast], with(sort.by, markers))
    markers = markers[i]
    if(!is.null(top)){markers = markers[,.SD[1:min(.N, top)],.(contrast)]}
    
    # Split markers by contrast
    if(split == TRUE){markers = split(markers, markers$contrast)}
    
    return(markers)
}


collapse_markers = function(m){

    # Collapse markers on (contrast, gene) using mean or [[1]]
    
    # track colnames
    a = colnames(m)
    
    # split columns by numeric
    j = sapply(m, is.numeric)
    u = m[, c('contrast', 'gene', names(which(j))), with=F]
    v = m[, names(which(!j)), with=F]

    # mean or [[1]] columns
    u = u[, lapply(.SD, mean), .(contrast, gene)]
    v = v[, .SD[1,], .(contrast, gene)]

    # merge lists
    setkey(u, contrast, gene)
    setkey(v, contrast, gene)
    m = u[v][,a,with=F]

    # order by contrast and pvalH
    m[order(contrast, pvalH)]
}

next_stats = function(markers, min_pval=1, min_padj=1, anno=NULL){

    # for each ident/gene, calculate next statistics using all non-overlapping idents
    # non-overlapping idents calculated using select_keys(anno, i, invert=T)
    # example: load_anno(key='name', value='ident')
    # optionally filter the marker matrix to keep only significant markers
    
    # filter marker lists (optionally)
    markers = markers[(pvalD <= min_pval | pvalC <= min_pval)]
    markers = markers[(padjD <= min_pval | padjC <= min_pval)]
    
    # make and test annotations
    idents = sort(unique(markers$ident))
    if(is.null(anno)){anno = structure(idents, names=idents)}
    if(any(! idents %in% names(anno))){stop('ident not in names(anno)')}
    
    # initialize columns
    markers[, c('next_alpha', 'next_mu', 'next_mean', 'next_log2fc') := as.numeric(NA)]
    setkeyv(markers, c('ident', 'gene'))
    # test if keys are unique
    if(nrow(markers) != uniqueN(markers[,.(ident, gene)])){stop('(ident, gene) not unique')}
    
    # iterate over idents
    idents = sort(unique(markers$ident))
    for(i in idents){
        j = select_keys(anno, i, invert=T)
	m = markers[ident %in% j][, list(next_alpha=max(alpha), next_mu=max(mu), next_mean=max(mean)), .(gene)]
	markers[.(i, m$gene), c('next_alpha', 'next_mu', 'next_mean') := m[,.(next_alpha, next_mu, next_mean)]]
    }
    markers$next_log2fc = markers$mean - markers$next_mean
    return(markers)
}
