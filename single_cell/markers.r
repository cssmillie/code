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


get_data = function(seur, data.use='tpm', tpm=NULL){
    
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
	    data = calc_tpm(seur)
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


select_cells = function(covariates, cells.use=NULL, max_cells=NULL){
    
    # Select cells from covariates matrix using cells.use and max_cells
    cat('\nSelecting cells\n')
    
    # Remove cells with missing data
    cells = rownames(covariates)[apply(covariates, 1, function(a) !any(is.na(a)))]
    
    # Intersect cells with cells.use
    if(!is.null(cells.use)){
        cells = intersect(cells, cells.use)
    }
    
    # Construct cell groups
    j = sapply(covariates, function(a) !is.numeric(a))
    groups = apply(covariates[cells,j,drop=F], 1, function(a) paste(a, collapse='.'))
    cat('\n')
    print(table(groups))
    
    # Select max_cells from each group
    if(!is.null(max_cells)){
        cells = unlist(tapply(cells, groups, function(a) sample(a, min(max_cells, length(a)))))
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


expression_stats = function(tpm, covariates, formula, lrt_regex, cells.use=NULL, invert_method='multi', invert_logic='last'){
    
    # Calculate expression statistics for groups specified by formula
    # each column of the model matrix is either a single term A or an interaction term A:B
    # if single term, then select cells A and ~A
    # if interaction, then select cells A:B and A:~B
    # for multi-level factors, ~A is the next highest level of A
    source('~/code/util/mm_utils.r')
    
    # Select cells
    total = rowSums(tpm)
    if(!is.null(cells.use)){
        tpm = tpm[,cells.use]
	covariates = covariates[cells.use,,drop=F]
    }
    
    # Model matrix
    print('Model matrix')
    mm = as.matrix(model.matrix(as.formula(formula), data=covariates))
    
    # Invert matrix
    print('Invert matrix')
    u = mm_logical_not(mm, formula, covariates, method=invert_method, invert=invert_logic)
    MM = u$x
    refs = structure(u$names, names=colnames(MM))
    
    # For every column that matches lrt_regex
    print('Expression stats')
    stats = lapply(grep(lrt_regex, colnames(mm), value=T), function(a){
        
	# Select cells in each group
	i = as.logical(mm[,a])
	j = as.logical(MM[,a])
	ref = refs[[a]]
	
	# Gene expression statistics
	n1 = rowSums(tpm[,i] > 0); n2 = rowSums(tpm[,j] > 0)
	s1 = rowSums(tpm[,i]); s2 = rowSums(tpm[,j])
	a1 = n1/sum(i); a2 = n2/sum(j)
	m1 = s1/n1; m2 = s2/n2
	u1 = s1/sum(i); u2 = s2/sum(j)
	t1 = s1/total; t2 = s2/total
	zero = min(min(u1[u1 > 0]), min(u2[u2 > 0]))
	log2fc = log2(u1 + .5*zero) - log2(u2 + .5*zero)
	res = data.frame(gene=rownames(tpm), contrast=a, ref=ref, n=n1, ref_n=n2, alpha=a1, ref_alpha=a2, mu=m1, ref_mu=m2, mean=log2(u1), ref_mean=log2(u2), total=t1, ref_total=t2, log2fc=log2fc)
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
    mm = as.matrix(model.matrix(as.formula(formula), data=covariates))
    
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


p.find_markers = function(seur, ident.1=NULL, ident.2=NULL, data.use='log2', genes.use=NULL, cells.use=NULL, test.use='roc', min_cells=3, min_alpha=.05, min_fc=1.25, max_cells=1000, dir='pos',
	                  tpm=NULL, covariates=NULL, formula='~ ident', lrt_regex='ident', gsea.boot=100, invert_method='multi', invert_logic='last', do.stats=FALSE, n.cores=1){
    
    # Build covariates
    print(c(ident.1, ident.2))
    if(!is.null(ident.1)){
        ident.use = seur@ident
	if(is.null(ident.2)){
	    ident.use = as.factor(ifelse(ident.use == ident.1, ident.1, 'Other'))
	    ident.use = relevel(ident.use, 'Other')
	} else {
	    ident.use = as.factor(ifelse(ident.use %in% c(ident.1, ident.2), as.character(ident.use), NA))
	    ident.use = relevel(ident.use, ident.2)
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
    cells.use = select_cells(covariates, cells.use=cells.use, max_cells=max_cells)
    
    # TPM for log fold changes
    tpm = get_data(seur, data.use='tpm', tpm=tpm)
    
    # Data for DE test
    data = get_data(seur, data.use=data.use, tpm=tpm)
    
    # Expression stats
    stats = expression_stats(tpm, covariates, formula, lrt_regex, cells.use=cells.use, invert_method=invert_method, invert_logic=invert_logic)
    
    # Select genes
    genes.use = select_genes(seur, stats, genes.use=genes.use, min_cells=min_cells, min_alpha=min_alpha, min_fc=min_fc, dir=dir)
    if(length(genes.use) <= 1){return(c())}
    
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
    covariates = covariates[cells.use, , drop=F]
    print(lapply(covariates, table))
    
    # Run marker tests
    labels = covariates[,1]
    if(test.use == 'f'){markers = de.rocr(data, labels, measures='f')}
    if(test.use == 'fdr'){markers = de.fdr(data, labels)}
    if(test.use == 'pr'){markers = de.rocr(data, labels, measures=c('prec', 'rec'))}
    if(test.use == 'mast'){markers = de.mast(data, covariates=covariates, formula=formula, lrt_regex=lrt_regex, n.cores=n.cores)}
    if(test.use == 'gsea'){markers = gsea.mast(data, covariates=covariates, formula=formula, lrt_regex=lrt_regex, gsea.boot=gsea.boot, n.cores=n.cores)}
    if(test.use == 'roc'){markers = de.rocr(data, labels, measures='auc')}
    
    # Add cluster information
    if(! 'contrast' %in% colnames(markers)){
        markers$contrast = paste0('ident', levels(labels)[nlevels(labels)])
    }
    
    # Merge results
    markers = as.data.table(markers)
    setkey(stats, gene, contrast)
    setkey(markers, gene, contrast)
    markers = markers[stats,,nomatch=0]
    
    # Sort markers
    if('auc' %in% colnames(markers)){
        markers = markers[order(contrast, -1*auc),]
    } else if('f1' %in% colnames(markers)){
        markers = markers[order(contrast, -1*f1),]
    } else if('pval' %in% colnames(markers)){
        markers = markers[order(contrast, pval),]
    } else {
        markers = markers[order(contrast, -1*alpha),]
    }
    
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
    library(pryr)
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


merge_markers = function(fns=NULL, path='.', pattern=NULL, order_by=NULL, cluster_regex=NULL, ref_regex=NULL){

    # Combine markers with option to extract cluster names with regex
    require(data.table)
    
    # Read files
    if(is.null(fns)){
        fns = list.files(path=path, pattern=pattern)
    }
    print(paste('Merging', length(fns), 'files'))

    # Test regex
    if(!is.null(cluster_regex)){print(paste('Cluster regex:', gsub(cluster_regex, '\\1', fns[[1]], perl=T)))}
    if(!is.null(ref_regex)){print(paste('Ref regex:', gsub(ref_regex, '\\1', fns[[1]], perl=T)))}

    # Get markers
    m = do.call(rbind, lapply(fns, function(fn){
	mi = as.data.frame(fread(fn, sep='\t', header=T, stringsAsFactors=F))
	if(!is.null(cluster_regex)){mi[['cluster']] = gsub(cluster_regex, '\\1', fn)}
	if(!is.null(ref_regex)){mi[['ref_cluster']] = gsub(ref_regex, '\\1', fn)}
	mi
    }))
    return(m)
}


collapse_markers = function(m, group_by=list(cluster, gene), collapse_by=which.max(mastfc), strict=FALSE){

    # Collapse pairwise marker list
    # Uses non-standard evaluation (see default arguments)
    
    require(data.table)
    group_by = substitute(group_by)
    collapse_by = substitute(collapse_by)
    
    m = setDT(m)[,{mi = .SD[eval(collapse_by)]; mi$N = .N; mi}, group_by]
    return(as.data.frame(m))
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
    res.d = res[res$component == 'D', .(primerid, contrast, coef)]
    res.c = res[res$component == 'C', .(primerid, contrast, coef)]
    res.h = res[res$component == 'H', .(primerid, contrast, `Pr(>Chisq)`)]
    
    # Combine results
    res = merge(res.d, res.c, by=c('primerid', 'contrast'), all=T, suffixes=c('D', 'C'))
    res = Reduce(function(...) merge(..., by=c('primerid', 'contrast'), all=T), list(res, res.f, res.h))
    res = data.frame(subset(res, !is.na(`Pr(>Chisq)`)), stringsAsFactors=F)
    
    # Cleanup results
    colnames(res) = c('gene', 'contrast', 'coefD', 'coefC', 'mastfc', 'pval')
    res = res[order(res$contrast, res$pval),]
    
    options(mc.cores=1)
    return(as.data.table(res))
}


gsea.fisher = function(gene_set1, gene_set2, reference='~/aviv/db/map_gene/hg19_genes.txt', n.cores=1){

    # GSEA with fisher test
    # gene_set1 = list of genes
    # gene_set2 = list of genes
    # return matrix of p-values

    # convert gene sets to list
    if(typeof(gene_set1) != 'list'){gene_set1 = list(gene_set1)}
    if(typeof(gene_set2) != 'list'){gene_set2 = list(gene_set2)}
    
    # pairwise fisher tests
    all_genes = readLines(reference)
    m = run_parallel(
        foreach(i=gene_set1, .combine=rbind) %:% foreach(j=gene_set2, .combine=c) %dopar% {
	    u = factor(all_genes %in% unlist(i), levels=c(FALSE, TRUE))
	    v = factor(all_genes %in% unlist(j), levels=c(FALSE, TRUE))
	    fisher.test(table(u,v))$p.value
	},
	n.cores=n.cores
    )
    rownames(m) = names(gene_set1)
    colnames(m) = names(gene_set2)
    return(m)
}

gsea.mast = function(data, covariates, formula=NULL, lrt_regex=TRUE, gsea.boot=100, n.cores=1){

    load_mast()
    options(mc.cores=n.cores)
    
    # Make single cell assay object
    fdata = data.frame(matrix(rep(1, nrow(data))))
    covariates = as.data.frame(covariates)
    sca = MAST::FromMatrix(as.matrix(data), covariates, fdata)    
    
    # Fit MAST hurdle model
    if(is.null(formula)){
        formula = paste('~' , paste(colnames(covariates), collapse=' + '))
    }
    zlm.obj = zlm(as.formula(formula), sca)
    
    # Calculate bootstraps
    boot = bootVcov1(zlm.obj, gsea.boot)
    
    # Get GO gene list
    genes = read.table('~/aviv/db/gopca/go_annotations_human.tsv', sep='\t', stringsAsFactors=F, row.names=1)
    sets = structure(strsplit(genes[,4], ','), names=rownames(genes))
    sets = lapply(sets, function(a){b = as.integer(na.omit(match(a, rownames(zlm.obj@coefC))))})
    sets = sets[sapply(sets, length) >= 5]
    
    # Get hypothesis columns
    names = colnames(zlm.obj@coefC)[grep(paste(lrt_regex, collapse='|'), colnames(zlm.obj@coefC))]
    
    # Perform GSEA
    res = lapply(names, function(a){
        gsea = summary(gseaAfterBoot(zlm.obj, boot, sets, CoefficientHypothesis(a)))
	gsea$name = genes[as.character(gsea$set), 3]
	gsea$genes = genes[as.character(gsea$set), 4]
	gsea$contrast = a
	return(gsea)
    })
    res = do.call(rbind, res)
    
    options(mc.cores=1)
    return(res)
}


plot_markers = function(markers, top=25, dir='pos', cols=1, base_size=12, regex=NULL, rm_regex=NULL){
    
    require(ggplot2)
    require(tidyr)
    require(naturalsort)
	     
    single_plot = function(m, cluster='cluster', gene='gene', value='auc', palette='OrRd', rm_regex=NULL){

        u = data.frame(tapply(1:nrow(m), list(m[,cluster]), function(a){b=a[1:top]; m[b,gene,drop=F]}))
	v = data.frame(tapply(1:nrow(m), list(m[,cluster]), function(a){b=a[1:top]; m[b,value,drop=F]}))
	
	if(!is.null(regex)){
	    u = u[,grepl(regex, colnames(u))]
	    v = v[,grepl(regex, colnames(v))]
	}
	
	if(!is.null(rm_regex)){
	    u = u[,!grepl(rm_regex, colnames(u))]
	    v = v[,!grepl(rm_regex, colnames(v))]
	}
	
	if(value == 'pval'){v = -1*log10(v)}
	
	u$index = rev(1:nrow(u))
	v$index = rev(1:nrow(v))
	
	u = gather(u, Cluster, Gene, -index)
	v = gather(v, Cluster, Value, -index)

	# Fix some formatting
	v$Value[is.infinite(v$Value)] = max(v$Value[is.finite(v$Value)])	
	u$Cluster = factor(u$Cluster, levels=naturalsort(unique(u$Cluster)), ordered=T)
	v$Cluster = factor(v$Cluster, levels=naturalsort(unique(v$Cluster)), ordered=T)
	legend_title = list('auc'='AUC', 'pval'='-log10(pval)')[value]
	
	p = ggplot(v) + geom_tile(aes(x=Cluster, y=index, fill=Value)) + scale_fill_distiller(palette=palette, trans='reverse', na.value='white') + geom_text(data=u, aes(x=Cluster, y=index, label=Gene)) + ylab('Gene') + xlab('Covariate') + guides(fill=guide_legend(title=legend_title)) + theme_minimal(base_size=base_size)
	
	return(p)
    }
    
    if('contrast' %in% colnames(markers)){
        cluster = 'contrast'
	value = 'pval'
	m.pos = markers[markers$log2fc > 0,]
	m.neg = markers[markers$log2fc < 0,]
    } else if('auc' %in% colnames(markers)) {
        cluster = 'cluster'
	value = 'auc'
	m.pos = markers[markers$auc >= .5,]
	m.neg = markers[markers$auc < .5,]
	m.neg[,'auc'] = 1 - m.neg[,'auc']
    } else if('prec_rec' %in% colnames(markers)){
        cluster = 'cluster'
	value = 'prec_rec'
	m.pos = markers[markers$prec_rec >= .5,]
	m.neg = markers[markers$prec_rec < .5,]
	m.neg['prec_rec'] = 1 - m.neg[,'prec_rec']
    } else {
        print('Unknown marker type')
    }
    
    if(dir == 'pos'){
        single_plot(m.pos, cluster=cluster, gene='gene', value=value, palette='OrRd', rm_regex=rm_regex)
    } else if(dir == 'neg'){
        single_plot(m.neg, cluster=cluster, gene='gene', value=value, palette='Blues', rm_regex=rm_regex)
    } else if(dir == 'both'){
        ps = list()
        ps[['pos']] = single_plot(m.pos, cluster=cluster, gene='gene', value=value, palette='OrRd', rm_regex=rm_regex)
        ps[['neg']] = single_plot(m.neg, cluster=cluster, gene='gene', value=value, palette='Blues', rm_regex=rm_regex)
    	
        source('~/code/single_cell/plot.r')
        multiplot(plotlist=ps, cols=cols)
    }

}


go_genes = function(seur, genes, ontology='BP'){

    require(topGO)
    all_genes = rownames(seur@data)
    GO2genes = readMappings(file='~/aviv/db/gopca/go.test.txt', sep='\t')
    gene_list = as.numeric(all_genes %in% genes)
    names(gene_list) = all_genes
    gene_list = factor(gene_list)
    
    # Run topGO tests
    GOdata = new('topGOdata', ontology=ontology, allGenes=gene_list, annot=annFUN.GO2genes, GO2genes=GO2genes)
    res = runTest(GOdata, algorithm='classic', statistic='ks')
    res = GenTable(GOdata, res, topNodes=min(1000, length(res@score)))
    res = res[res$result1 <= .05,]
    return(res)
}


go_markers = function(m, top=NULL, pval=NULL, auc=NULL, ontology='BP', n.cores=1){

    require(topGO)
    require(naturalsort)
    GO2genes = readMappings(file='~/aviv/db/gopca/go.test.txt', sep='\t')
    
    # Ontology can be: BP (biological process), MF (molecular function), or CC (cellular component)

    # Get column names
    if('contrast' %in% colnames(m)){
	score = 'pval'
	cluster = 'contrast'
	m = m[m$log2fc > 0,]
    } else {
	score = 'auc'
	cluster = 'cluster'
	m = m[m$auc > .5,]
    }
    
    # Test each cluster
    all_genes = unique(m$gene)
    clusters = naturalsort(as.character(unique(m[,cluster])))
    
    go_terms = run_parallel(foreach(a=clusters, .packages=c('topGO')) %dopar% {

        # Filter marker genes
        mi = m[m[,cluster] == a,]
	if(!is.null(top)){
	    print('Using top genes')
	    mi = mi[1:top,]
	}
	if(!is.null(pval)){
	    print('Filtering by pval')
	    mi = mi[mi$pval <= pval,]
	}
	if(!is.null(auc)){
	    print('Filtering by AUC')
	    mi = mi[mi$auc <= auc,]
	}
	
	# Construct gene list
	gene_list = as.numeric(all_genes %in% mi$gene)
	names(gene_list) = all_genes
	gene_list = factor(gene_list)

	# Run topGO tests
	GOdata = new('topGOdata', ontology=ontology, allGenes=gene_list, annot=annFUN.GO2genes, GO2genes=GO2genes)
	res = runTest(GOdata, algorithm='classic', statistic='ks')
	res = GenTable(GOdata, res, topNodes=min(1000, length(res@score)))
	res = res[res$result1 <= .05,]
	return(res)
    }, n.cores=n.cores)
    
    names(go_terms) = clusters
    return(go_terms)
}

update_signatures = function(seur){    
    
    # Predict host (human or mouse)
    genes = rownames(seur@data)
    
    # Calculate signatures
    seur@data.info = cbind(seur@data.info, get_scores(seur, file='cell_cycle'))
    seur@data.info = cbind(seur@data.info, get_scores(seur, file='early'))
    seur@data.info$nGene = colSums(seur@data > 0)
    
    return(seur)
}

marker_table = function(mlist, pval=1, padj=1, top=100){
    mlist$padj = p.adjust(mlist$pval, method='fdr')
    if(! 'cluster' %in% colnames(mlist)){
    	 print('Could not find cluster column in marker list')
    	 if('contrast' %in% colnames(mlist)){
	     print('Using contrast column instead')
	     mlist$cluster = mlist$contrast
	 } else {stop()}
    }
    mlist = mlist[mlist$pval <= pval & mlist$padj <= padj,]
    mlist = mlist[order(mlist$cluster, mlist$pval),]
    mtab = tapply(mlist$gene, list(mlist$cluster), function(a){
        as.character(a[1:top])
    })    
    mtab = sapply(mtab, function(a){if(!is.null(a)){a}})    
    return(mtab)
}
