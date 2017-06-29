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


select_cells = function(seur, ident.1, ident.2=NULL, cells.use=NULL, max_cells=NULL){
    
    # Select cells from Seurat object by identity, cells.use, and max_cells
    
    # Select cells by ident
    cells.1 = na.omit(colnames(seur@data)[as.character(seur@ident) == as.character(ident.1)])
    if(is.null(ident.2)){
        cells.2 = na.omit(colnames(seur@data)[as.character(seur@ident) != as.character(ident.1)])
    } else {
        cells.2 = na.omit(colnames(seur@data)[as.character(seur@ident) == as.character(ident.2)])
    }
    
    # Intersect with cells.use
    if(!is.null(cells.use)){
        cells.1 = intersect(cells.1, cells.use)
        cells.2 = intersect(cells.2, cells.use)
    }

    # Subsample by max_cells
    if(!is.null(max_cells)){
        cells.1 = sample(cells.1, min(max_cells, length(cells.1)))
        cells.2 = sample(cells.2, min(max_cells, length(cells.2)))
    }
    
    return(list(cells.1=cells.1, cells.2=cells.2))
}


select_genes = function(seur, cells.1, cells.2, tpm=NULL, data=NULL, genes.use=NULL, min_cells=3, min_pct=.05, fc.use='tpm', min_fc=1.25, dir='both'){

    # Select genes from Seurat object by cells.1, cells.2, genes.use, min_cells, min_pct, and min_fc
    # fc.use = 'tpm' or 'data' (for the latter, 'data' must be supplied)

    # Intersect genes with genes.use
    if(is.null(genes.use)){genes.use = rownames(seur@data)}
    genes.use = intersect(genes.use, rownames(seur@data))

    # Calculate number of cells
    ncells.1 = rowSums(seur@raw.data[genes.use, cells.1, drop=F] > 0)
    ncells.2 = rowSums(seur@raw.data[genes.use, cells.2, drop=F] > 0)

    # Select genes by min_cells
    i1 = (ncells.1 >= min_cells) | (ncells.2 >= min_cells)
    
    # Select genes by min_pct
    i2 = (ncells.1/length(cells.1) >= min_pct) | (ncells.2/length(cells.2) >= min_pct)
    
    # Calculate log fold change
    genes.use = genes.use[i1 & i2]
    if(fc.use == 'tpm'){
        u = rowMeans(tpm[genes.use, cells.1, drop=F])
	v = rowMeans(tpm[genes.use, cells.2, drop=F])
	logfc = log2(u + 1) - log2(v + 1)
    } else if(fc.use == 'data') {
        u = rowMeans(data[genes.use, cells.1, drop=F])
	v = rowMeans(data[genes.use, cells.2, drop=F])
	logfc = log2(u - min(u) + 1) - log2(v - min(v) + 1)
    } else {
        stop('Invalid argument: fc.use not in [tpm, data]')
    }
    
    # Select by log fold change
    if(min_fc < 1){stop('Invalid argument: min_fc < 1')}
    if(dir == 'pos'){
        i3 = logfc >= log2(min_fc)
    } else if(dir == 'neg'){
        i3 = logfc <= -1*log2(min_fc)
    } else if(dir == 'both'){
        i3 = abs(logfc) >= log2(min_fc)
    } else {
        stop('Invalid argument: dir not in [pos, neg, both]')
    }
    genes.use = genes.use[i3]
    
    return(list(genes.use=genes.use, logfc=logfc))    
}


p.find_markers = function(seur, ident.1, ident.2=NULL, data.use='log2', genes.use=NULL, cells.use=NULL, test.use='roc', min_cells=3, min_pct=.05, min_fc=1.25, max_cells=1000, dir='pos',
	                  fc.use='tpm', tpm=NULL, covariates=NULL, formula='~ label', lrt_regex='label', gsea.boot=100, n.cores=1){
    
    # Select cells
    cells = select_cells(seur, ident.1, ident.2=ident.2, cells.use=cells.use, max_cells=max_cells)
    cells.1 = cells$cells.1
    cells.2 = cells$cells.2
    cells.use = c(cells.1, cells.2)
    
    # Fix variables
    if(length(cells.1) <= 5 | length(cells.2) <= 5){return(c())}
    if(is.null(ident.2)){ident.2 = 'other'}
    
    # TPM for log fold changes
    if(fc.use == 'tpm'){tpm = get_data(seur, data.use='tpm', tpm=tpm)}
    
    # Data for DE test
    data = get_data(seur, data.use=data.use, tpm=tpm)
    
    # Cell identities
    labels = factor(c(rep(ident.1, length(cells.1)), rep(ident.2, length(cells.2))), levels=c(ident.2, ident.1))
    
    # Select genes
    genes = select_genes(seur, cells.1, cells.2, tpm=tpm, data=data, genes.use=genes.use, min_cells=min_cells, min_pct=min_pct, fc.use=fc.use, min_fc=min_fc, dir=dir)
    genes.use = genes$genes.use
    if(length(genes.use) <= 1){return(c())}
    
    # Subset data
    print(paste('Testing', length(genes.use), 'genes in', length(cells.use), 'cells'))
    data = data[genes.use, cells.use, drop=F]

    # Build covariates
    if(is.null(covariates)){
        covariates = data.frame(label=labels)
	rownames(covariates) = cells.use
    } else {
        covariates = covariates[cells.use,]
	covariates$label = labels
    }
    
    # Run marker tests
    if(test.use == 'f'){markers = de.rocr(data, labels, measures='f')}
    if(test.use == 'fdr'){markers = de.fdr(data, labels)}
    if(test.use == 'pr'){markers = de.rocr(data, labels, measures=c('prec', 'rec'))}
    if(test.use == 'mast'){markers = de.mast(data, covariates=covariates, formula=formula, lrt_regex=lrt_regex, n.cores=n.cores)}
    if(test.use == 'gsea'){markers = gsea.mast(data, covariates=covariates, formula=formula, lrt_regex=lrt_regex, gsea.boot=gsea.boot, n.cores=n.cores)}
    if(test.use == 'roc'){markers = de.rocr(data, labels, measures='auc')}
    
    # Add cluster information
    markers$cluster = ident.1
    markers$ref_cluster = ident.2

    # Log fold change
    markers$log2fc = genes$logfc[markers$gene]
    
    # Return marker genes
    return(markers)
}


p.find_all_markers = function(seur, data.use='log2', genes.use=NULL, cells.use=NULL, test.use='roc', min_cells=3, min_pct=.05, min_fc=1.25, max_cells=1000, dir='pos',
		              fc.use='tpm', tpm=NULL, covariates=NULL, formula=NULL, lrt_regex='label', gsea.boot=100, n.cores=1){
    # Get cell identities
    idents = as.character(levels(seur@ident))
    
    # Pre-calculate TPM
    tpm = get_data(seur, data.use='tpm', tpm=tpm)
    
    # Pre-calculate data
    data = get_data(seur, data.use=data.use, tpm=tpm)
    
    # Get functions to export
    export = c('p.find_markers', 'get_data', 'select_cells', 'select_genes', 'calc_rocr', 'de.rocr', 'calc_fdr', 'de.fdr', 'de.mast', 'gsea.mast', 'format_mast_output')
    
    # Find marker genes
    run_parallel(
        foreach(i=idents, .export=export, .combine=rbind) %dopar% {
	    print(i)
	    m = p.find_markers(seur, ident.1=i, tpm=tpm, data.use=data, genes.use=genes.use, cells.use=cells.use, test.use=test.use, min_cells=min_cells, min_pct=min_pct, min_fc=min_fc, max_cells=max_cells, dir=dir, fc.use=fc.use, covariates=covariates, formula=formula, lrt_regex=lrt_regex, gsea.boot=gsea.boot, n.cores=n.cores)
	    return(m)
	},
	n.cores = n.cores
    )
}


p.pairwise_markers = function(seur, ident.1, data.use='log2', genes.use=NULL, cells.use=NULL, test.use='roc', min_cells=3, min_pct=.05, min_fc=1.25, max_cells=1000, dir='pos',
		              fc.use='tpm', tpm=NULL, covariates=NULL, formula=NULL, lrt_regex='label', gsea.boot=100, n.cores=1){
    # Get cell identities
    idents = as.character(levels(seur@ident))
    
    # Pre-calculate TPM
    tpm = get_data(seur, data.use='tpm', tpm=tpm)
    
    # Pre-calculate data
    data = get_data(seur, data.use=data.use, tpm=tpm)
    
    # Get functions to export
    export = c('p.find_markers', 'get_data', 'select_cells', 'select_genes', 'calc_rocr', 'de.rocr', 'calc_fdr', 'de.fdr', 'de.mast', 'gsea.mast', 'format_mast_output')
    
    # Find marker genes
    m = run_parallel(
        foreach(i=setdiff(idents, ident.1), .export=export, .combine=rbind) %dopar% {
	    print(c(ident.1, i))
	    m = p.find_markers(seur, ident.1=ident.1, ident.2=i, tpm=tpm, data.use=data, genes.use=genes.use, cells.use=cells.use, test.use=test.use, min_cells=min_cells, min_pct=min_pct, min_fc=min_fc, max_cells=max_cells, dir=dir, fc.use=fc.use, covariates=covariates, formula=formula, lrt_regex=lrt_regex, gsea.boot=gsea.boot, n.cores=n.cores)
	    if(!is.null(m)){m$ref_cluster = i}
	    return(m)
	},
	n.cores = n.cores
    )
}


p.pairwise_all_markers = function(seur, data.use='log2', genes.use=NULL, cells.use=NULL, test.use='roc', min_cells=3, min_pct=.05, min_fc=1.25, max_cells=1000, dir='pos',
		              fc.use='tpm', tpm=NULL, covariates=NULL, formula=NULL, lrt_regex='label', gsea.boot=100, n.cores=1){
    # Get cell identities
    idents = as.character(levels(seur@ident))
    
    # Pre-calculate TPM
    tpm = get_data(seur, data.use='tpm', tpm=tpm)
    
    # Pre-calculate data
    data = get_data(seur, data.use=data.use, tpm=tpm)
    
    # Get functions to export
    export = c('p.find_markers', 'get_data', 'select_cells', 'select_genes', 'calc_rocr', 'de.rocr', 'calc_fdr', 'de.fdr', 'de.mast', 'gsea.mast', 'format_mast_output')
    
    # Find marker genes
    m = run_parallel(
        foreach(i=idents, .combine=rbind) %:% foreach(j=setdiff(idents, i), .export=export, .combine=rbind) %dopar% {
	    print(c(i,j))
	    m = p.find_markers(seur, ident.1=i, ident.2=j, tpm=tpm, data.use=data, genes.use=genes.use, cells.use=cells.use, test.use=test.use, min_cells=min_cells, min_pct=min_pct, min_fc=min_fc, max_cells=max_cells, dir=dir, fc.use=fc.use, covariates=covariates, formula=formula, lrt_regex=lrt_regex, gsea.boot=gsea.boot, n.cores=n.cores)
	    if(!is.null(m)){m$ref_cluster = j}
	    return(m)
	},
	n.cores = n.cores
    )
}


merge_markers = function(fns=NULL, path='.', pattern=NULL, order_by=NULL, cluster_regex=NULL, ref_regex=NULL){

    # Combine markers with option to extract cluster names with regex
    library(data.table)
    
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
    
    library(data.table)
    group_by = substitute(group_by)
    collapse_by = substitute(collapse_by)
    
    m = setDT(m)[,{mi = .SD[eval(collapse_by)]; mi$N = .N; mi}, group_by]
    return(as.data.frame(m))
}


calc_rocr = function(predictions, labels, measures='auc', retx=FALSE){
    library(ROCR)
    
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
    return(res)
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
	    u = factor(all_genes %in% i, levels=c(FALSE, TRUE))
	    v = factor(all_genes %in% j, levels=c(FALSE, TRUE))
	    fisher.test(table(u,v))$pval
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
        formula = as.formula(paste('~' , paste(colnames(covariates), collapse=' + ')))
    }
    zlm.obj = zlm(formula, sca)
    
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
    
    library(ggplot2)
    library(tidyr)
    library(naturalsort)
	     
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

    library(topGO)
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

    library(topGO)
    library(naturalsort)
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
    seur@data.info$nGene = colSums(seur@raw.data > 0)
    seur@data.info$nUMI = colSums(seur@raw.data)
    
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
