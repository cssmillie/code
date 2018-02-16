require(data.table)
source('~/code/single_cell/tpm.r')


fast_gsea = function(ranks, pathways='~/aviv/db/gsea/pathways.rds', minSize=10, maxSize=500, nperm=1000){
    require(fgsea)
    
    # inputs:
    # ranks = list(gene1 = pval, gene2 = pval)
    # pathways = list(db1 = list(sig1 = genes, sig2 = genes, ...), db2 = list(sig1 = genes, sig2 = genes, ...))

    # notes:
    # this function should be run on ranks from full gene universe, e.g.
    # markers = markers[pvalH < .05]
    # ranks = structure(markers$log2fc, names=markers$gene)
    
    # load default pathways
    if(!is.list(pathways)){pathways = readRDS(pathways)}
    
    gsea = list()

    # run gsea on each set of pathways
    for(name in names(pathways)){
        print(paste(name, 'testing', length(pathways[[name]]), 'gene sets'))

	# intersect pathway with gene universe
	pathway = pathways[[name]]
	pathway = sapply(pathway, function(a) intersect(a, names(ranks)))

	# run gsea and sort by pvalue
        res = fgsea(pathways=pathways[[name]], stats=ranks, nperm=nperm, minSize=minSize, maxSize=maxSize)
	res = res[order(res$pval),]
	gsea[[name]] = res
    }
    
    return(gsea)
}


fast_enrich = function(genes, regex='GO_.*2017$|KEGG.*2016|Reactome.*2016|Panther_2016', collapse=FALSE){
    
    require(enrichR)
    require(data.table)
    
    # For full list of databases: listEnrichrDbs()
    # Computes enrichment with Fisher exact test
    # Also uses random gene sets to correct the Fisher test
    
    # Select databases to use
    dbs = grep(regex, listEnrichrDbs()[,1], value=T)
    
    # Run enrichment test
    res = enrichr(genes, dbs)

    # Fix each list
    res = sapply(res, function(a) {

        # Sort by adjusted p-value
        a = as.data.table(a)[,.(Term, Overlap, P.value, Adjusted.P.value, Genes)][order(Adjusted.P.value)]

	# Get overlap statistics
	a[, NM := as.integer(gsub('/.*', '', Overlap))]
	a[, NQ := length(genes)]
	a[, NT := as.integer(gsub('.*/', '', Overlap))]

    }, simplify=F)

    # Return results
    if(collapse == TRUE){res = do.call(rbind, res)[order(P.value)]}
    res
}


gsea_heatmap = function(terms=NULL, names=NULL, genes=NULL, max_names=50, min_genes=2, show_tree=FALSE, fix_names=TRUE, out=NULL, ...){
    require(NMF)
    
    if(!is.null(terms)){
        if('Term' %in% colnames(terms) & 'Genes' %in% colnames(terms)){
            names = terms$Term
	    genes = terms$Genes
	}
	if('pathway' %in% colnames(terms) & 'leadingEdge' %in% colnames(terms)){
	    names = terms$pathway
	    genes = terms$leadingEdge
	}
    }
    
    # Split genes by [,; ]
    if(! is.list(genes)){genes = strsplit(genes, ',|;| ')}
    
    # Filter terms
    i = sapply(genes, length) >= min_genes
    names = names[i]
    genes = genes[i]
    
    if(length(names) > max_names){
        names = names[1:max_names]
	genes = genes[1:max_names]
    }

    # Get all genes
    all_genes = sort(unique(unlist(genes)))

    # Fix names
    if(fix_names == TRUE){
        names = gsub('\\(GO[^\\)]*\\)', '', names)
	names = gsub('positive', 'pos', names)
	names = gsub('negative', 'neg', names)
	names = gsub('regulation', 'reg', names)
	names = gsub('response', 'resp', names)
	names = gsub('signaling', 'sig', names)
	names = gsub('interaction', 'ix', names)
	names = substr(names, 1, 40)
    }
    
    # Incidence matrix
    x = t(sapply(genes, function(gi) all_genes %in% gi))
    rownames(x) = names
    colnames(x) = all_genes
    
    # Make heatmap
    aheatmap(ifelse(x == TRUE, 1, 0), scale='none', hclustfun='complete', Rowv=show_tree, Colv=show_tree, color='Blues:100', border='black', filename=out, ...)
}


pair_enrich = function(pairs, gene_sets, universe='~/aviv/db/map_gene/hg19_genes.txt', require_ix=TRUE){
    require(data.table)
    
    # Calculate enrichment scores for all gene pairs (e.g. ligand-receptor interactions) and gene sets (e.g. GO terms)
    # This function uses the minimum overlap found in each column of the gene pair
    # If require_ix, then both genes in the gene pair must be found in the gene set
    
    # Read input data
    universe = readLines(universe)
    pairs = as.data.table(pairs)
    colnames(pairs) = c('a', 'b')
    
    # Select gene universe
    pairs = pairs[a %in% universe & b %in% universe]
    gene_sets = sapply(gene_sets, function(a) intersect(a, universe), simplify=F)

    # Pre-calculate sizes
    nu = length(unique(universe))
    np = length(unique(unlist(pairs)))
    
    # Fisher test on gene pairs
    res = sapply(names(gene_sets), function(name){        
        gene_set = gene_sets[[name]]

	if(require_ix == TRUE){
	    p = pairs[a %in% gene_set & b %in% gene_set]
	} else {
	    p = pairs
	}
	
	# Calculate intersections
	a = unique(intersect(p$a, gene_set))
	b = unique(intersect(p$b, gene_set))

	# (2x2) contingency table
	u = 2*min(length(a), length(b))
	v = length(gene_set) - u
	w = np - u
	x = nu - u - v - w

	# Fisher test
	if(u < 2){return(NULL)}
	data.frame(Term=name, P.value=fisher.test(matrix(c(u,v,w,x), nrow=2, ncol=2))$p.value, Ligand=paste(a, collapse=','), Receptor=paste(b, collapse=','), N=u)
    })
    
    as.data.table(do.call(rbind, res))[order(P.value)]
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


go_genes = function(genes, universe, ontology='BP'){

    require(topGO)
    GO2genes = readMappings(file='~/aviv/db/gopca/go.test.txt', sep='\t')
    gene_list = as.numeric(universe %in% genes)
    names(gene_list) = universe
    gene_list = factor(gene_list)
    
    # Run topGO tests
    GOdata = new('topGOdata', ontology=ontology, allGenes=gene_list, annot=annFUN.GO2genes, GO2genes=GO2genes)
    res = runTest(GOdata, algorithm='classic', statistic='fisher')
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

