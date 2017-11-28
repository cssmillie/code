require(data.table)
source('~/code/single_cell/tpm.r')


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

