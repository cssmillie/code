source('~/code/single_cell/map_gene.r')


psi_log = function(x, base=2, zero=NULL){
    if(is.null(zero)){zero = .5*min(x[x > 0])}
    log(x + zero, base=base)
}

geom_mean = function(x, base=2, zero=NULL){
    base**sum(psi_log(x, base=base, zero=zero))/length(x)
}

colGmeans = function(x, base=2){
    zero = .5*min(x[x > 0])
    apply(x, 2, geom_mean, base=base, zero=zero)
}

predict_dge_type = function(x, bases.use=c(2, exp(1), 10), tol=1e-4){
    # Returns counts, tpm, or log base

    # Select data
    u = as.numeric(x[,1])
    v = as.numeric(x[,2])

    # Check for integers
    if(abs(sum(u - round(u))) <= tol & abs(sum(v - round(v))) <= tol){
        return('counts')
    }
    
    # Check for constant sum
    if(abs(sum(u) - sum(v)) <= tol){
        return('tpm')
    }
    
    # Check for big numbers
    if(max(u) > 32 | max(v) > 32){
        print('predict_log_base: guessing data type = imputed TPM')
        return('tpm')
    }

    # Check for logTPM
    u = sapply(bases.use, function(base){
        abs(sum(base**u) - sum(base**v))
    })
    i = which.min(u)
    
    # Test logTPM deviance
    if(u[[i]] <= tol){
        return(bases.use[[i]])
    } else {
        stop('predict_log_base: unknown data type')
    }
}


calc_tpm = function(seur=NULL, data=NULL, genes.use=NULL, cells.use=NULL, total=1e4){    
    # Calculate TPM from @raw.data or data argument
    
    # Select counts data
    if(is.null(data)){data = seur@raw.data}
    if(is.null(genes.use)){genes.use = rownames(data)}
    if(is.null(cells.use)){cells.use = colnames(data)}
    
    # Get genes and cells
    genes.use = intersect(genes.use, rownames(data))
    cells.use = intersect(cells.use, colnames(data))
    
    # Predict DGE type
    type = predict_dge_type(data[, sample(1:ncol(data), 2)])
    
    if(type == 'counts'){
        require(wordspace)    
        numi = colSums(data[,cells.use])
	data = scaleMargins(data[genes.use, cells.use, drop=F], cols=total/numi)
    }
    if(type == 'tpm'){
        # Do nothing
    }
    if(is.numeric(type)){
        data = type**data[genes.use, cells.use, drop=F]
	data = data - min(data)
    }
    
    return(data)
}


get_data = function(seur, data.use='tpm', tpm=NULL, genes.use=NULL, cells.use=NULL){
    
    # Retrieve data from a Seurat object
    # data.use can be: counts, tpm, log2, data, any matrix
    # optionally, pre-calculate tpm for speed

    if(!is.character(data.use)){
        data = data.use
    }
    if(data.use == 'counts'){
        data = seur@raw.data
    }
    if(data.use == 'tpm'){
        if(is.null(tpm)){
	    data = calc_tpm(seur, genes.use=genes.use, cells.use=cells.use)
	} else {
            data = tpm
        }
    }
    if(data.use == 'log2'){
        if(predict_dge_type(seur@data, bases.use=c(2)) != 2){
	    stop('Error: seur@data log base != 2')
	}
	data = seur@data
    }
    if(data.use == 'scale'){
        if(!is.null(genes.use)){
	    genes.use = intersect(genes.use, rownames(seur@data))
	} else {
	    genes.use = rownames(seur@data)
	}
	data = t(scale(t(seur@data[genes.use,])))
    }
    
    # Subset genes and cells
    if(!is.null(genes.use)){data = data[intersect(genes.use, rownames(data)),]}
    if(!is.null(cells.use)){data = data[,intersect(cells.use, colnames(data))]}

    return(data)
}


map_names = function(seur=NULL, data=NULL, meta=NULL, names=NULL, regex=NULL, files=NULL, file.cols=NULL, file.regex=NULL, top=NULL, source='auto', target='auto'){
    
    # Map input arguments to genes and feats
    
    # Initialize variables
    names = as.list(names)
    genes = feats = c()
    regex = as.list(regex)
    
    # Get data and metadata
    if(is.null(data)){data = seur@data}
    if(is.null(meta)){meta = seur@data.info}
    
    # Map names
    if(!is.null(names)){
        if(target == 'auto'){target = predict_organism(rownames(data)[1:100])}
	genes = sapply(names, function(a){
	    i = a %in% rownames(data)
	    u = a[i]
	    if(sum(!i) > 0){u = c(u, intersect(map_gene(a[!i], target=target), rownames(data)))}
	    unique(u)
	}, simplify=F)
	feats = sapply(names, function(a) a[a %in% colnames(meta)])
    }
    
    # Map regex
    if(!is.null(regex)){ # wrap in structure? names = regex
        genes = c(genes, sapply(regex, function(a) grep(a, rownames(data), value=T, perl=T), simplify=F))
	feats = c(feats, sapply(regex, function(a) grep(a, colnames(meta), value=T, perl=T), simplify=F))
    }
    
    # Map files
    if(!is.null(files)){
        sig = do.call(c, lapply(files, function(file) load_signature(file, file.regex=file.regex, file.cols=file.cols)))
	genes = c(genes, sig)
    }

    # Filter genes and feats
    genes = genes[sapply(genes, length) > 0]
    feats = feats[sapply(feats, length) > 0]
    if(!is.null(top)){
        genes = sapply(genes, function(a) as.character(na.omit(a[1:top])), simplify=F)
	feats = sapply(feats, function(a) as.character(na.omit(a[1:top])), simplify=F)
    }
    
    # Fix gene names
    genes = genes[sapply(genes, length) > 0]
    if(length(genes) > 0){
        if(is.null(names(genes))){names(genes) = rep('', length(genes))}
	names(genes)[names(genes) == ''] = sapply(genes[names(genes) == ''], paste, collapse='.')
    }
    
    # Fix feat names
    feats = feats[sapply(feats, length) > 0]
    if(length(feats) > 0){
        if(is.null(names(feats))){names(feats) = rep('', length(feats))}
	names(feats)[names(feats) == ''] = sapply(feats[names(feats) == ''], paste, collapse='.')
    }
    
    return(list(genes=genes, feats=feats))
}


score_cells = function(seur=NULL, data=NULL, meta=NULL, names=NULL, regex=NULL, files=NULL, file.cols=NULL, file.regex=NULL, top=NULL, source='auto', target='auto', scores=NULL,
                       data.use='tpm', combine_genes='mean', groups=NULL, group_stat='mean', do.log=FALSE, log_zero='auto', genes_first=TRUE, cells.use=NULL, make.names=TRUE){
    require(Matrix)
    require(Matrix.utils)
    
    # Score gene expression across cells and optionally aggregate
    # The default steps are:
    # 1) Select data (data.use = 'tpm', 'log2', or 'scaled')
    # 2) Calculate mean expression across genes (combine_genes = 'sum', 'mean', or 'gmean')
    # 3) Calculate mean expression within each cell type (group_by, group_stat = 'mean', 'alpha', 'mu')
    # 4) Log transform results (do.log)
    # If genes_first == FALSE, then calculate the group means *before* combining across genes
    
    # Fix input arguments
    if(is.null(data)){data = get_data(seur, data.use=data.use, cells.use=cells.use)}
    if(is.null(meta)){meta = seur@data.info}
    if(!is.null(scores)){old_scores = as.data.frame(scores)} else {old_scores=c()}
    if(!is.null(groups)){names(groups) = colnames(seur@data)}
       
    # Get genes and feats
    res = map_names(seur=seur, data=data, meta=meta, names=names, regex=regex, files=files, file.cols=file.cols, file.regex=file.regex, top=top, source=source, target=target)
    genes = res$genes
    feats = res$feats
        
    # Subset cells
    if(!is.null(cells.use)){
        data = data[,cells.use]
        meta = meta[cells.use,]
        if(!is.null(groups)){groups = groups[cells.use]}
        if(!is.null(scores)){scores = scores[cells.use,,drop=F]}
    }
    
    group_genes = function(x, method){
        
        # combine expression data across genes within a signature
	# x = [genes x cells] matrix
	# method = 'sum', 'mean', or 'gmean'
	# returns [cells]-vector
	
	if(nrow(x) == 1){return(x[1,])}
	
	if(method == 'sum'){
	    colSums(x)
	} else if(method == 'mean'){
	    colMeans(x)
	} else if(method == 'gmean'){
	    colGmeans(x)
	} else {
	    stop('Error: invalid combine_genes method')
	}
    }
    
    group_cells = function(x, groups, method){
        
        # combine expression data across cells
	# x = [genes x cells] matrix
	# group_stat = 'alpha', 'mu', or 'mean'
	# returns [genes x groups] matrix
	
	if(is.null(groups)){return(x)}
	if(method == 'alpha'){x = x > 0}
	if(method == 'mu'){x[x == 0] = NA}
	x = t(data.frame(aggregate(t(x), list(groups), mean, na.rm=T), row.names=1))
	x[is.na(x)] = 0
	x
    }
    
    # Calculate scores
    names.use = unique(c(names(genes), names(feats)))
    scores = sapply(names.use, function(name){
        
        # Combine data and metadata
	if(name %in% names(genes)){si = data[genes[[name]],,drop=F]} else {si = c()}
	if(name %in% names(feats)){if(is.null(si)){si = t(meta[,feats[[name]],drop=F])} else {si = rBind(si, t(meta[,feats[[name]],drop=F]))}}
	si = as.matrix(si)
	
	if(genes_first == TRUE){
	    si = group_genes(si, method=combine_genes)
	    si = group_cells(t(si), groups=groups, method=group_stat)[1,]
	} else {
	    si = group_cells(si, groups=groups, method=group_stat)
	    si = group_genes(si, method=combine_genes)
	}
	si
    })

    # Fix names
    colnames(scores) = names.use
    if(make.names == TRUE){colnames(scores) = make.names(colnames(scores))}
    
    # Log transform
    if(do.log == TRUE){scores = psi_log(scores, zero=.5*min(data@x))}
    return(as.data.frame(cbind(scores, old_scores)))
}

