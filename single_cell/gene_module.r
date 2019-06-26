get_module = function(seur, gene){
    u = as.numeric(seur@data[gene,])
    v = t(as.matrix(seur@data))
    sort(cor(u, v)[1,])
}

gene_modules = function(seur, ident, data, gene, num_pcs=10, method='pearson', scale=FALSE){
    
    require(rsvd)
    require(ppcor)
    
    # subsample data
    cells.use = intersect(rownames(data), colnames(seur@data)[seur@ident == ident])
    ngene = scale(seur@data.info[cells.use, 'nGene'])
    
    # memory efficient scaling
    data = data[cells.use,]
    if(scale == TRUE){
        g = split(1:ncol(data), cut(1:ncol(data), 10))
	for(j in g){
	    data[,j] = scale(data[,j])
	}
    }
    data = as.data.frame(data)
    
    # type 1 = correlation
    mod1 = sort(cor(data[,gene], data, method=method)[1,])
    
    # type 2 = partial correlation
    mod2 = sort(sapply(data, function(a) tryCatch(pcor.test(data[,gene], a, ngene, method=method)$estimate, error=function(e) NA)))
        
    # calculate pca
    pca.obj = rpca(data, center=FALSE, scale=FALSE, retx=TRUE, k=20)
    
    # type 3 = gene loadings
    mod3 = t(pca.obj$rotation)
    mod3 = sort(cor(mod3[,gene], mod3, method=method)[1,])
    
    # set up linear model
    x = as.data.frame(pca.obj$x[,1:num_pcs])
    colnames(x) = paste0('PC', 1:ncol(x))
    x$nGene = ngene
    x = as.matrix(x)
    
    # type 4 = correlation of residuals
    r = sapply(data, function(a) .lm.fit(x, a)$residuals)
    mod4 = sort(cor(r[,gene], r, method=method)[1,])
    
    # type 5 = linear models
    mod5 = sort(sapply(data, function(a) .lm.fit(cbind(a, x), data[,gene])$coefficients[[1]]))

    list(mod1, mod2, mod3, mod4, mod5)
}


mi_modules = function(data, genes.use, cells.use=NULL, nbins=8){
    library(data.table)
    library(entropy)
    
    # fix data
    data = as.data.frame(data)    
    if(!is.null(cells.use)){
        cells.use = intersect(cells.use, rownames(data))
        data = data[cells.use,]
    }
    if(nrow(data) == 0){
        stop('error: cell names')
    }
    if(sum(is.na(data)) > 0){
        stop('error: missing values')
    }
        
    # discretize
    data = data[,sapply(data, sd) > 0]
    data = lapply(data, function(a) cut(a, 8, include.lowest=T))

    # pairwise mutual information
    sapply(genes.use, function(a) sapply(data, function(b) tryCatch({mi.empirical(table(data[[a]], b))}, error=function(e){NA})))
}
