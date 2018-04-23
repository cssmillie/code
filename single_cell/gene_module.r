gene_modules2 = function(seur, ident, data, gene, num_pcs=10, method='pearson', scale=FALSE){
    
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

gene_modules = function(data, gene, ngene, num_pcs=10, method='pearson', scale=FALSE){
    
    require(rsvd)
    require(ppcor)
    
    # fix input arguments
    data = as.data.frame(data)
    ngene = scale(ngene)
    
    # type 1 = correlation
    mod1 = sort(cor(data[,gene], data, method=method)[1,])
    
    # type 2 = partial correlation
    mod2 = sort(sapply(data, function(a) tryCatch(pcor.test(data[,gene], a, ngene, method=method)$estimate, error=function(e) NA)))
    
    # calculate pca
    pca.obj = rpca(data, center=FALSE, scale=FALSE, retx=TRUE, k=20)
    
    # type 3 = gene loadings
    mod3 = t(pca.obj$rotation)
    mod3 = sort(cor(mod3[,gene], mod3, method=method)[1,])
    
    # type 4 = linear model
    x = as.data.frame(pca.obj$x[,1:num_pcs])
    colnames(x) = paste0('PC', 1:ncol(x))
    x$nGene = ngene
    x = as.matrix(x)
    r = sapply(data, function(a) .lm.fit(x, a)$residuals)
    mod4 = sort(cor(r[,gene], r, method=method)[1,])

    list(mod1=mod1, mod2=mod2, mod3=mod3, mod4=mod4)
}
