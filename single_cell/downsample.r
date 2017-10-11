

require(rsvd)
require(cccd)
require(expm)

knn_downsample = function(seur = NULL, counts = NULL, ident = NULL, ngene = NULL, pca.rot = NULL,
	       		  k1 = 10, dist.use='cosine', transitions = 2, num_pcs = 10, uniform_rates = FALSE,
	       		  k2 = 50, total_cells = NULL, pct_cells = NULL, cells_per_ident = NULL){
    
    # set variables
    if(!is.null(seur) & is.null(counts)){counts = seur@raw.data}
    if(is.null(ngene)){ngene = structure(colSums(counts > 0), names=colnames(counts))}
    if(is.null(ident)){ident = structure(gsub('\\..*', '', colnames(counts)), names=colnames(counts))}
    if(length(ident) == 1){ident = rep(ident, ncol(counts))}
    ident = as.factor(ident)
    
    cat('\n\n----------\ndownsample\n----------\n')
    cat('\n\nSizes before:\n')
    print(sort(table(ident)))
    
    # total cells
    if(!is.null(pct_cells)){total_cells = pct_cells*ncol(counts)/100.}

    # select number of cells per group
    num_cells = num_keep = sort(table(ident))
    if(!is.null(cells_per_ident)){
        num_keep[num_keep > cells_per_ident] = cells_per_ident
    } else {
        n = sort(table(ident))
	u = c(0, cumsum(n)[1:(length(n)-1)])
	i = (total_cells - u)/seq(length(n), 1, -1) < n
	num_keep[i] = as.integer((total_cells - sum(n[!i]))/sum(i))
    }
    cat('\n\nSizes after:\n')
    print(num_keep)
    
    # cells to keep
    cells = data.frame(rep(NA, ncol(counts)), row.names=colnames(counts))
    rates = data.frame(rep(NA, ncol(counts)), row.names=colnames(counts))
    
    # calculate pca
    if(is.null(pca.rot)){
        cat('Calculating PCA\n')
        var_genes = select_var_genes(data=counts, method='karthik', num_genes=1500)
        data = log2(calc_tpm(data=data) + 1)
        pca.rot = as.data.frame(rpca(t(data[var_genes,]), k=num_pcs, retx=TRUE)$x)
    } else {
        pca.rot = pca.rot[,1:num_pcs]
    }
    
    # downsample each group separately
    ident = droplevels(ident)
    for(group in levels(ident)){cat(paste0('\n\nDownsampling ', group, ':\n'))
        
        # skip small groups
	if(num_cells[[group]] <= num_keep[[group]]){
	    cat(paste('Selected', sum(ident == group), 'cells\n'))
	    rates[colnames(data)[ident == group], 1] = NA
	    cells[colnames(data)[ident == group], 1] = 1
	    next
	}
	
        # subset data
        x = data[,ident == group]
	i = ifelse(apply(x, 1, sd) %in% c(NA, 0), FALSE, TRUE)
	x = x[i,]
	
	if(uniform_rates == FALSE){
	    
	    # calculate k1 nearest neighbors
	    cat(paste('Calculating', k1, 'nearest neighbors\n'))
	    g = nng(pca.rot[colnames(x),], k=k1, method=dist.use, use.fnn=FALSE)
	    
	    # calculate transition matrix
	    cat(paste('Calculating', transitions, 'transitions\n'))
	    G = t(as_adj(g))	    
	    G = scale(G, center=F, scale=colSums(G))
	    G = G %^% transitions
	    
	    # number of neighbors for each cell
	    p = rowSums(G > 0) + 1
	    
	    # solve optimal downsampling rates
	    cat('Solving optimal downsampling rates\n')
	    num_remove = num_cells[[group]] - num_keep[[group]]
	    f = function(a){(sum((exp(a)/p)^(1/(p-1))) - num_remove)^2}
	    a = optim(0, f, method='L-BFGS-B', upper=log(min(p)))$par
	    cat(paste0('lambda = ', signif(exp(a),3), ', f(lambda) = ', signif(f(a),3), '\n'))
	    d = structure((exp(a)/p)^(1/(p-1)), names=colnames(x))
	} else {
	    d = (num_cells[[group]] - num_keep[[group]])/num_cells[[group]]
	    d = structure(rep(d, ncol(x)), names=colnames(x))
	}
	    
	cat('Deletion rates:\n')
	print(quantile(d))
	
	# calculate k2 nearest neighbors
	cat(paste('Local subsampling with', k2, 'nearest neighbors'))
	g = nng(pca.rot[colnames(x),], k=k2, method=dist.use, use.fnn=FALSE)
	
	# downsample cells
	keep = sapply(1:ncol(x), function(i){
	    ci = colnames(x)[[i]]
	    nn = colnames(x)[g[[i]][[1]]]
	    return(sum(ngene[[ci]] >= ngene[nn]) >= d[[i]]*k2)
	})
	cells[colnames(x),1] = keep
	rates[colnames(x),1] = d
    }
    cat(paste('\n\nSelected', sum(cells[,1]), 'cells\n\n\n'))
    
    return(list(cells=cells, rates=rates))
}
