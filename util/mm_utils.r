# Utilities for working with model matrices


mm_reverse = function(mm, form, data){

    # Map a model matrix to its associated columns and factors in data
    # Returns list(1=df(columns, factors), 2=df(columns, factors), ...)

    # Find columns
    ttf = attr(terms(as.formula(form), data=data), 'factors')
    mmi = attr(mm, 'assign')
    ttf = ttf[,mmi,drop=F]
    
    # Find factors
    # for each column of the model matrix...
    mm.cols = grep('Intercept', colnames(mm), invert=T, value=T)
    
    res = lapply(1:ncol(ttf), function(j){

        # retrieve colnames of data
        cj = rownames(ttf)[which(ttf[,j] == 1)]
	
	# find the factors of data
        fj = strsplit(mm.cols[j], ':')[[1]]
        fj = sapply(1:length(cj), function(i) gsub(cj[[i]], '', fj[[i]]))

	# return matrix with 1=columns, 2=factors
	cbind(cj, fj)
    })
    return(res)
}


mm_logical_not = function(mm, form, data, method='multi', invert='last'){

    # Calculate the logical inverse of a model matrix
    # Each column of the model matrix is either a single term A or an interaction A:B
    # If single term, then logical not is ~A
    # If interaction, then logical not is A:~B (invert = last) or ~A:~B (invert = all)
    # For multi-level factors, ~A is the next level of A (method = multi) or !A (method = binary)
    
    # find columns and factors associated with the model matrix
    vars = mm_reverse(mm, form, data)
    
    # for each model matrix column...
    mm.cols = grep('Intercept', colnames(mm), invert=T, value=T)
    res = sapply(1:length(vars), function(j){

	# setup the logical not vector
	u = rep(TRUE, nrow(data))
	n = ''
	
	# iterate through all (column,factor) pairs
	cat(paste0('\nInverse for ', mm.cols[[j]], ': '))
	vj = vars[[j]]
	for(i in 1:nrow(vj)){

	   # get the column, factor, and next highest level
	   ci = vj[i,1]
	   fi = vj[i,2]
	   li = levels(data[,ci])[which(levels(data[,ci]) == fi) - 1]
	   
	   # inversion logic
	   if(i < nrow(vj)){
	       if(invert == 'all'){
	           n = paste0(n, fi, '-')
		   u = u & data[,ci] != fi
	       } else if(invert == 'last'){
		   n = paste0(n, fi, '+ ')
		   u = u & data[,ci] == fi
	       }
	   } else {
	       if(method == 'binary'){
	           n = paste0(n, fi, '- \n')
		   u = u & data[,ci] != fi
	       } else if(method == 'multi'){
	           n = paste0(n, li, '+')
		   u = u & data[,ci] == li
	       } else {
	           stop(paste('Error: invalid method', method))
	       }
	   }
	}
	cat(paste0(n, '\n'))
	
	# return the logical not vector
	q = list()
	q[[n]] = u
	return(q)
    })
    res = as.data.frame(res, check.names=F)
    refs = gsub(' ', '', colnames(res))
    colnames(res) = mm.cols
    return(list(x=res, names=refs))
}
