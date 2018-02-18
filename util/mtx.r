

sparse_cbind = function(M){
    require(Matrix)
    
    # Initialize matrix
    rows = unique(sort(unlist(sapply(M, rownames))))
    x = Matrix(0, nrow=length(rows), sparse=T)
    rownames(x) = rows

    # Fast merge with rbind/cbind
    for(name in names(M)){

        # Align data
        m = M[[name]]
	n = Matrix(0, nrow=length(rows)-nrow(m), ncol=ncol(m), sparse=T)
	rownames(n) = setdiff(rows, rownames(m))
	m = rbind(m, n)
	m = m[rows,,drop=F]

	# Combine data
	x = cbind(x, m)
    }

    # Remove first column
    x = x[,2:ncol(x),drop=F]
    x
}


sparse_rbind = function(M){
    require(Matrix)
    
    # Initialize matrix
    cols = unique(sort(unlist(sapply(M, colnames))))
    x = Matrix(0, ncol=length(cols), sparse=T)
    colnames(x) = cols

    # Fast merge with rbind/cbind
    for(name in names(M)){

        # Align data
	m = M[[name]]
	n = Matrix(0, ncol=length(cols)-ncol(m), nrow=nrow(m), sparse=T)
	colnames(n) = setdiff(cols, colnames(m))
	m = cbind(m, n)
	m = m[,cols,drop=F]

	# Combine data
	x = rbind(x, m)
    }

    # Remove first row
    x = x[2:nrow(x),,drop=F]
    x
}


fast_rnorm = function(x){
    require(wordspace)
    scaleMargins(x, rows=1/rowSums(x))
}


fast_cnorm = function(x){
    require(wordspace)
    scaleMargins(x, cols=1/colSums(x))
}



mtx_filenames = function(prefix, data='matrix.mtx', rows='genes.tsv', cols='barcodes.tsv'){

    # Get sparse matrix filenames from prefix and patterns
    
    # Get prefix
    if(dir.exists(prefix)){prefix = paste0(prefix, '/')}
    if(!file.exists(prefix)){prefix = paste0(prefix, '.')}
    
    # Fix names
    prefix = gsub('//', '/', prefix)
    prefix = gsub('\\.\\.', '.', prefix)
    prefix = gsub('/\\.', '/', prefix)
    
    # Get filenames
    if(!file.exists(data)){data = paste0(prefix, data)}
    if(!file.exists(rows)){rows = paste0(prefix, rows)}
    if(!file.exists(cols)){cols = paste0(prefix, cols)}
    
    return(list(data=data, rows=rows, cols=cols))
}

read_mtx = function(prefix, data='matrix.mtx', rows='genes.tsv', cols='barcodes.tsv', fix_duplicates=FALSE){

    # Get filenames
    fns = mtx_filenames(prefix=prefix, data=data, rows=rows, cols=cols)
    data = fns$data
    rows = fns$rows
    cols = fns$cols
    
    # Read data
    data = readMM(data)
    rows = read.table(rows)
    rows = rows[,ncol(rows)]
    cols = read.table(cols)
    cols = cols[,ncol(cols)]
    
    # Set names
    if(fix_duplicates == TRUE){data = as(rowsum(as.matrix(data), rows), 'sparseMatrix')} else {rownames(data) = rows}
    colnames(data) = cols
    
    return(data)
}

write_mtx = function(x, prefix='.', data='matrix.mtx', rows='genes.tsv', cols='barcodes.tsv', temp=FALSE){

    # Get filenames
    if(temp == TRUE){prefix = tempfile(tmpdir='~/tmp')}
    fns = mtx_filenames(prefix=prefix, data=data, rows=rows, cols=cols)
    data = fns$data
    rows = fns$rows
    cols = fns$cols
    
    # Write data
    writeMM(x, data)
    writeLines(rownames(x), rows)
    writeLines(colnames(x), cols)

    # Return filenames
    list(data=data, rows=rows, cols=cols)
}
