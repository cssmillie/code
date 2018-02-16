source('~/code/util/mtx.r')


impute_magic = function(data, num_pcs=20, k=30, t=6, ka=10, eps=1, rescale=99, sparse=TRUE, do.log=FALSE, out=NULL, return=TRUE){
    require(data.table)

    # Run magic imputation on TPM or log2TPM (authors use TPM)
    # see: https://github.com/KrishnaswamyLab/magic
    # data = [genes x cells] matrix
    # k = #nn, t = #transitions, ka = autotune, eps = epsilon, rescale = %rescale
    
    # keep track of cells and filenames
    cells = colnames(data)
    files = c()
    
    # write data
    if(sparse == FALSE){
        
        # write data
        tmp = tempfile(pattern='magic.', tmpdir='~/tmp', fileext='.txt')
	if(is.null(out)){out = tmp}
	files = c(files, tmp)
	fwrite(as.data.frame(t(as.matrix(data))), file=tmp, sep=',', quote=F)
	
	# magic command
	cmd = paste('python ~/code/single_cell/run_magic.py --txt', tmp)
    
    } else {

        # write data
        tmp = gsub('.txt', '', tempfile(pattern='magic.', tmpdir='~/tmp', fileext='.txt'))
	fns = write_mtx(t(data), prefix=tmp)
	files = c(files, fns$data, fns$rows, fns$cols)
	if(is.null(out)){out = fns$data}
	
	# magic command
	cmd = paste('python ~/code/single_cell/run_magic.py --mtx', fns$data, '--genes', fns$cols)
    }
    
    # run magic
    cmd = paste(cmd, '-p', num_pcs, '-k', k, '-t', t, '--ka', ka, '-e', eps, '-r', rescale, '--out', out)
    system(cmd)
    
    # load results
    if(return == TRUE){
        data = ffread(out, sep=',')
        data = t(data)
        colnames(data) = cells
    } else {
        data = NULL
    }
    
    # cleanup files
    for(fn in fns){system(paste0('rm ', fn))}
    
    # return results
    return(data)
}


