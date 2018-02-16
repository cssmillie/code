source('~/code/util/mtx.r')

# impute data with magic

impute_magic = function(data, num_pcs=20, k=30, sparse=TRUE, do.log=FALSE, out=NULL, return=TRUE){
    
    # data = [genes x cells] matrix
    # num_pcs = number of PCs to use
    # k = number of nearest neighbors

    # trackers
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
	cmd = paste0('python ~/code/single_cell/run_magic.py --txt ', tmp, ' -n ', num_pcs, ' -k ', k, ' --out ', out)
    
    } else {

        # write data
        tmp = gsub('.txt', '', tempfile(pattern='magic.', tmpdir='~/tmp', fileext='.txt'))
	fns = write_mtx(t(data), path=tmp)
	files = c(files, fns$data, fns$rows, fns$cols)
	if(is.null(out)){out = fns$data}
	
	# magic command
	cmd = paste0('python ~/code/single_cell/run_magic.py --mtx ', fns$data, ' --genes ', fns$cols, ' -n ', num_pcs, ' -k ', k, ' --out ', out)
    }

    # run magic
    system(cmd)
    
    # load results
    if(return == TRUE){
        data = ffread(out, sep=',')
        data = t(data)
        colnames(data) = cells
    }
    
    # cleanup files
    for(fn in fns){system(paste0('rm ', fn))}

    # return results
    if(return == TRUE){
        return(data)
    } else {
        return(NULL)
    }
}


