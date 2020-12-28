nice_rsub = function(f, n, prefix, ...){

    # save current environment to disk
    # write r script that:
    # - loads r environment    
    # - input = integer in 1:n
    # - calls fxn(i)
    # - saves output as .rds
    # - merges output

    # fix input arguments
    if(!grepl('\\.$', prefix)){prefix = paste0(prefix, '.')}
    
    # save environment to disk
    env = tempfile(pattern=prefix, tmpdir='~/tmp', fileext='.RData')
    print(paste('Saving R environment to', env))
    save.image(env)
    
    # write r script
    rscript = tempfile(pattern=prefix, tmpdir='~/tmp', fileext='.r')
    print(paste('Writing parallel Rscript to', rscript))
    res = tempfile(pattern=prefix, tmpdir='~/tmp', fileext='.res')
    print(paste('Results prefix:', res))
    fxn = as.list(match.call())[['f']]
    txt = paste0(' \
        print("loading image") \
	load("', env, '") \    
        i = commandArgs(trailingOnly=T)[[1]] \
	print(i) \
	res = ', fxn, '(i) \
	print(res) \
	saveRDS(res, file=paste0("', res, '", i)) \
    ')
    writeLines(txt, rscript)
    
    # run parallel code
    commands = sapply(1:n, function(i) paste('Rscript', rscript, i))
    rsub(commands, ...)
    
    # load results
    out = sapply(1:n, function(i) readRDS(paste0(res, i)), simplify=F)
    out

}


rsub = function(commands, ...){
    
    # submit unix-style commands to cluster and wait for them to finish
    # example:
    # rsub(c('python a.py', 'python b.py'), m=64, u='csmillie, t='12:00:00')
    
    # write commands to file
    temp = tempfile(tmpdir='~/tmp', fileext='.txt')
    writeLines(commands, temp)
    args = list(...)
    
    # make ssub command
    command = paste('cat', temp, '| ~/code/sge/ssub3 -o rsub -w 5')
    for(k in names(args)){
        v = args[[k]]
	if(length(k) == 1){
	    command = paste0(command, ' -', k)
	} else{
	    command = paste0(command, ' --', k)
	}
	if(!is.null(v)){
	    if(is.numeric(v)){
	        command = paste(command, v)
	    } else {
	        command = paste0(command, " '", v, "'")
	    }
	}
    }
    
    # submit with system
    system(command)
}
