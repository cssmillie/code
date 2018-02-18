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
