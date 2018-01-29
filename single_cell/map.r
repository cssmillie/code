load_map = function(fn, key, value){

    # Construct a map between two comma-separated columns of a mapping file

    # Load mapping file
    x = read.table(fn, stringsAsFactors=F, header=T)

    # Map names to indices
    n2i = structure(1:ncol(x), names=colnames(x))
    
    # Construct map
    x = unique(do.call(rbind, apply(x, 1, function(xi){
        k = unlist(strsplit(xi[[n2i[[key]]]], ','))
	v = unlist(strsplit(xi[[n2i[[value]]]], ','))
	expand.grid(k,v)
    })))
    colnames(x) = c('k', 'v')

    # Aggregate by key
    x = data.frame(aggregate(v ~ k, x, function(a){c(as.character(a))}), row.names=1, stringsAsFactors=F)
    x = structure(x[,1], names=rownames(x))
    
    # Return a named list of keys -> values
    return(x)
}


select_keys = function(map, keys, invert=FALSE){
    
    # Get indices of matching (or mismatching) values
    i = sapply(map, function(v){
        length(intersect(v, unlist(map[keys]))) > 0
    })
    if(invert == TRUE){i = !i}

    # Return matching (or mismatching) keys
    return(names(map)[i])
}


flatten_keys = function(map, keys){

    # Given a map and a list of keys, this function returns:
    # keys + all non-overlapping "base" keys (length = 1)
    # Example:
    # > anno = load_anno(key='name', value='ident')
    # > good = c('Tcell', 'B', 'N', 'F_Fib', 'M')
    # > flatten_keys(anno, good)
    # > c('B', 'N', ..., 'E_Enterocytes', 'E_Tuft', ..., 'F_Glial')

    unmapped = keys[! keys %in% names(map)]
    if(length(unmapped) > 0){stop(unmapped)}
    others = select_keys(map, keys, invert=TRUE)
    others = others[sapply(map[others], length) == 1]
    return(c(keys, others))
}


load_anno = function(fn='~/Gut_Human/current/anno/all.anno.txt', key='ident', value='name', fix=c('description')){
    
    # Load annotation map
    x = load_map(fn, key, value)
    
    # Reformat text
    fix_text = function(a){
        a = gsub('\\\\n', '\n', a)
	a = gsub('_', ' ', a)
	a
    }
    
    if(key %in% fix){
        names(x) = fix_text(names(x))
    }
    
    if(value %in% fix){
        x = sapply(x, fix_text)
    }
    
    return(x)
}
    
