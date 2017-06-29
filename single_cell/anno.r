source('~/code/single_cell/map.r')

load_anno = function(fn='~/Gut_Human/csmillie/anno/tree_annotations.txt', key='ident', value='name'){
    x = load_map(fn, key, value)
}

flat_names = function(names){
    # return list of non-overlapping names
    # e.g. flat_names('Tcell') = c('Tcell', list of non-Tcells)
    n2i = load_anno(key='name', value='ident')
    o = select_keys(n2i, names, invert=T)
    c(names, o[sapply(n2i[o], length) == 1])
}

