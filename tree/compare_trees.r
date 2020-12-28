source('~/code/tree/treedist.r')

# compare two large phylogenetic trees
compare_trees = function(g1, g2, org='Escherichia_coli', elen=T){
    # load trees
    t1 = read.tree(paste0('~/Gut_Human/hyper/', org, '/tree/RAxML_bestTree.', g1))
    t2 = read.tree(paste0('~/Gut_Human/hyper/', org, '/tree/RAxML_bestTree.', g2))
    # align trees
    res = align_trees(t1, t2, label_fxn=function(a) gsub('\\..*', '', a), ntip=1e5)
    print(res)
    t1 = res[[1]]
    t2 = res[[2]]
    
    # midpoint root
    t1 = midpoint.root(t1)
    t2 = midpoint.root(t2)
    print(tdist_pair(t1, t2, label_fxn=function(a) gsub('\\..*', '', a), type='quartet', normalize=T, ntip=450))
    print(tdist_pair(t1, t2, label_fxn=function(a) gsub('\\..*', '', a), type='wRF', normalize=T, ntip=450))
    
    # color map
    cmap = setNames(material.heat(length(t1$tip.label)), t1$tip.label)
    c1 = cmap[t1$tip.label]
    c2 = cmap[t2$tip.label]
    t1$tip.label[] = 'x'
    t2$tip.label[] = 'x'
    par(mfrow=c(1,2))
    plot(t1, tip.color=c1, use.edge.length=elen)
    plot(t2, tip.color=c2, use.edge.length=elen)
}
