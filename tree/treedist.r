library(ape)
library(phangorn)
library(phytools)

remove_duplicate_tips = function(tree, groups=NULL, group_fxn=NULL){
    
    # groups for tip labels
    if(is.null(groups) & is.null(group_fxn)){stop('error: must specify groups')}
    if(any(duplicated(tree$tip.labels))){stop('error: tree has duplicate labels')}
    if(is.null(groups)){
       if(!is.null(group_fxn)){
           groups = group_fxn(tree$tip.label)
       } else {
           groups = tree$tip.label
       }
    }
        
    # select unique tips from each group
    #tips = tapply(tree$tip.label, groups, sample, 1)
    tips = tree$tip.label[!duplicated(groups)]
    
    # drop duplicate tips
    tree = drop.tip(tree, setdiff(tree$tip.label, tips))

    return(tree)
}

rename_tips = function(tree, tips=NULL, rename_fxn=NULL){
    if(!is.null(tips)){
        tree$tip.label = tips
    } else if(!is.null(rename_fxn)) {
        tree$tip = rename_fxn(tree$tip)
        tree$tip.label = rename_fxn(tree$tip.label)
    } else {
        # pass
    }
    return(tree)
}

intersect_tips = function(tree1, tree2, retn=FALSE){
    tips = intersect(tree1$tip.label, tree2$tip.label)
    if(retn == TRUE){return(length(tips))}
    if(length(tips) > 3){
        tree1 = drop.tip(tree1, setdiff(tree1$tip.label, tips))
        tree2 = drop.tip(tree2, setdiff(tree2$tip.label, tips))
        return(list(tree1, tree2))
    } else {
        return(list(NULL, NULL))
    }
}

align_trees = function(tree1, tree2, label_fxn=NULL, ntip=NULL, retn=FALSE){

    # remove duplicate tips
    tree1 = remove_duplicate_tips(tree1, group_fxn=label_fxn)
    tree2 = remove_duplicate_tips(tree2, group_fxn=label_fxn)
    
    # rename tip labels
    tree1 = rename_tips(tree1, rename_fxn=label_fxn)
    tree2 = rename_tips(tree2, rename_fxn=label_fxn)
 
    # intersect tips
    if(retn == TRUE){return(intersect_tips(tree1, tree2, retn=TRUE))}
    trees = intersect_tips(tree1, tree2)
    tree1 = trees[[1]]
    tree2 = trees[[2]]
        
    # subsample tips
    if(!is.null(tree1) & ntip > 0){
        tips = intersect(tree1$tip.label, tree2$tip.label)
	if(length(tips) > ntip){
            i = sample(tips, ntip)
            tree1 = drop.tip(tree1, setdiff(tree1$tip.label, i))
            tree2 = drop.tip(tree2, setdiff(tree2$tip.label, i))
	}
    }
    
    return(c(tree1, tree2))
}

tdist_pair = function(tree1, tree2, label_fxn=NULL, type='RF', normalize=FALSE, ntip=0, ret_ntip=FALSE){

    if(type == 'quartet'){
        if(ntip == 0 | ntip > 450){
	    ntip = 450
	    print('warning: quartet method requires ntip < 477. setting to 450')
	}
    }
    
    trees = align_trees(tree1, tree2, label_fxn=label_fxn, ntip=ntip)
    
    if(is.null(trees[[1]])){
        if(ret_ntip == FALSE){return(NA)} else {return(list(dist=NA, ntip=0))}
    } else {
        ntip = length(trees[[1]]$tip.label)
    }
    
    # tree similarity metrics
    tryCatch({
    if(type == 'RF'){
        di = RF.dist(trees[[1]], trees[[2]], normalize=normalize)
    } else if(type == 'wRF'){
        di = wRF.dist(trees[[1]], trees[[2]], normalize=normalize)
    } else if(type == 'KF'){
        di = KF.dist(trees[[1]], trees[[2]])
    } else if(type == 'path'){
        di = path.dist(trees[[1]], trees[[2]])
    } else if(type == 'wpath'){
        di = path.dist(trees[[1]], trees[[2]], use.weight=TRUE)
    } else if(type == 'spr'){
        di = sprdist(trees[[1]], trees[[2]])
    } else if(type == 'quartet'){
        library(Quartet)
	if(normalize == FALSE){
	    di = TQDist(trees)[1,2]
	} else {
	    di = TQDist(trees)[1,2]/choose(ntip, 4)
	}
    } else if(type == 'nye'){
        library(TreeDist)
	di = NyeTreeSimilarity(trees[[1]], trees[[2]], similarity=FALSE, normalize=T)
    } else if(type == 'mutual'){
        library(TreeDist)
	di = ClusteringInfoDistance(trees[[1]], trees[[2]], normalize=T)
    } else if(type == 'match'){
        library(TreeDist)
	di = MatchingSplitInfoDistance(trees[[1]], trees[[2]], normalize=T)
    } else {
        print('invalid type')
        NULL
    }
    }, except=function(e){
        di = NA
    })

    # return dist and ntip
    if(ret_ntip == FALSE){
        return(di)
    } else {
        return(list(dist=di, ntip=ntip))
    }
}

tdist = function(trees, label_fxn=NULL, type='RF', normalize=FALSE, ntip=0, ret_ntip=FALSE){

    # initialize data
    D = matrix(NA, nrow=length(trees), ncol=length(trees))
    rownames(D) = colnames(D) = names(trees)
    diag(D) = 0
    N = D

    if(type == 'quartet'){
        if(ntip == 0 | ntip > 450){
	    ntip = 450
	    print('warning: quartet method requires ntip < 477. setting to 450')
	}
    }

    for(i in 1:(length(trees)-1)){
        for(j in (i+1):length(trees)){
	    res = tdist_pair(trees[[i]], trees[[j]], label_fxn=label_fxn, type=type, normalize=normalize, ntip=ntip, ret_ntip=TRUE)
	    D[i,j] = D[j,i] = res[[1]]
	    N[i,j] = N[j,i] = res[[2]]
	}
    }
    if(ret_ntip == FALSE){
        return(D)
    } else {
        return(list(dist=D, ntip=N))
    }
}



