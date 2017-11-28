
merge_clusters = function(seur, min_auc=.65, k=10){
    
    finished = c()
    iter = 0
    pca.rot = seur@pca.rot
    batch = seur@data.info$orig.ident
    
    while(TRUE){iter = iter + 1; print(paste('Iteration', iter))

        ident = as.character(seur@ident)
        print(table(ident))
	
	# distance matrix
	data = data.frame(aggregate(seur@pca.rot, list(ident), mean), row.names=1)
	dist = as.matrix(dist(data))
	
	# batch effects
	batch_remove = tapply(batch, ident, function(a){b=sort(table(a), decreasing=T); sum(b[1:3])/sum(b) >= .85})
	print(paste(sum(batch_remove), 'batch idents'))
	
	# nearest neighbors
    	nn = apply(dist, 1, function(di) colnames(dist)[order(di)[2]])
	nn = t(apply(cbind(names(nn), nn), 1, sort))
	
	# order by distance
	min_dist = apply(dist, 1, function(di) min(di[di > 0]))
	nn = unique(nn[order(min_dist),])
	print(nn)

	# test finished
	uids = apply(nn, 1, paste, collapse=' ')
	if(sum(!uids %in% finished) == 0){break}
	
	# iterate through pairs
	for(i in 1:nrow(nn)){
	
	    # cluster names
	    ci = nn[i,1]
	    cj = nn[i,2]
	    cn = paste(ci, cj, sep='_')
	    print(paste('Testing', ci, cj))
	    
	    # check finished
	    uid = paste(ci, cj)
	    if(uid %in% finished){print(paste('Skipping', ci, cj)); next}
	    finished = c(finished, uid)
	    
	    # test batch
	    if(batch_remove[[ci]] == TRUE | batch_remove[[cj]] == TRUE){
	        print(paste('Batch merge', ci, cj))
	        ident[ident %in% c(ci, cj)] = cn
		seur = set.ident(seur, ident.use=ident)
		break
	    }
	    
	    # test markers
	    m = p.find_markers(seur, ci, cj, test.use = 'roc', dir='both', min_fc=2)
	    auc_hi = sort(na.omit(m$auc), decreasing=T)[k]
	    auc_lo = sort(na.omit(m$auc), decreasing=F)[k]
	    print(paste('AUC', auc_hi, auc_lo))
	    print(paste('CUT', min_auc, 1-min_auc))
	    if(auc_hi < min_auc & auc_lo > (1 - min_auc)){
	        print(paste('Marker merge', ci, cj))
	        ident[ident %in% c(ci, cj)] = cn
		seur = set.ident(seur, ident.use=ident)
		break
	    }
	}
    }
    return(seur)
}
