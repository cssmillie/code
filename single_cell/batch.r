require(Seurat)

load_scran = function(){
    library(BiocGenerics, pos=length(search()))
    library(S4Vectors, pos=length(search()))
    library(DelayedArray, pos=length(search()))
    library(scran)
}

batch_correct = function(seur, batch, design=NULL, method='combat', genes.use='var', liger.data=NULL, ndim=30){
    
    # note: design matrix does not include batch
    if(!is.null(design)){
	if(is.null(rownames(design))){
            if(nrow(design) == ncol(seur@data)){
	        rownames(design) = colnames(seur@data)
	    } else {
	        stop('Error: check dimensions')
	    }
	}
        # subset design matrix
        design = design[colnames(seur@data),,drop=F]
    }
    
    # Get genes.use
    if(length(genes.use) == 1){
        if(genes.use == 'all'){
            print('Running batch correction on all genes')
	    genes.use = rownames(seur@data)
    	} else if(genes.use == 'var'){
            print('Running batch correction on variable genes')
	    if(length(seur@var.genes) == 0){
	        seur@var.genes = get_var_genes(seur, method='loess', num_genes=1500)
	    }
	    genes.use = seur@var.genes
    	} else {
            stop("Error: can only use 'all' or 'var' genes")
    	}
    }
    
    if(method == 'combat'){
        require('sva')	
	if(is.null(design)){
	    print('Running ComBat')
    	    print(table(batch))
	    new.data = ComBat(dat=as.matrix(seur@data[genes.use,]), batch, par.prior=T, prior.plots=F)
	} else {
	    print(sprintf('Fixing batches while controlling for %d covariates', ncol(design)))
	    model = model.matrix(~ ., data=design)
	    print(dim(model))
	    new.data = ComBat(dat=as.matrix(seur@data[genes.use,]), batch=batch, mod=model, par.prior=T, prior.plots=F)
	}	
    }
    
    if(method == 'cca'){
        require(PMA)
	
	# Split data into batches
	data = lapply(split(as.data.frame(t(as.matrix(seur@data[genes.use,]))), batch), t)

	# Batch correct data with lambda=1
    	u = MultiCCA(data, penalty=1, ws=pout$ws.init)$ws

	# Unsplit data into original format
	new.data = unsplit(lapply(u, as.vector), batch)

	new.data = lapply(1:length(u), function(i) data.frame(u[[i]], row.names=paste(i, 1:nrow(u[[i]]), sep='.')))
	new.data = t(unsplit(new.data, batch))
	
    }
    
    if(method == 'multicca'){
    
	# Check parameters
	if(is.null(ndim)){print('ndim is null. hmmm.... setting to 30?'); ndim = 30}
    	
        # Select batch
	batch = as.character(batch)
	batch[is.na(batch)] = 'Other'
	u = sort(table(batch))
	i = names(which(u <= 50))
	batch[batch %in% i] = 'Other'
	u = sort(table(batch))
	print(u)
	if(any(u <= 50)){
	    i = max(which(cumsum(u) <= 50)) + 1
	    i = names(u)[1:min(i, length(u))]
	    batch[batch %in% i] = 'Other'
	}
	batch = as.factor(batch)
	print(table(batch))
	
	# Run MultiCCA
	new.data = RunMultiCCA2(seur=seur, groups=batch, num.ccs=ndim)
	return(as.data.frame(new.data))
    }

    if(method == 'mnn'){
        load_scran()
	
	# Split data into batches
	data = lapply(split(as.data.frame(t(as.matrix(seur@data[genes.use,]))), batch), t)
	
	# Batch correct data
	u = do.call(mnnCorrect, data)$corrected
	
	# Unsplit data into original format
	new.data = lapply(names(u), function(a) data.frame(t(u[[a]]), row.names=paste(a, 1:ncol(u[[a]]), sep='.')))
	new.data = t(unsplit(new.data, batch))
	
	# Fix row and column names
	rownames(new.data) = genes.use
	colnames(new.data) = colnames(seur@data)
    }

    if(method == 'liger'){
        library(liger)

	# Split data into batches
	sparse_split = function(x, g){
	    g = as.factor(g)
	    sapply(levels(g), function(gi){
	        x[g == gi,,drop=F]
	    }, simplify=F)
	}
	
	# Select batch
	batch = as.character(batch)
	batch[is.na(batch)] = 'Other'
	u = sort(table(batch))
	i = names(which(u <= 50))
	batch[batch %in% i] = 'Other'
	u = sort(table(batch))
	print(u)
	if(any(u <= 50)){
	    i = max(which(cumsum(u) <= 50)) + 1
	    i = names(u)[1:min(i, length(u))]
	    batch[batch %in% i] = 'Other'
	}
	batch = as.factor(batch)
	print(table(batch))
	
	new.data = seur@raw.data
	new.data = sapply(sparse_split(t(new.data), batch), t, simplify=F)
	for(i in 1:length(new.data)){
	    new.data[[i]][,1] = new.data[[i]][,1] + runif(nrow(new.data[[i]]), min=0, max=.01)
	}
	
	new.data = createLiger(new.data)
	new.data = normalize(new.data)
	new.data@var.genes = seur@var.genes
	new.data = scaleNotCenter(new.data)
	
	# Batch correct data
	new.data = optimizeALS(new.data, k=ndim)
	new.data = quantileAlignSNF(new.data)
	new.data = new.data@H.norm[colnames(seur@data),]
    }
    
    return(as.data.frame(new.data))
}


# ---------------
# Seurat MultiCCA
# ---------------

GetCors <- function(mat.list, ws, num.sets){
    cors <- 0
    for(i in 2:num.sets){
        for(j in 1:(i-1)){
            thiscor  <-  cor(mat.list[[i]]%*%ws[[i]], mat.list[[j]]%*%ws[[j]])
            if(is.na(thiscor)) thiscor <- 0
            cors <- cors + thiscor
        }
        return(cors)
    }
}

UpdateW <- function(mat.list, i, num.sets, ws, ws.final){
    tots <- 0
    for(j in (1:num.sets)[-i]){
        diagmat <- (t(ws.final[[i]])%*%t(mat.list[[i]]))%*%(mat.list[[j]]%*%ws.final[[j]])
        diagmat[row(diagmat)!=col(diagmat)] <- 0
        tots <- tots + t(mat.list[[i]])%*%(mat.list[[j]]%*%ws[[j]]) - ws.final[[i]]%*%(diagmat%*%(t(ws.final[[j]])%*%ws[[j]]))
    }
    w <- tots/l2n(tots)
    return(w)
}

GetCrit <- function(mat.list, ws, num.sets){
    crit <- 0
    for(i in 2:num.sets){
        for(j in 1:(i-1)){
            crit <- crit + t(ws[[i]])%*%t(mat.list[[i]])%*%mat.list[[j]]%*%ws[[j]]
        }
    }
    return(crit)
}

l2n <- function(vec){
    a <- sqrt(sum(vec^2))
    if(a==0){
        a <- .05
    }
    return(a)
}

RunMultiCCA2 = function(seur, groups, niter=25, num.ccs=10){
    library(irlba)
    
    num.sets = nlevels(groups)
    
    mat.list = sapply(levels(groups), function(gi){print(paste('Scaling', gi))
        get_data(seur, data.use='scale', genes.use=seur@var.genes, cells.use=colnames(seur@data)[groups == gi])
    }, simplify=F)
    
    ws = list()
    for (i in 1:num.sets){
        ws[[i]] <- irlba(mat.list[[i]], nv = num.ccs)$v[, 1:num.ccs, drop = F]
    }
    
    ws.init <- ws
    ws.final <- list()

    cors <- NULL
    for(i in 1:length(ws)){
        ws.final[[i]] <- matrix(0, nrow=ncol(mat.list[[i]]), ncol=num.ccs)
    }
    for (cc in 1:num.ccs){
        print(paste0("Computing CC ", cc))
        ws <- list()
        for (i in 1:length(ws.init)){
            ws[[i]] <- ws.init[[i]][, cc]
        }
        cur.iter <- 1
        crit.old <- -10
        crit <- -20
        storecrits <- NULL
        while(cur.iter <= niter && abs(crit.old - crit)/abs(crit.old) > 0.001 && crit.old !=0){print(paste('Iteration', cur.iter))
            crit.old <- crit
            crit <- GetCrit(mat.list, ws, num.sets)
            storecrits <- c(storecrits, crit)
            cur.iter <- cur.iter + 1
            for(i in 1:num.sets){
                ws[[i]] <- UpdateW(mat.list, i, num.sets, ws, ws.final)
            }
        } 
	for(i in 1:length(ws)){
	    ws.final[[i]][, cc] <- ws[[i]]
        }
        cors <- c(cors, GetCors(mat.list, ws, num.sets))
    }
    results <- list(ws=ws.final, ws.init=ws.init, num.sets = num.sets, cors=cors)
    
    cca.data <- results$ws[[1]]
    for(i in 2:length(mat.list)){
        cca.data <- rbind(cca.data, results$ws[[i]])
    }
    rownames(cca.data) = unname(unlist(sapply(mat.list, colnames)))
    cca.data <- apply(cca.data, MARGIN = 2, function(x){
    if(sign(x[1]) == -1) {
        x <- x * -1
    }
    return(x)
    })
    cca.data[colnames(seur@data),]
}


