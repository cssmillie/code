require(Seurat)

load_scran = function(){
    library(BiocGenerics, pos=length(search()))
    library(S4Vectors, pos=length(search()))
    library(DelayedArray, pos=length(search()))
    library(scran)
}

batch_correct = function(seur, batch, design=NULL, method='combat', genes.use='var'){

    # note: design matrix does not include batch
    
    # Get genes.use
    if(genes.use == 'all'){
        print('Running batch correction on all genes')
	genes.use = rownames(seur@data)
    } else if(genes.use == 'var'){
        print('Running batch correction on variable genes')
	if(length(seur@var.genes) == 0){
	    seur@var.genes = select_var_genes(seur, method='karthik', num_genes=1500)
	}
	genes.use = seur@var.genes
    } else {
        stop("Error: can only use 'all' or 'var' genes")
    }
    
    if(method == 'combat'){
        require('sva')	
	if(is.null(design)){
	    print('Running ComBat')
    	    print(table(batch))
	    new.data = ComBat(dat=as.matrix(seur@data), batch, par.prior=T, prior.plots=F)
	} else {
	    print(sprintf('Fixing batches while controlling for %d covariates', ncol(design)-1))
	    model = model.matrix(~ ., data=design)
	    print(dim(model))
	    new.data = ComBat(dat=as.matrix(seur@data), batch=batch, mod=model, par.prior=T, prior.plots=F)
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
    
    return(as.data.frame(new.data))
}
