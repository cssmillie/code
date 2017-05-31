library(Seurat)
   
batch_correct = function(seur, batch, design=NULL, method='combat'){
    
	if(method == 'combat'){library('sva');
		
		if(is.null(design)){
			print('Running ComBat')
    		print(table(batch))
			new.data = ComBat(dat=seur@data, batch, par.prior=T, prior.plots=F)
		} else {
			print(sprintf('Fixing batches while controlling for %d covariates', ncol(design)-1))
			model = model.matrix(~ ., data=design)
			new.data = ComBat(dat=seur@data, batch, par.prior=T, prior.plots=F, mod=model)
		}
		
		seur@data = as.data.frame(new.data)
		seur@scale.data = t(scale(t(seur@data)))
	}
	
	return(seur)
}
