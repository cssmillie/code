library(Seurat)
library(data.table)
library(foreach)
library(parallel)
library(Rtsne)


make_seurat = function(name, dge=NULL, min_cells=10, min_genes=500, max_genes=4000){

    print('Loading DGE')
    counts = data.frame(fread(paste('zcat', dge)), row.names=1)
    
    print('Filtering DGE')
    i = rowSums(counts > 0) >= min_cells
    j1 = colSums(counts > 0) >= min_genes
    j2 = colSums(counts > 0) <= max_genes
    counts = counts[i, (j1 & j2)]
    
    print('Transforming data')
    data = 10000*scale(counts, center=FALSE, scale=colSums(counts))
    data = data.frame(log2(data + 1))

    print('Making Seurat object')
    seur = new('seurat', raw.data = counts)
    seur = setup(seur, project=name, min.cells=0, min.genes=0, calc.noise=F, is.expr=0, names.delim='\\.', names.field=1)
    seur@data = data
    seur@scale.data = t(scale(t(seur@data)))
    
    return(seur)
}


run_seurat = function(name, dge=NULL, min_cells=10, min_genes=500, max_genes=4000, do.batch=F, perplexity=25, max_iter=1000, cluster='infomap', k=c(), ncores=1){

    seur = make_seurat(name=name, dge=dge, min_cells=min_cells, min_genes=min_genes, max_genes=max_genes)

    print('Selecting variable genes')
    var_genes = select_var_genes(seur@raw.data, diffCV.num_genes=num_genes)
    seur@var.genes = intersect(var_genes, rownames(seur@data))
    
    if(do.batch){
	seur = batch_correct(seur, seur@ident)
    }
    
    print('PCA')
    num_pcs = sig.pcs.perm(t(seur@scale.data[seur@var.genes,]), randomized=T, n.cores=ncores)$r

	seur@data.info$num_pcs = num_pcs
	msg(name, sprintf('Found %d significant PCs', num_pcs), verbose)
	seur = run_rpca(seur, k=25, genes.use=seur@var.genes)

	msg(name, 'TSNE', verbose)
	if(tsne_cor == F){
	    seur = run_tsne(seur, dims.use=1:num_pcs, do.fast=T, max_iter=max_iter, verbose=T, perplexity=perplexity)
        } else {
	    d = as.dist(1 - cor(t(seur@pca.rot[,1:num_pcs])))
	    q = Rtsne(d, is_distance=T, do.fast=T, max_iter=max_iter, verbose=T, perplexity=perplexity)$Y
	    rownames(q) = colnames(seur@data)
	    colnames(q) = c('tSNE_1', 'tSNE_2')
	    seur@tsne.rot = as.data.frame(q$Y)
	}
	
    }
    
    msg(name, 'PC loadings', verbose)
    loaded_genes = get.loaded.genes(seur@pca.obj[[1]], components=1:num_pcs, n_genes=20)

    msg(name, 'Calculate signatures', verbose)
    seur = update_signatures(seur)
    
    if(length(k) > 0){
        
	msg(name, 'Saving backup Seurat object', verbose) 
	if(do.backup){
	    out = name
	    saveRDS(seur, file=paste0(name, '.seur.rds'))
	} else {
	    out = NULL
	}
	
	msg(name, 'Clustering cells', verbose)
	k = k[k < ncol(seur@data)]
	u = paste('Cluster.Infomap.', k, sep='')
	if(do.backup){prefix = name} else {prefix = NULL}
	v = run_cluster(seur@pca.rot[,1:num_pcs], k, method=cluster, weighted=FALSE, n.cores=min(length(k), ncores), dist='cosine', do.fast=T, prefix=prefix)
	seur@data.info[,u] = v
    }

    if(write_out){
	
	write.table(loaded_genes, paste0(name, '.loaded_genes.txt'), sep='\t', quote=F)

	png(paste0(name, '.tsne.png'), width=1200, height=900)
	tsne.plot(seur, pt.size=1, label.cex.text=.25)
	dev.off()
	
	if(length(k) > 0){
	pdf(paste0(name, '.clusters.pdf'), width=9, height=9)
	plot_clusters(seur)
	dev.off()
	}

	msg(name, 'Plot summary statistics', verbose)
	png(paste0(name, '.summary_stats.png'), width=1200, height=800)
	plot_tsne(seur, subset(seur@data.info, select=c(orig.ident, G1S, G2M, nGene, nUMI)), do.label=F)
	dev.off()

	saveRDS(seur, file=paste0(name, '.seur.rds'))
    }

    return(seur)
}


select_var_genes = function(count.data, reads.use=FALSE, do.text=FALSE, diffCV.cutoff=NULL, diffCV.num_genes=NULL, do.spike=FALSE, main.use=NULL, prefix=NULL, do.plot=FALSE, ret.diffCV=FALSE){
    
    # Empirical mean, var and CV
    mean_emp = apply(count.data, 1, mean)
    var_emp = apply(count.data, 1, var)
    cv_emp = sqrt(var_emp) / mean_emp
    
    # NB sampling
    a=colSums(count.data)
    size_factor =  a/ mean(a)
    fit=fitdistr(size_factor, "Gamma")
    if (do.spike) spike.genes=grep("^ERCC", rownames(count.data), value=TRUE)
    print(fit)
    
    if(do.plot==T){
    par(mfrow=c(2,2))
    if (!reads.use){
        hist(size_factor, 50, probability=TRUE, xlab="N_UMI/<N_UMI>", main = main.use)
    } else {
        hist(size_factor, 50, probability=TRUE, xlab="N_Reads/<N_Reads>", main = main.use)
    }
    curve(dgamma(x, shape=fit$estimate[1], rate=fit$estimate[2]),from=0, to=quantile(size_factor, 0.999), add=TRUE, col="red", main="Gamma dist fit for size factor")
    text(5,0.6, paste("shape = ", round(fit$estimate[1],2)))
    text(5,0.5, paste("rate = ", round(fit$estimate[2],2)))
    }
    
    # Gamma distributions of individual genes are just scaled versions. If X ~ Gamma(a,b)
    # then cX ~ Gamma(a, b/c)
    a_i = rep(fit$estimate[1], length(mean_emp)); names(a_i) = names(mean_emp)
    b_i = fit$estimate[2] / mean_emp; names(b_i) = names(mean_emp)
    mean_NB = a_i / b_i; var_NB = a_i*(1+b_i) / (b_i^2)
    cv_NB = sqrt(var_NB)/mean_NB
    
    diffCV = log(cv_emp) - log(cv_NB)
    if(ret.diffCV == TRUE){
        return(diffCV)
    }
    
    if(!is.null(diffCV.cutoff)){
        pass.cutoff=names(diffCV)[which(diffCV > diffCV.cutoff & (mean_emp > 0.005 & mean_emp < 100))]
    } else if (!is.null(diffCV.num_genes)){
        pass.cutoff=names(sort(diffCV,decreasing=T)[1:diffCV.num_genes])
    }
    
    if(do.plot == T){
    plot(mean_emp,cv_emp,pch=16,cex=0.5,col="black",xlab="Mean Counts",ylab="CV (counts)", log="xy", main = main.use)
    if (do.spike) points(mean_emp[spike.genes],cv_emp[spike.genes],pch=16,cex=0.5,col="red")
    curve(sqrt(1/x), add=TRUE, col="red", log="xy", lty=2, lwd=2)
    or = order(mean_NB)
    lines(mean_NB[or], cv_NB[or], col="magenta", lwd=2)
    if(do.text) text(mean_emp[pass.cutoff],cv_emp[pass.cutoff],pass.cutoff,cex=cex.text.use)
    hist(diffCV, 50, probability=TRUE, xlab="diffCV")
    par(mfrow=c(1,1))
    }
    
    return(as.character(pass.cutoff))
}


batch_correct = function(seur, batch){
    library(sva)
    print('Running ComBat with the following batches:')
    print(table(batch))
    new.data = ComBat(dat=seur@data, batch, par.prior=T, prior.plots=F)
    seur@data = as.data.frame(new.data)
    seur@scale.data = t(scale(t(seur@data)))
    return(seur)
}


sig.pcs.perm <- function (dat, B = 100, threshold = 0.05, randomized=F, verbose=TRUE, seed = NULL, max.pc=100, n.cores=1, center=T, scale=T) {
    
    ptm <- proc.time()
    if(B %% n.cores != 0){stop("Permutations must be an integer multiple of n.cores")}
    cat(sprintf("Scaling input matrix [center=%s, scale=%s]\n", center, scale))
    dat = as.matrix(scale(dat, center=center, scale=scale))
    #dat = as.matrix(t(scale(t(dat), center=center, scale=scale)))
    if (!is.null(seed)) set.seed(seed)
    n <- min(max.pc, ncol(dat))
    m <- nrow(dat)
    print(paste0("Considering only the top ", n, " PCs. Supply max.pc if you wish to change"))
    cat(sprintf("Running initial PCA\n"))
    if(randomized){
        library(rsvd)
        uu <- rsvd(as.matrix(dat), k=max.pc)
    }else{
        uu <- corpcor::fast.svd(dat, tol = 0)
    }
    
    ndf <- n - 1
    dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
    dstat0 <- matrix(0, nrow = B, ncol = ndf)
    if(verbose==TRUE) message("Estimating number of significant principal components. Permutation: ")
    
    #permutations
    if(n.cores==1){
        for (i in 1:B) {
            if(verbose==TRUE) cat(paste(i," "))
            dat0 <- t(apply(dat, 1, sample, replace = FALSE))
            if(randomized){
                library(rsvd)
                uu0 <- rsvd(as.matrix(dat0), k=max.pc)
            }else{
                uu0 <- corpcor::fast.svd(dat0, tol = 0)
            }
            dstat0[i, ] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
        }
    }else{
        library(parallel)
        library(foreach)
        library(doParallel)
        cl<-makePSOCKcluster(n.cores, outfile="")
        registerDoParallel(cl, n.cores)
        chunksize = B/n.cores
        vals = split(1:B, ceiling(seq_along(1:B)/chunksize))
        dstat0 = foreach(run.id=1:n.cores, .packages="corpcor", .combine=cbind) %dopar% {
            v = vals[[run.id]]
            #cat(sprintf("Core %s will run perms: %s \n", run.id, paste(v, collapse=",")))
            do.call(rbind, lapply(v, function(i) {
                if(verbose==TRUE) cat(paste(i," "))
                dat0 <- t(apply(dat, 1, sample, replace = FALSE))
                
                if(randomized){
                    library(rsvd)
                    uu0 <- rsvd(as.matrix(dat0), k=max.pc)
                }else{
                    uu0 <- corpcor::fast.svd(dat0, tol = 0)
                }
                uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
            }))
            
        }
        cat("\nUnregistering parallel backend..")
        stopCluster(cl)
        registerDoSEQ()
        cat(" done\n");
    }
    p <- rep(1, n)
    for (i in 1:ndf) {
      p[i] <- mean(dstat0[, i] >= dstat[i])
    }
    for (i in 2:ndf) {
      p[i] <- max(p[(i - 1)], p[i])
    }
    r <- sum(p <= threshold)
    y = proc.time() - ptm
    cat(sprintf("\n\n PC permutation test completed. \n %s PCS significant (p<%s, %s bootstraps)\n Runtime: %s s\n ", r,  threshold, B,signif(y[["elapsed"]], 3)))
    
    return(list(r = r, p = p))
}
