library(MASS)
source('~/dev/adam/rna_seq/r/samples.r')
source('~/dev/adam/rna_seq/r/util.r')


select_var_genes = function(data, method='karthik', vcut=NULL, num_genes=NULL, min.cv2=.25, ret.diffCV=FALSE, do.plot=F){
        
    if(method == 'adam'){
        
	# Get TPM
	data = 10000*scale(data, center=FALSE, scale=colSums(data))
	
	# Select variable genes
	var_genes = get.variable.genes(data, min.cv2=min.cv2)
	i = (!is.na(var_genes[,4])) & (var_genes[,4] <= .05)
	var_genes = var_genes[i, 1]

    }

    if(method == 'karthik'){

	# Get variable genes
	var_genes = meanCVfit(data, diffCV.cutoff=vcut, diffCV.num_genes=num_genes, ret.diffCV=ret.diffCV, do.plot=do.plot)
    
    }
    
    # Return variable genes
    return(var_genes)
}


meanCVfit = function(count.data, reads.use=FALSE, do.text=FALSE, diffCV.cutoff=NULL, diffCV.num_genes=NULL, do.spike=FALSE, main.use=NULL, prefix=NULL, do.plot=FALSE, ret.diffCV=FALSE){

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
