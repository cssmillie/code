paramSweep_v3 <- function(seu, PCs=1:10, num_pcs) {
  seurat(); require(fields);
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)
  #pN = c(0.10, 0.20, 0.30)
  
  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@data.info)/(1-0.05) - nrow(seu@data.info))
  pK.test <- round(pK*min.cells)
  pK <- pK[which(pK.test >= 1)]
  
  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (nrow(seu@data.info) > 10000) {
    real.cells <- rownames(seu@data.info)[sample(1:nrow(seu@data.info), 10000, replace=FALSE)]
    data <- seu@raw.data[,real.cells]
    n.real.cells <- ncol(data)
  }
  
  if (nrow(seu@data.info) <= 10000){
    real.cells <- rownames(seu@data.info)
    data <- seu@raw.data
    n.real.cells <- ncol(data)
  }
  
  ## Iterate through pN, computing pANN vectors at varying pK
  output2 <- lapply(as.list(1:length(pN)),
                    FUN = parallel_paramSweep_v3,
                    n.real.cells,
                    real.cells,
                    pK,
                    pN,
                    data,
                    PCs,
		    num_pcs,
		    seu@var.genes)
  
  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for(i in 1:length(output2)){
    for(j in 1:length(output2[[i]])){
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }

  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)

}
