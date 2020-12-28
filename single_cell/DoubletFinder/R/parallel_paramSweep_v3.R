parallel_paramSweep_v3 <- function(n, n.real.cells, real.cells, pK, pN, data, PCs, num_pcs, var_genes)  {

  sweep.res.list = list()
  list.ind = 0
  
  ## Make merged real-artifical data
  print(paste("Creating artificial doublets for pN = ", pN[n]*100,"%",sep=""))
  n_doublets <- round(n.real.cells/(1 - pN[n]) - n.real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[, real.cells1] + data[, real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data, doublets)
  
  ## Pre-process Seurat object
  seu_wdoublets = run_seurat(name='test', write_out=F, dge=data_wdoublets, num_pcs=num_pcs, max_iter=0, var_genes=var_genes)
  
  ## Compute PC distance matrix
  print("Calculating PC distance matrix...")
  nCells <- nrow(seu_wdoublets@data.info)
  pca.coord <- seu_wdoublets@pca.rot[,PCs]
  rm(seu_wdoublets)
  gc()
  dist.mat <- fields::rdist(pca.coord)[,1:n.real.cells]
  
  ## Pre-order PC distance matrix prior to iterating across pK for pANN computations
  print("Defining neighborhoods...")
  for (i in 1:n.real.cells) {
    dist.mat[,i] <- order(dist.mat[,i])
  }

  ## Trim PC distance matrix for faster manipulations
  ind <- round(nCells * max(pK))+5
  dist.mat <- dist.mat[1:ind, ]

  ## Compute pANN across pK sweep
  print("Computing pANN across all pK...")
  for (k in 1:length(pK)) {
    print(paste("pK = ", pK[k], "...", sep = ""))
    pk.temp <- round(nCells * pK[k])
    pANN <- as.data.frame(matrix(0L, nrow = n.real.cells, ncol = 1))
    colnames(pANN) <- "pANN"
    rownames(pANN) <- real.cells
    list.ind <- list.ind + 1

    for (i in 1:n.real.cells) {
      neighbors <- dist.mat[2:(pk.temp + 1),i]
      pANN$pANN[i] <- length(which(neighbors > n.real.cells))/pk.temp
    }

    sweep.res.list[[list.ind]] <- pANN

  }

  return(sweep.res.list)
}
