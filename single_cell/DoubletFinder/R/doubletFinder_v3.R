doubletFinder_v3 <- function(seu, PCs, pN = 0.25, pK, nExp, num_pcs=25){
  seurat(); require(fields); require(KernSmooth)

  ## Make merged real-artifical data
  real.cells <- rownames(seu@data.info)
  data <- seu@raw.data[,real.cells]
  n_real.cells <- length(real.cells)
  n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
  print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[, real.cells1] + data[, real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data, doublets)

  ## Run Seurat workflow
  seu_wdoublets = run_seurat(name='test', write_out=F, dge=data_wdoublets, num_pcs=num_pcs, max_iter=0, var_genes=seu@var.genes)
  pca.coord = seu_wdoublets@pca.rot[, PCs]
  cell.names <- rownames(seu_wdoublets@data.info)
  nCells <- length(cell.names)
  
  ## Compute PC distance matrix
  print("Calculating PC distance matrix...")
  dist.mat <- fields::rdist(pca.coord)
  
  ## Compute pANN
  print("Computing pANN...")
  pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
  rownames(pANN) <- real.cells
  colnames(pANN) <- "pANN"
  k <- round(nCells * pK)
  for (i in 1:n_real.cells) {
    neighbors <- order(dist.mat[, i])
    neighbors <- neighbors[2:(k + 1)]
    neighbor.names <- rownames(dist.mat)[neighbors]
    pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
  }
  
  print("Classifying doublets..")
  classifications <- rep("Singlet",n_real.cells)
  classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
  seu@data.info[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@data.info), 1]
  seu@data.info[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
  return(seu)
}
