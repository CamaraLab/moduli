library(Seurat)

test_that("Distance between embeddings is computed correctly", {
  set.seed(123)
  
  data("pbmc_small")
  n.clusters <- 4
  npcs <- 3
  
  oldw <- getOption("warn")
  options(warn = -1)
  pbmc_small <- SCTransform(pbmc_small, verbose = F, variable.features.n = 100)
  options(warn = oldw)
  
  vf <- VariableFeatures(pbmc_small)
  
  membership <- ceiling(runif(length(vf), min = 0, max = n.clusters))
  names(membership) <- vf
  
  # straightforward calculation
  get_metric_R <- function(seurat, membership, npcs, emb.metric, moduli.metric){
    n.clusters <- length(unique(membership))
    n.cells <- length(Cells(seurat))
    n.cell.pairs <- n.cells*(n.cells -1)/2
    
    points <- NULL
    for(n in 1:n.clusters){
      comb <- combn(1:n.clusters, n, simplify = F)
      points <- append(points, comb)
    }
    n.pts <- length(points)
    emb.mat <- matrix(nrow = n.cells*npcs, ncol = n.pts)
    dist.mat <- matrix(nrow = n.cell.pairs, ncol = n.pts)
    
    scaled.data <- GetAssayData(seurat, slot = "scale.data")
  
    for(i in 1:n.pts){
      feature.subset <- names(membership)[membership %in% points[[i]]]
      cell.embeddings <- .pca_aux(scaled.data, feature.subset, npcs, scale = emb.metric)
      emb.mat[,i] <- c(cell.embeddings)
      if(emb.metric == "euclidean"){
        dist.mat[,i] <- dist(cell.embeddings)
      } else if(emb.metric == "cosine"){
        temp <- cell.embeddings %*% t(cell.embeddings)
        temp <- temp[lower.tri(temp)]
        dist.mat[,i] <- c(acos(temp))
      }
    }
    if(moduli.metric == "l1"){
      metric.R <- dist(t(dist.mat), method = "minkowski", p = 1 )/n.cell.pairs
    } else if(moduli.metric == "correlation"){
      metric.R <- as.dist(1 - cor(dist.mat))
    }
    return(metric.R)
  }
  
  for(emb.metric in c("cosine", "euclidean")){
    for(moduli.metric in c("correlation", "l1")){
      metric.R <- get_metric_R(pbmc_small, membership, npcs, emb.metric, moduli.metric)
      metric <- get_moduli(pbmc_small, gene.membership = membership, npcs = npcs,
                           emb.metric = emb.metric, moduli.metric = moduli.metric,
                           filebacked = F, verbose = F)$metric
      expect_equal(metric, metric.R, ignore_attr = T)
    }
  }
})