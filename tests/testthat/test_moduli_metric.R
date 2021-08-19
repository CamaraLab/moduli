
test_that("Distance between embeddings is computed correctly", {
  set.seed(123)
  
  data("pbmc_small_moduli")
 
  pbmc_small <- pbmc_small_moduli$seurat
  gene.clusters <- pbmc_small_moduli$gene.clusters
  npcs <- 3
  
  # straightforward calculation
  get_metric_R <- function(seurat, gene.clusters, npcs, emb.metric, moduli.metric){
    n.clusters <- nrow(gene.clusters)
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
      feature.subset <- unique(unlist(gene.clusters$genes[gene.clusters$id %in% points[[i]]]))
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
      metric.R <- get_metric_R(pbmc_small, gene.clusters, npcs, emb.metric, moduli.metric)
      metric <- get_moduli(pbmc_small, gene.clusters = gene.clusters$genes, npcs = npcs,
                           emb.metric = emb.metric, moduli.metric = moduli.metric,
                           filebacked = F, verbose = F)$metric
      expect_equal(metric, metric.R, ignore_attr = T)
    }
  }
})