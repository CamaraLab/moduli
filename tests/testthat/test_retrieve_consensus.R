
test_that("Consesus metric is correct", {
  data("pbmc_small_moduli")
  set.seed(123)
  pts <- which_point(pbmc_small_moduli, c(1), filter = "is_subset")
  
  consensus_metric_R <- function(moduli, point.ids){
    cells <- Cells(moduli$seurat)
    scaled.data <- GetAssayData(moduli$seurat[[moduli$assay]], slot = "scale.data")
    cons <- numeric(length(cells)*(length(cells) -1)/2)
    
    for(i in seq_along(point.ids)){
      pt.pos <- which(moduli$points$id == point.ids[i])
      cls.in.pt <- moduli$points$clusters[[pt.pos]]
      feature.subset <- unique(unlist(moduli$gene.clusters$genes[moduli$gene.clusters$id %in% cls.in.pt]))
      cell.embeddings <- .pca_aux(scaled.data, feature.subset, moduli$npcs, scale = moduli$emb.metric)
      if(moduli$emb.metric == "euclidean"){
        cons <- cons + dist(cell.embeddings)
      } else if(moduli$emb.metric == "cosine"){
        temp <- cell.embeddings %*% t(cell.embeddings)
        temp <- temp[lower.tri(temp)]
        cons <- cons + acos(temp)
      }
    }
    return(cons/length(point.ids))
  }
  
  R.result <- consensus_metric_R(pbmc_small_moduli, pts)
  result <- retrieve_consensus(pbmc_small_moduli, point.ids = pts, verbose = F)
  expect_equal(result, R.result, ignore_attr = T)
  
  
})


