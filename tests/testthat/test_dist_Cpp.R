
test_that("Distance between embeddings is computed correctly", {
  set.seed(123)
  
  n.pts <- 30
  npcs <- 40
  n.cells <- 50
  stride <- 27
  n.cores <- 3
  
  emb.mat <- matrix(nrow = n.cells*npcs, ncol = n.pts)
  dist.mat <- matrix(nrow = (n.cells^2 - n.cells)/2, ncol = n.pts)
  
  for(i in 1:n.pts){
    emb <- matrix(data = runif(n.cells*npcs), nrow = n.cells)
    emb.mat[,i] <- c(emb)
    dist.mat[,i] <- dist(emb)
  }
  metric.R <- dist(t(dist.mat), method = "minkowski", p = 1 )/((n.cells^2 - n.cells)/2)

  emb <- bigmemory::as.big.matrix(emb.mat)
  emb.descr <- bigmemory::describe(emb)
  
  pairs.per.core <- ceiling( n.cells*(n.cells - 1)/(2*n.cores) )
  metric <- foreach(i = 1:n.cores,
                    .noexport = c("seurat", "gene.clusters", "points"),
                    .packages = c("moduli", "bigmemory"),
                    .inorder = F,
                    .combine =  "+",
                    .init = numeric((n.pts^2 - n.pts)/2))%do%{
                      
      start <- pairs.per.core*(i - 1) + 1
      end <- min(pairs.per.core*i, n.cells*(n.cells - 1)/2)
      emb <- attach.big.matrix(emb.descr)
      partial_dist_Cpp(emb@address, start, end, npcs, stride)
  }
  
  metric.matrix <- matrix(nrow = n.pts, ncol = n.pts)
  metric.matrix[lower.tri(metric.matrix)] <- metric/((n.cells^2 - n.cells)/2)
  metric <- as.dist(metric.matrix)
 
  expect_equal(metric, metric.R, ignore_attr = T)
})