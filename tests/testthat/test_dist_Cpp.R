
test_that("Distance between embeddings is computed correctly", {
  set.seed(123)
  
  n.pts <- 30
  npcs <- 40
  n.cells <- 50
  stride <- 53
   
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
  
  metric <- foreach(i = 1:(ceiling((n.cells^2 - n.cells)/(2*stride))),
                    .inorder = F,
                    .combine =  "+",
                    .init = numeric((n.pts^2 - n.pts)/2))%do%{
    start <- stride*(i - 1) + 1
    end <- min(stride*i, (n.cells^2 - n.cells)/2)
    emb <- bigmemory::attach.big.matrix(emb.descr)
    partial_dist_Cpp(emb@address, start, end, npcs)
  }
  
  metric.matrix <- matrix(nrow = n.pts, ncol = n.pts)
  metric.matrix[lower.tri(metric.matrix)] <- metric/((n.cells^2 - n.cells)/2)
  metric <- as.dist(metric.matrix)
 
  expect_equal(metric, metric.R, ignore_attr = T)
})