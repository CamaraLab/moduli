library("moduli")

test_that("Distance between embeddings is computed correctly", {
  set.seed(123)
  
  nemb <- 30
  npcs <- 40
  npts <- 50
  
  emb.mat <- matrix(nrow = npts*npcs, ncol = nemb)
  dist.mat <- matrix(nrow = (npts^2 - npts)/2, ncol = nemb)
  
  for(i in 1:nemb){
    emb <- matrix(data = runif(npts*npcs), nrow = npts)
    emb.mat[,i] <- c(emb)
    dist.mat[,i] <- dist(emb)
  }
  dist.moduli.R <- dist(t(dist.mat), method = "minkowski", p = 1 )/nrow(dist.mat)

  bm.emb <- bigmemory::as.big.matrix(emb.mat)
  dist.moduli.Cpp <- matrix(nrow = nemb, ncol = nemb)
  for(i in 1:nemb){
    dist.moduli.Cpp[,i] <- dist_Cpp(bm.emb@address, i, npcs)
  }
  expect_equal(dist.moduli.Cpp + t(dist.moduli.Cpp), as.matrix(dist.moduli.R), ignore_attr = T)
})