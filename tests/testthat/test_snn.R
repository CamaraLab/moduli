
test_that("SNN graph is constructed correctly",{
  
  d <- as.dist(rbind(
    c(0,  0,  0, 0),
    c(1,  0,  0, 0),
    c(2, 1/2, 0, 0),
    c(3,  3, 1/4, 0)
  ))
  # knn matrix fro k = 1 is 1 -- 2 -- 3 -- 4
  
  exp.adj <- rbind(
    c(0,   2/3,  1/4,  0),
    c(2/3,   0,  1/2,  1/4),
    c(1/4, 1/2,  0,    2/3),
    c(0,   1/4,  2/3,  0)
  )
  
  graph <- .create_snn(d, 1)
  expect_equal(igraph::as_adjacency_matrix(graph, attr = "weight", sparse = F),
               exp.adj)
})