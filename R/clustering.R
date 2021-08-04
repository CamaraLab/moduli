library(igraph)
library(cluster)
#library(leiden)
library(uwot)

create_knn <- function(dist, k){
  dist <- as.matrix(dist)
  n.vertices <- nrow(dist)
  edge_list <- matrix(data = -1, nrow = k*n.vertices, ncol = 2)
  for(i in 1:n.vertices){
    neighbors <- order(dist[i,])[2:(k+1)]
    for(j in 1:k){
      # not very efficient but works
      if(all(edge_list[,1] != min(i, neighbors[j]) | edge_list[,2] != max(i, neighbors[j]))){
        edge_list[(i-1)*k + j, 1] <- min(i, neighbors[j])
        edge_list[(i-1)*k + j, 2] <- max(i, neighbors[j])
      }
    }
  }
  graph <- graph_from_edgelist(edge_list[edge_list[,1] > 0, ], directed = F)
  return(graph)
}

modularity_partition <- function(seurat, k = 10, slot = "scale.data",
                                 features = VariableFeatures(seurat),
                                 resolution = 1,
                                 verbose = T){
  data <- FetchData(seurat, vars = features, slot = slot)
  dist.mat <- 1 - cor(data)
  knn.graph <- create_knn(dist.mat,k)
  jaccard <- similarity(knn.graph, method = "jaccard", loops = T)
  jaccard.graph <- graph_from_adjacency_matrix(jaccard, mode = "undirected", weighted = T, diag = F)
  
  membership <- igraph::cluster_louvain(jaccard.graph)$membership
  #membership <- leiden(jaccard.graph, partition_type = "RBConfigurationVertexPartition",
  #                     resolution_parameter = resolution)
  names(membership) <- features
  
  if(verbose){
    print("Modularity:")
    print(modularity(jaccard.graph, membership, weights = E(jaccard.graph)$weight))
  }
  return(membership)
}

medioid_partition <- function(seurat, n.clusters, slot = "scale.data", 
                             features = VariableFeatures(seurat)){
  data <- FetchData(seurat, vars = features, slot = slot)
  dist.mat <- 1 - cor(data)
  part <- cluster::pam(dist.mat, n.clusters, diss = T)
  return(part)
}

visualize_partition <- function(seurat, membership, slot = "scale.data", title = NULL){
  set.seed(123)
  data <- t(FetchData(seurat, vars = names(membership), slot = slot))
  umap.coords <- umap(data, metric = "correlation")
  plot(umap.coords, col = membership, main = title)
  
}

