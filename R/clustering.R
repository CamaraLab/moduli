
create_snn <- function(dist, k){
  snn <- dbscan::kNN(dist, k)
  id <- snn$id
  n.vertices <- nrow(id)
  rows <- 1:n.vertices
  neighbors <- matrix(0, nrow = n.vertices, ncol = n.vertices)
  
  for (col in 1:k) neighbors[cbind(rows, id[,col])] <- 1
  
  diag(neighbors) <- 1
  
  adjacency <- neighbors %*% t(neighbors)
  adjacency[adjacency > 0] <- adjacency[adjacency > 0] / (2*k + 2 - adjacency[adjacency > 0])
  snn.graph <- igraph::graph_from_adjacency_matrix(
    adjacency,
    mode = "undirected",
    weighted = T,
    diag = F
  )
  
  return(snn.graph)
}

#' @export
gene_medioid_clustering <- function(seuratObject, n.clusters, slot = "scale.data", 
                             features = VariableFeatures(seuratObject)){
  data <- FetchData(seuratObject, vars = features, slot = slot)
  dist.mat <- 1 - cor(data)
  part <- cluster::pam(dist.mat, n.clusters, diss = T)
  return(part)
}

#' @export
cluster_moduli_space <- function(moduli, k, partition_type = "RBConfigurationVertexPartition",
                                 resolution_parameter = 1, save.snn = T, enrich = T, seed = 123){
  snn.graph <- create_snn(moduli$metric, k)
  membership <- leiden::leiden(
    snn.graph,
    partition_type = partition_type,
    resolution_parameter = resolution_parameter,
    seed = seed
  )
  out <- moduli
  out$analysis.clusters <- data.frame(id = sort(unique(membership)))
  out$analysis.clusters$points <- lapply(out$analysis.clusters$id,
                                           function(p) out$points$id[membership == p])
  
  if(save.snn) out$snn.graph <- snn.graph
  
  if(enrich) out <- enrich_moduli_clusters(out)
  
  return(out)
}


#' @export
enrich_moduli_clusters <- function(moduli, thld = 0.05){
  
  g.cluster.table <- matrix(data = 0, nrow = nrow(moduli$analysis.clusters),
                            ncol = nrow(moduli$gene.clusters))
  a.cluster.sizes <- integer(nrow(moduli$analysis.clusters))
  
  for(i in 1:nrow(moduli$analysis.clusters)){
    members <- moduli$analysis.clusters$points[[i]]
    tab <- table(unlist(moduli$points$clusters[moduli$points$id %in% members]))
    for(j in seq_along(tab)){
      g.cluster.table[i, moduli$gene.clusters$id %in% names(tab)[j]] <- tab[j]
    }
    a.cluster.sizes[i] <- length(members)
  }
  
  g.cluster.freq <- apply(g.cluster.table, 2, sum)/nrow(moduli$points)
  exp.gene.clusters <- NULL
  for(i in 1:nrow(moduli$analysis.clusters)){
    exp.gene.clusters[[i]] <- integer()
    for(j in 1:ncol(g.cluster.table)){
      if((g.cluster.freq[j] == 0) || (g.cluster.freq[j] == 1)) next
      p.val <- pbinom(g.cluster.table[i,j] -1, a.cluster.sizes[i], g.cluster.freq[j], lower.tail = F)
      if(p.val <= thld){
        exp.gene.clusters[[i]] <- append(exp.gene.clusters[[i]], moduli$gene.clusters$id[j])
      }
    }
  }
  out <- moduli
  out$analysis.clusters$exp.gene.clusters <- exp.gene.clusters
  return(out)
}



