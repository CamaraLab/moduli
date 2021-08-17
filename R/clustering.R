
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
                                 resolution_parameter = 1, save.snn = T, enrich = T,
                                 thld = 0.05, seed = 123){
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
  
  if(enrich) out <- enrich_moduli_clusters(out, thld)
  
  return(out)
}


#' @export
enrich_moduli_clusters <- function(moduli, thld = 0.05){
  # total occurrence of each gene cluster
  g.clt.totals <- integer(nrow(moduli$gene.clusters))
  tab <- table(unlist(moduli$points$clusters))
  ord <- order(moduli$gene.clusters$id)
  g.clt.totals[ord[ord %in% names(tab)]] <- tab
  
  
  exp.gene.clusters <- NULL
  enrichment.p.vals <- NULL
  for(i in 1:nrow(moduli$analysis.clusters)){
    
    # counts of occurrences of gene clusters in the analysis cluster
    g.clt.counts <- integer(nrow(moduli$gene.clusters))
    tab <- table(
      unlist(moduli$points$clusters[moduli$points$id %in% moduli$analysis.clusters$points[[i]]])
    )
    g.clt.counts[ord[ord %in% names(tab)]] <- tab
    
    a.clt.size <- length(moduli$analysis.clusters$points[[i]])
    p.vals <- numeric(nrow(moduli$gene.clusters))
    
    for(j in 1:nrow(moduli$gene.clusters)){
      # table         g.clst in pt
      #                 Y   N
      # pt in a.clst Y
      #              N
      cont.table <- rbind(
        c(g.clt.counts[j]                  , a.clt.size - g.clt.counts[j]),
        c(g.clt.totals[j] - g.clt.counts[j], nrow(moduli$points) - g.clt.totals[j] - a.clt.size + g.clt.counts[j])
      )
      p.vals[j] <- fisher.test(cont.table, alternative = "greater")$p.value
    }
    p.vals <- p.adjust(p.vals, method = "BH")
    exp.gene.clusters[[i]] <- moduli$gene.clusters$id[p.vals < thld]
    enrichment.p.vals[[i]] <- p.vals[p.vals < thld]
    # ordering by p-value
    exp.gene.clusters[[i]] <- exp.gene.clusters[[i]][order(enrichment.p.vals[[i]])]
    enrichment.p.vals[[i]] <- sort(enrichment.p.vals[[i]])
  }
  out <- moduli
  out$analysis.clusters$exp.gene.clusters <- exp.gene.clusters
  out$enrichment.p.vals <- enrichment.p.vals
  return(out)
}



