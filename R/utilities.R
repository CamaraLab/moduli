
# Remove eventually
#' @export
update_moduli <- function(moduli){
  out <- moduli
  out$gene.clusters <- data.frame(id = seq_along(moduli$gene.clusters))
  out$gene.clusters$genes <- moduli$gene.clusters
  
  out$points <- data.frame(id = seq_along(moduli$points))
  out$points$clusters <- moduli$points
  return(out)
}

#' Retrieve ids of points
#' 
#' Retrieves ids of points using a set of gene clusters.
#' 
#' @param moduli A moduli object
#' @param gene.clusters vector with ids of gene clusters
#' @param filter a choice of filter, default value is \code{"equal"}. See details.
#' 
#' @return A vector with point ids
#' 
#' @details The filter \code{"equal"} gets the point that represents the point representing
#' the clusters in \code{gene.clusters}. The option \code{"contains"} will get all points that
#' that representing subsets of \code{gene.clusters}, while \code{"is_subset"} will get all
#' points that have all elements of \code{gene.clusters}. The filter \code{"intersetcs"} retrieves
#' ids of points that have any gene cluster in \code{gene.clusters}, and \code{"in_complement"}
#' does the opposite, retrieving ids of points that only represent clusters not in \code{gene.clusters}.
#' 
#' @export
which_point <- function(moduli, gene.clusters,
                        filter = c("equal", "contains", "is_subset", "intersects", "in_complement")){
  filter <- match.arg(filter)
  f <- function(pt.clusters){
    switch (filter,
            equal = all(pt.clusters %in% gene.clusters) && all(gene.clusters %in% pt.clusters),
            contains = all(pt.clusters %in% gene.clusters),
            is_subset = all(gene.clusters %in% pt.clusters),
            intersects = any(pt.clusters %in% gene.clusters),
            in_complement = !any(pt.clusters %in% gene.clusters)
    )
  }
  return(moduli$points$id[sapply(moduli$points$clusters, f)])
}

#' Retrieve point metadata
#' 
#' Retrieves metadata of moduli point
#' 
#' @param moduli A moduli object
#' @param points vector with ids of points to be retrieved, will retrieve all points by default
#' 
#' @return A data.frame with point ids, gene clusters represented and, if available,
#' the analysis cluster the point belongs to
#' 
#' @export
point_metadata <- function(moduli, points = moduli$points$id){
  out <- moduli$points[moduli$points$id %in% points,]
  if(!is.null(moduli$analysis.clusters)){
    out$analysis.cluster <- integer(nrow(analysis.cluster))
    for(i in 1:nrow(moduli$analysis.clusters)){
      out$analysis.cluster[out$id == moduli$analysis.clusters$points[i]] <- moduli$analysis.clusters$id[i]
    }
  }
  return(out)
}

#' Enriches moduli with its shared nearest neighbors graph
#' 
#' @param moduli A moduli object
#' @param k Number of nearest neighbors
#' @param thld Minimum of edges weights to be kept
#' 
#' @return A moduli object with an snn graph in the \code{snn.graph} slot.
#' @export
get_snn <- function(moduli, k, thld = 0.0){
  snn.graph <- .create_snn(moduli$metric, k, thld)
  out <- moduli
  out$snn.graph <- snn.graph
  return(out)
}



# keep?
#' @export
run_pca <- function(moduli, point){
  seuratObject <- moduli$seurat
  gene.clusters <- unlist(modili$points$clusters[moduli$points$id == point])
  features <- unique(unlist(moduli$gene.clusters$genes[moduli$gene.clusters$id %in% gene.clusters]))
  approx  <- (length(features) > 2*moduli$npcs)
  seuratObject <- RunPCA(seuratObject, assay = moduli$assay, features = features,
                    npcs = moduli$npcs, approx = approx, verbose = F)
  return(suratObject)
}


.pca_aux <- function(scaled.data, feature.subset, npcs, weight.by.var= T, scale = NULL, seed = 42){
  set.seed(seed)
  # based on Seurat's RunPCA
  if (length(feature.subset) > npcs*2){
    pca.results <- irlba(t(scaled.data[feature.subset,]), nv = npcs)
    if(weight.by.var){
      cell.embeddings <- pca.results$u %*% diag(pca.results$d)
    } else {
      cell.embeddings <- pca.results$u
    }
  } else {
    rank <- min(npcs, length(feature.subset))
    pca.results <- prcomp(t(scaled.data[feature.subset,]), rank. = rank)
    if(weight.by.var){
      cell.embeddings <- pca.results$x
    } else {
      cell.embeddings <- pca.results$x / (pca.results$sdev[1:rank] * sqrt(nrow(cell.embeddings) - 1))
    }
    # add zeros to make dimension equal to npcs if necessary
    cell.embeddings <- cbind(
      cell.embeddings,
      matrix(data = 0, nrow = nrow(cell.embeddings), ncol = npcs - rank)
    )
  }
  if(is.null(scale)) return(cell.embeddings)
  if(scale == "euclidean"){
    pc.means <- apply(cell.embeddings, 2, mean)
    cell.embeddings <- t(t(cell.embeddings) - pc.means)
    total.var <- sum(cell.embeddings^2)/(nrow(cell.embeddings) - 1)
    cell.embeddings <- cell.embeddings*sqrt(npcs/total.var)
    return(cell.embeddings)
  }
  if(scale == "cosine"){
    cell.embeddings <- cell.embeddings / sqrt(apply(cell.embeddings^2, 1, sum))
    return(cell.embeddings)
  }
}

.create_snn <- function(dist, k, thld = 0.0){
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


