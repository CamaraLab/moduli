
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

#' Enriches moduli object with its shared nearest neighbors graph
#' 
#' Finds a k-nearest neighbors graph representing the moduli space, computes the
#' Jaccard similarity matrix of this graph and uses this matrix as edge weights 
#' to construct a second graph that is saved in the moduli object.
#' 
#' @param moduli A moduli object
#' @param k Number of nearest neighbors
#' @param thld Weight threhold for graph pruning, edges with weight less than the given
#' value are removed before forming the snn graph. Default value is 0.0
#' 
#' @return A moduli object with an snn graph in the \code{snn.graph} slot.
#' 
#' @export
get_snn <- function(moduli, k, thld = 0.0){
  snn.graph <- .create_snn(moduli$metric, k, thld)
  moduli$snn.graph <- snn.graph
  return(moduli)
}

#' Retrieves Seurat object associated to a point in the moduli space
#' 
#' Runs PCA with features associated to the given point, also enriches
#' Seurat object with expression levels of each cluster in each cell
#' 
#' @param moduli A moduli object
#' @param point.id Id of point
#' 
#' @returns A Seurat object with the PCA analysis associated to the point in the
#' \code{reductions} slot and the average quantize-normalized expression level of
#' each gene cluster in each cell saved in the slot "gene_cliuster_id_exp" of the object's
#' meta data, where "id" is the id of the associated gene cluster.
#'
#' @export
retrieve_point <- function(moduli, point.id){
  seuratObject <- moduli$seurat
  gene.clusters <- unlist(modili$points$clusters[moduli$points$id == point.id])
  features <- unique(unlist(moduli$gene.clusters$genes[moduli$gene.clusters$id %in% gene.clusters]))
  approx  <- (length(features) > 2*moduli$npcs)
  seuratObject <- RunPCA(seuratObject, assay = moduli$assay, features = features,
                    npcs = moduli$npcs, approx = approx,
                    weight.by.var = moduli$weight.by.var, verbose = F)
  
  # enriching with gene cluster expression
  data <- GetAssayData(seuratObject[[moduli$assay]], slot = "scale.data")
  n.cells <- ncol(data)
  gene.expression.levels <- t(apply(data, 1, rank))
  gene.expression.levels <- gene.expression.levels/n.cells
  
  cls.expression.levels <- sapply(
    moduli$gene.clusters$genes,
    function(genes) colMeans(gene.expression.levels[genes, ,drop = F])
  )
  
  names <- paste0("gene_cluster_", moduli$gene.clusters$id,"_expr")
  for(i in 1:nrow(moduli$gene.clusters)){
    seuratObject[[names[i]]] <- cls.expression.levels[,i]
  }
  return(seuratObject)
}

#' Retries consensus metric
#' 
#' Computes the average distance between cell pairs across a set of embeddings
#' 
#' @param moduli A moduli object
#' @param point.ids Ids of moduli points to have their associated metrics averaged, default value is NULL
#' @param analysis.cluster.ids Ids of analysis clusters to have associated metrics averaged,
#' default value is NULL. The user must give either point ids or cluster analysis ids
#' @param desc.file Optional bigmatrix descriptor file with moduli embeddings as produced by
#' \link{get_moduli}. See 'Details' 
#' @param n.cores Numbers of cores to use, default value is 1
#' @param seed Random seed to use, default value is 42
#' @param filebacked Use a filebacked bigmatrix to store embeddings, only has an effect if
#' \code{desc.file} is NULL. Default value is FALSE
#' @param persist Keep bigmatrix with embeddings, only has effect if \code{filebacked} is TRUE.
#' Default value is FALSE
#' @param verbose Print time stamps, default value is TRUE
#' 
#' @returns A "dist" object with avereged cell distances and cell labels
#' 
#' @details 
#' If a descriptor file is given, pre-calculated embeddings will be used.
#' In this case the only cells represented in the output are the cells in \code{moduli$cells}.
#' If no descriptor file is given, embedding will be calculated again, but all cells in the
#' moduli's Seurat object will be represented in the output.
#' 
#' @export
retrieve_consesus <- function(moduli, point.ids = NULL, analysis.cluster.ids = NULL,
                              desc.file = NULL, n.cores = 1, seed = 42,
                              filebacked = F, persist = F, verbose = T){
  
  if(sum(c(is.null(point.ids), is.null(analysis.cluster.ids))) != 1){
    stop("Error: please give exatcly one, a set of point ids or a set of analysis cluster ids")
  }
  if(!is.null(analysis.cluster.ids)){
    point.ids <- unlist(
      moduli$analysis.clusters$points[moduli$analysis.clusters$id %in% analysis.cluster.ids]
    )
  }
  n.pts <- length(point.ids)
  emb.metric <- moduli$emb.metric
  npcs <- moduli$npcs
  weight.by.var <- moduli$weight.by.var
  
  if(is.null(desc.file)){
    
    cells <- Cells(moduli$seurat)
    scaled.data <- GetAssayData(moduli$seurat[[moduli$assay]], slot = "scale.data")
    gene.clusters <- moduli$gene.clusters
    points <- moduli$points
    
    if(filebacked){
      embeddings <- filebacked.big.matrix(
        nrow = length(cells)*npcs,
        ncol = n.pts,
        backingfile = "cons_embeddings.bin",
        descriptorfile = "cons_embeddings.desc"
      )
      if(!persist) on.exit(file.remove("cons_embeddings.bin", "cons_embeddings.desc"))
    } else {
      embeddings <- big.matrix(
        nrow = length(cells)*npcs,
        ncol = n.pts,
        shared = T
      )
    }
    emb.descr <- describe(embeddings)
    
    cl <- parallel::makeCluster(n.cores)
    registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = T)
    
    # Building embeddings
    if(verbose) message(paste(format(Sys.time(), "%a %b %d %X"),"Building embeddings"))
    
    foreach(
      i = 1:n.pts,
      .noexport = c("moduli", "cells"),
      .packages = c("bigmemory", "irlba")
    )%dopar%{
      pt.pos <- which(points$id == point.ids[i])
      feature.subset <- unique(unlist(gene.clusters$genes[points$clusters[[pt.pos]]]))
      cell.embeddings <- .pca_aux(scaled.data, feature.subset, npcs, weight.by.var = weight.by.var,
                                  scale = emb.metric, seed = seed)
      emb <- attach.big.matrix(emb.descr)
      emb[,i] <- c(cell.embeddings)
    }
    parallel::stopCluster(cl)
    
  } else {
    emb.descr <- desc.file
    cells <- moduli$cells
  }
  
  
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  
  # Building consensus
  if(verbose) message(paste(format(Sys.time(), "%a %b %d %X"),"Computing consensus metric"))
  n.cells <- length(cells)
  n.cell.pairs <- n.cells*(n.cells - 1)/2
  pairs.per.core <- ceiling(n.cell.pairs/n.cores)
  
  consensus <- foreach(
    i = 1:n.cores,
    .noexport = c("seurat", "gene.clusters", "points"),
    .packages = c("moduli", "bigmemory"),
    .inorder = T,
    .combine =  "c"
  )%dopar%{
    start <- pairs.per.core*(i - 1) + 1
    end <- min(pairs.per.core*i, n.cell.pairs)
    emb <- attach.big.matrix(emb.descr)
    consensus_dist(emb@address, emb.metric, start, end, npcs)
  }
  
  # formatting output
  consensus.metric.matrix <- matrix(nrow = n.cells, ncol = n.cells)
  consensus.metric.matrix[lower.tri(consensus.metric.matrix)] <- consensus
  rownames(consensus.metric.matrix) <- cells
  colnames(consensus.metric.matrix) <- cells
  
  return(as.dist(consensus.metric.matrix))
}

# Runs PCA and scales result depending on the embedding metric that will be used 
# (a translation and homotety for "euclidean" and a projection to the sphere if "cosine")

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
    pc.means <- colMeans(cell.embeddings)
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

# creates a knn graph, computes its Jaccard similarity matrix and uses this matrix
# to define a new graph

.create_snn <- function(dist, k, thld = 0.0){
  knn <- dbscan::kNN(dist, k)
  id <- knn$id
  n.vertices <- nrow(id)
  rows <- 1:n.vertices
  neighbors <- matrix(0, nrow = n.vertices, ncol = n.vertices)
  
  for(col in 1:k){
    neighbors[cbind(rows, id[,col])] <- 1
    neighbors[cbind(id[,col], rows)] <- 1
  }
  
  knn.graph <- igraph::graph_from_adjacency_matrix(neighbors, mode = "undirected")
  jaccard.sim <- igraph::similarity.jaccard(knn.graph, loops = T)
  jaccard.sim[jaccard.sim < thld] <- 0.0

  snn.graph <- igraph::graph_from_adjacency_matrix(
    jaccard.sim,
    mode = "undirected",
    weighted = T,
    diag = F
  )
  return(snn.graph)
}


