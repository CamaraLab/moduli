#' Computes moduli space
#' 
#' Here, by moduli space we mean a metric space where each point represents an embedding
#' of cells by PCA with a different choice of features, and each choice of features is a
#' combination of gene clusters. By default, computes PCA embeddings of all nonempty 
#' combinations of gene clusters.
#' 
#' 
#' @param seuratObject A Seurat object
#' @param gene.membership An integer vector with gene cluster assignments ans genes as names
#' @param gene.clusters Clusters of genes as a list of character vectors, these clusters may
#' intersect (the user must give exactly one the arguments gene.membership or gene.clusters)
#' @param assay Assay in which to run PCA (DefaultAssay(seuratObject) by default). The assay must have a
#' "scale.data" slot
#' @param npcs Number of principal components (defaults to 30)
#' @param weight.by.var Weight the cell embedding by variance of each PC (defaults to True)
#' @param point Points to be considered as a list of integer vectors. If NULL all nonempty subsets of the
#' set of clusters will be considered (defaults to NULL)
#' @param n.cores Number of cores to use
#' @param n.pts Number of points to subsample, if 0 all points are taken (default is 0)
#' @param n.cells Number of cells to subsample when computing distances in the moduli space, if 0 all
#' cells are used (default is 0).
#' @param stride Number of pairs of cells to use in one iteration when computing the moduli metric. Taking
#' larger strides is more efficient but requires more memory (default is 10000)
#' @param seed Random seed to use (default is 42)
#' @param filebacked Use a filebacked bigmatrix to store embeddings (default is TRUE)
#' @param persist Keep bigmatrix with embeddings (default is FALSE)
#' @param verbose Print time stamps (default is TRUE)
#' 
#' @return An object of class moduli with the following components:
#' \item{gene.clusters}{A data.frame with gene clusters}
#' \item{points}{A data.frame with point metadata}
#' \item{metric}{A 'dist' object with the moduli metric}
#' \item{seurat}{The Seurat object used to contruct the moduli}
#' \item{assay}{Name of the assay used to construct the moduli}
#' \item{cells}{Character vector with names of cells used to compute metric}
#' \item{npcs}{Number of principal components used in moduli construction}
#' 
#' @details Computing the moduli space is very memory and processing intensive, specially
#' computing the distances in the moduli space. The time complexity of the algorithm is
#' \deqn{O(n.cells^2 n.pts^2 npcs + n.cells^2 n.pts^2)}
#' and by default \eqn{n.pts = 2^{n.clusters} - 1}, where \eqn{n.clusters} is the number of  gene clusters.
#' The matrix with embeddings has about \eqn{8*n.pts*npcs*n.cells*10^{-9}} GB, and in
#' computing distances, each core will use need about \eqn{8*stride*n.pts*10^{-9}} GB.
#' 
#' @examples
#' library(Seurat)
#' data("pbmc_small")
#' 
#' # normalize data
#' pbmc_small <- SCTransform(pbmc_small, verbose = F)
#' 
#' # cluster genes
#' pam <- gene_medioid_clustering(pbmc_small, 4)
#' 
#' # compute moduli
#' moduli <- get_moduli(pbmc_small, membership = pam$clustering, npcs = 5, verbose = F, filebacked = F)
#' 
#' 
#' @export
get_moduli <- function(seuratObject, gene.membership = NULL, gene.clusters = NULL,
                       assay = DefaultAssay(seuratObject), npcs = 30, weight.by.var = T,
                       points = NULL, n.cores = 1, n.pts = 0, n.cells = 0, stride = 10000,
                       seed = 42, filebacked = T, persist = F, verbose = T){
  # gene.clusters : clusters of genes as a list of vectors of gene names
  if(sum(c(is.null(gene.membership), is.null(gene.clusters))) != 1){
    stop("Error: provide exatcly one of the arguments gene.membership or gene.clusters")
  }
  
  if(!is.null(gene.membership)){
    gene.clusters <- list()
    cluster.ids <- sort(unique(gene.membership))
    gene.clusters <- list()
    for(k in cluster.ids){
      gene.clusters[[k]] <- names(gene.membership[gene.membership == k])
    }
  }else{
    cluster.ids <- seq_along(gene.clusters)
  }
  n.clusters <- length(gene.clusters)
  temp <- gene.clusters
  gene.clusters <- data.frame(id = cluster.ids)
  gene.clusters$genes <- temp
  
  set.seed(seed)
  if(is.null(points) && (n.pts == 0 || n.pts >= 2^n.clusters -1 )){
    for(n in 1:n.clusters){
      comb <- combn(gene.clusters$id, n, simplify = F)
      points <- append(points, comb)
    }
  } else if(is.null(points)){
    ints <- sample(1:(2^n.clusters -1), n.pts)
    points <- lapply(ints, function(d) gene.clusters$id[intToBits(d) == 1])
  }
  
  if( n.pts > 0 && length(points) > n.pts) points <- sample(points, n.pts)
  n.pts <- length(points)
  temp <- points
  points <- data.frame(id = 1:n.pts)
  points$clusters <- temp
  
  if(n.cells > 0 && length(Cells(seuratObject)) > n.cells){
    cells <- sample(Cells(seuratObject), n.cells) 
  } else {
    cells <- Cells(seuratObject)
    n.cells <- length(cells)
  }
  
  scaled.data <- GetAssayData(seuratObject[[assay]], slot = "scale.data" )[,cells]
  
  if(filebacked){
    embeddings <- filebacked.big.matrix(
      nrow = length(cells)*npcs,
      ncol = n.pts,
      backingfile = "embeddings.bin",
      descriptorfile = "embeddings.desc"
    )
    if(!persist) on.exit(file.remove("embeddings.bin", "embeddings.desc"))
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
  
  if(verbose) message(paste(format(Sys.time(), "%a %b %d %X"),"Building embeddings"))
  foreach(i = 1:n.pts,
          .noexport = c("seuratObject", "cells"),
          .packages = c("bigmemory", "irlba"))%dopar%{
    feature.subset <- unique(unlist(gene.clusters$genes[points$clusters[[i]]]))
    
    #based on Seurat's RunPCA
    if (length(feature.subset) > npcs*2){
      set.seed(seed)
      pca.results <- irlba(t(scaled.data[feature.subset,]), nv = npcs)
      if(weight.by.var){
        cell.embeddings <- pca.results$u %*% diag(pca.results$d)
      } else {
        cell.embeddings <- pca.results$u
      }
    } else {
      rank <- min(npcs, nrow(scaled.data))
      pca.results <- prcomp(t(scaled.data[feature.subset,]), rank. = rank)
      if(weight.by.var){
        cell.embeddings <- pca.results$x
      } else {
        cell.embeddings <- pca.results$x / (pca.results$sdev[1:rank] * sqrt(nrow(cell.embeddings) - 1))
      }
      
      cell.embeddings <- cbind(
        cell.embeddings,
        matrix(data = 0, nrow = nrow(cell.embeddings), ncol = npcs - rank)
      )
    }
    #standardizing embedding
    pc.means <- apply(cell.embeddings, 2, mean)
    cell.embeddings <- apply(cell.embeddings, 1, function(row) row - pc.means)
    total.var <- sum(cell.embeddings^2)/(nrow(cell.embeddings) - 1)
    cell.embeddings <- cell.embeddings*sqrt(npcs/total.var)
    
    #saving to bigmatrix
    emb <- attach.big.matrix(emb.descr)
    emb[,i] <- c(cell.embeddings)
  }
  
  parallel::stopCluster(cl)
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  
  if(verbose) message(paste(format(Sys.time(), "%a %b %d %X"), "Computing distances"))
  pairs.per.core <- ceiling( n.cells*(n.cells - 1)/(2*n.cores) )
  metric <- foreach(i = 1:n.cores,
                      .noexport = c("seuratObject", "gene.clusters", "points"),
                      .packages = c("moduli", "bigmemory"),
                      .inorder = F,
                      .combine =  "+",
                      .init = numeric((n.pts^2 - n.pts)/2))%dopar%{
    
    start <- pairs.per.core*(i - 1) + 1
    end <- min(pairs.per.core*i, n.cells*(n.cells - 1)/2)
    emb <- attach.big.matrix(emb.descr)
    partial_dist_Cpp(emb@address, start, end, npcs, stride)
  }
  if(verbose) message(paste(format(Sys.time(), "%a %b %d %X"), "Finished computing distances"))
  
  metric.matrix <- matrix(nrow = n.pts, ncol = n.pts)
  metric.matrix[lower.tri(metric.matrix)] <- metric/((n.cells^2 - n.cells)/2)
  metric <- as.dist(metric.matrix)
  
  out <- list(
    gene.clusters = gene.clusters,
    points = points,
    metric = metric,
    seurat = seuratObject,
    assay = assay,
    cells = cells,
    npcs = npcs
  )
  
  class(out) <- "moduli"
  return(out)
}



