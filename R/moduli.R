#' @useDynLib moduli
#' @importFrom Rcpp sourceCpp
#' @import Seurat
#' @import foreach
#' @import doParallel
#' @import irlba
#' @import bigmemory
#' 
#'
#' @export
get_moduli <- function(seurat, membership = NULL, gene.clusters = NULL,
                       assay = DefaultAssay(seurat), npcs = 30, weight.by.var = T,
                       points = NULL, n.cores = 1, n.pts = 0, n.cells = 0, stride = 10000,
                       seed = 123, filebacked = T, persist = T){
  # gene.clusters : clusters of genes as a list of vectors of gene names
  
  if(!is.null(membership)){
    gene.clusters <- list()
    cluster.labels <- sort(unique(membership))
    gene.clusters <- list()
    for(k in cluster.labels){
      gene.clusters[[k]] <- names(membership[membership == k])
    }
  }
  
  cluster.labels <- names(gene.clusters)
  if(is.null(cluster.labels)) cluster.labels <- 1:length(gene.clusters)
  
  set.seed(seed)
  if(is.null(points)){
    for(n in 1:length(gene.clusters)){
      comb <- combn(cluster.labels, n, simplify = F)
      points <- append(points, comb)
    }
  }
  
  if( n.pts > 0 && length(points) > n.pts) sample(points, n.pts)
  n.pts <- length(points)
  
  if(n.cells > 0 && length(Cells(seurat)) > n.cells){
    cells <- sample(Cells(seurat), n.cells) 
  } else {
    cells <- Cells(seurat)
    n.cells <- length(cells)
  }
  
  scaled.data <- GetAssayData(seurat[[assay]], slot = "scale.data" )[,cells]
  
  if(filebacked){
    embeddings <- filebacked.big.matrix(
      nrow = length(cells)*npcs,
      ncol = length(points),
      backingfile = "embeddings.bin",
      descriptorfile = "embeddings.desc"
    )
    if(!persist) on.exit(file.remove("embeddings.bin", "embeddings.desc"))
  } else {
    embeddings <- big.matrix(
       nrow = length(cells)*npcs,
       ncol = length(points),
       shared = T
    )
  }
  
  emb.descr <- describe(embeddings)
  
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = T)
  
  message(paste(format(Sys.time(), "%a %b %d %X"),"Building embeddings"))
  foreach(i = seq_along(points),
          .noexport = c("seurat", "cells"),
          .packages = c("bigmemory", "irlba"))%dopar%{
    feature.subset <- unique(unlist(gene.clusters[points[[i]]]))
    
    #based on Seurat's RunPCA
    if (length(feature.subset) > npcs*2){
      set.seed(42) #seed used in Seurat
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
  
  message(paste(format(Sys.time(), "%a %b %d %X"), "Computing distances"))
  pairs.per.core <- ceiling( n.cells*(n.cells - 1)/(2*n.cores) )
  metric <- foreach(i = 1:n.cores,
                      .noexport = c("seurat", "gene.clusters", "points"),
                      .packages = c("moduli", "bigmemory"),
                      .inorder = F,
                      .combine =  "+",
                      .init = numeric((n.pts^2 - n.pts)/2))%dopar%{
    
    start <- pairs.per.core*(i - 1) + 1
    end <- min(pairs.per.core*i, n.cells*(n.cells - 1)/2)
    emb <- attach.big.matrix(emb.descr)
    partial_dist_Cpp(emb@address, start, end, npcs, stride)
  }
  message(paste(format(Sys.time(), "%a %b %d %X"), "Finished computing distances"))
  
  metric.matrix <- matrix(nrow = n.pts, ncol = n.pts)
  metric.matrix[lower.tri(metric.matrix)] <- metric/((n.cells^2 - n.cells)/2)
  metric <- as.dist(metric.matrix)
  
  out <- list(
    gene.clusters = gene.clusters,
    points = points,
    metric = metric,
    seurat = seurat,
    assay = assay,
    cells = cells,
    npcs = npcs
  )
  names(out$gene.clusters) <- cluster.labels
  return(out)
}


#' @export
visualize_moduli <- function(moduli, clusters = NULL, title = NULL, n_neighbors = 15, seed = 123){
  set.seed(seed)
  umap.coords <- uwot::umap(moduli$metric, n_neighbors = n_neighbors)
  if(is.null(clusters)){
    plot(umap.coords, main = title)
  }else{
    has.cluster <- sapply(moduli$points, function(pt) all(clusters %in% pt))
    color <- character(length(has.cluster))
    color[has.cluster] <- "red"
    color[!has.cluster] <- "black"
    plot(umap.coords, col = color, main = title)
  }
}

#' @export
restrict_moduli <- function(moduli, clusters){
  res <- list()
  res$gene.clusters <- moduli$gene.clusters[clusters]
  
  keep <- sapply(moduli$points, function(pt)all( pt %in% clusters))
  res$points <- moduli$points[keep]
  res$metric <- as.dist(as.matrix(moduli$metric)[keep, keep])
  res$cells <- moduli$cells
  
  res$seurat <- moduli$seurat
  res$assay <- moduli$assay
  res$npcs <- moduli$npcs
  
  return(res)
}



#' @export
get_freq_table <- function(moduli, partition){
 cluster.labels <- names(moduli$gene.clusters)
 partition.labels <- sort(unique(partition))
 freq.table <- matrix(data = 0, nrow = length(partition.labels),
                     ncol = length(cluster.labels))
 row.names(freq.table) <- partition.labels
 colnames(freq.table) <- cluster.labels
 
  for(i in seq_along(partition.labels)){
    part <- which(partition == partition.labels[i])
    tab <- table(unlist(moduli$points[part] ))
    freq.table[i, cluster.labels %in% names(tab)] <- tab/length(part)
  }
 return(freq.table)
}

#' @export
partition_moduli <- function(moduli, k = 5, verbose = T){
  knn_moduli <- create_knn(moduli$metric, 5)
  jaccard <- igraph::similarity(knn_moduli, method = "jaccard", loops = T)
  jaccard.graph <- igraph::graph_from_adjacency_matrix(jaccard, mode = "undirected", weighted = T, diag = F)

  clusters <- igraph::cluster_louvain(jaccard.graph)
  
  if(verbose){
    print(igraph::modularity(jaccard.graph, clusters$membership, weights = igraph::E(jaccard.graph)$weight))
    print(table(clusters$membership))
  }
  return(clusters$membership)
}


#' @export
find_genes <- function(moduli, genes, ignore.case = T){
  out <- list()
  for(i in seq_along(moduli$gene.clusters)){
    if(ignore.case){
      found <- tolower(moduli$gene.clusters[[i]]) %in% tolower(genes)
    } else
    {
      found <- moduli$gene.clusters[[i]] %in% genes
    }
    out[[i]] <- moduli$gene.clusters[[i]][found]
  }
  names(out) <- names(moduli$gene.clusters)
  out <- out[sapply(out, function(x) length(x) > 0)]
  return(out)
}


#' @export
get_idents <- function(moduli, point){
  seurat <- moduli$seurat
  features <- unique(unlist(moduli$gene.clusters[point]))
  approx  <- (length(features) >= moduli$npcs*3)
  seurat <- RunPCA(seurat, assay = moduli$assay, features = features,
                    npcs = moduli$npcs, approx = approx, verbose = F)
  
  seurat <- FindNeighbors(seurat, verbose = F)
  seurat <- FindClusters(seurat, verbose = F)
  return(Idents(seurat))
}

#' @export
integrate_idents <- function(moduli){
  n.pts <- length(moduli$points)
  all_idents <- sapply(1:n.pts, function(i) get_idents(moduli, i) ) #each column is a membership
  n.cells <- nrow(all_idents)
  out <- matrix(data = 0, nrow = n.cells, ncol = n.cells)
  for(i in 1:n.pts){
    for(k in unique(all_idents[,i])){
      coords <- which(all_idents[,i] == k)
      for(j in coords){
        out[j, coords] <- out[j, coords] + 1
      }
    }
  }
  out <- out/n.pts
  row.names(out) <- row.names(all_idents)
  colnames(out) <- row.names(all_idents)
  return(out)
}




