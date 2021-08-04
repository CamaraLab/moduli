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
                       assay = DefaultAssay(seurat), npcs = 30, points = NULL,
                       pow = 1, n.cores = 1, size = 0, ncells = 0, seed = 123){
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

  is.large.enougth <- vapply(
      points,
      function(set) length(unique(unlist(gene.clusters[set]))) >= npcs,
      FUN.VALUE = rep_len(T, length(comb))
  )
  points <- points[is.large.enougth]
  if( size > 0 && length(points) > size) sample(points, size)
  n.pts <- length(points)
  
  if(ncells > 0 && length(Cells(seurat)) > ncells){
    cells <- sample(Cells(seurat), ncells) 
  } else {
    cells <- Cells(seurat)
  }
  scaled.data <- GetAssayData(seurat[[assay]], slot = "scale.data" )[,cells]
  
  embeddings <- bigmemory::filebacked.big.matrix(
    nrow = length(cells)*npcs,
    ncol = length(points),
    backingfile = "embeddings.bin",
    descriptorfile = "embeddings.desc"
  )
  
  emb.descr <- bigmemory::describe(embeddings)
  
  cl <- parallel::makeCluster(n.cores)
  on.exit(parallel::stopCluster(cl))
  registerDoParallel(cl)
  
  
  #clusterEvalQ( cl = cl,{
  #  library(bigmemory)
  #  library(Seurat)
  #})
  
  
  message(paste(format(Sys.time(), "%a %b %d %X"),"Building embeddings"))
  foreach(i = seq_along(points),
          .noexport = c("seurat"),
          .packages = c("bigmemory", "irlba"))%dopar%{
    feature.subset <- unique(unlist(gene.clusters[points[[i]]]))
    if (length(feature.subset) > npcs*2){
      pca.results <- irlba(t(scaled.data[feature.subset,]), nv = npcs)
      cell.embeddings <- pca.results$u %*% diag(pca.results$d) #weighting by variance
    } else {
      pca.results <- prcomp(t(scaled.data[feature.subset,]), rank. = npcs)
      cell.embeddings <- pca.results$x
    }
    
    emb <- attach.big.matrix(emb.descr)
    emb[,i] <- c(cell.embeddings)
  }
  #message(paste(format(Sys.time(), "%a %b %d %X"),"Finished building embeddings"))
  
  
  message(paste(format(Sys.time(), "%a %b %d %X"), "Computing distances"))
  metric.matrix <- foreach(i = seq_along(points),
                           .noexport = c("seurat", "gene.clusters", "points"),
                           .packages = c("moduli"),
                           .combine =  cbind)%dopar%{
    
    emb <- attach.big.matrix(emb.descr)
    dist_Cpp(emb@address, i, npcs, pow)
    #d <- numeric(n.pts)
    #if(i < n.pts){
    #  m1 <- dist(matrix(data = emb[,i], ncol = npcs))
    #  m1 <- m1*(length(m1)/sum(m1))
    #
    #  for(j in (i+ 1):length(points)){
    #    m2 <- dist(matrix(data = emb[,j], ncol = npcs))
    #    m2 <- m2*(length(m2)/sum(m2))
    #    d[j] <- (sum(abs(m1 - m2)^pow)/length(m2))^(1/pow)
    #  }
    #}
    #d
  }
  message(paste(format(Sys.time(), "%a %b %d %X"), "Finished computing distances"))

  metric.matrix <- metric.matrix + t(metric.matrix)
  
  out <- list(
    gene.clusters = gene.clusters,
    points = points,
    metric.matrix = metric.matrix,
    seurat = seurat,
    assay = assay,
    cells = cells,
    npcs = npcs
  )
  names(out$gene.clusters) <- cluster.labels
  return(out)
}


#' @export
visualize_moduli <- function(moduli, clusters = NULL, title = NULL, n_neighbors = 15){
  set.seed(123)
  umap.coords <- uwot::umap(as.dist(moduli$metric.matrix), n_neighbors = n_neighbors)
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
  res$metric.matrix <- moduli$metric.matrix[keep, keep]
  
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
  knn_moduli <- create_knn(moduli$metric.matrix, 5)
  jaccard <- igraph::similarity(knn_moduli, method = "jaccard", loops = T)
  jaccard.graph <- igraph::graph_from_adjacency_matrix(jaccard, mode = "undirected", weighted = T, diag = F)

  clusters <- igraph::cluster_louvain(jaccard.graph)
  
  if(verbose){
    print(igraph::modularity(jaccard.graph, clusters$membership, weights = E(jaccard.graph)$weight))
    print(table(clusters$membership))
  }
  return(clusters$membership)
}


#' @export
find_genes <- function(moduli, genes, ignore.case = T){
  out <- list()
  for(i in seq_along(moduli$gene.clusters)){
    if(ignore.case){
      genes.found <- genes[tolower(genes) %in% tolower(moduli$gene.clusters[[i]])]
    }else
    {
      genes.found <- genes[genes %in% moduli$gene.clusters[[i]]]
    }
    for(j in seq_along(genes.found)){
      out[[genes.found[j]]] <- append(out[[genes.found[j]]], names(moduli$gene.clusters)[i])
    }
  }
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




