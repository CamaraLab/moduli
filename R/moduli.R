library(Seurat)
library(parallel)
library(bigmemory)
library(doParallel)
library(foreach)

source("clustering.R")
Rcpp::sourceCpp("compute_embedding_distances.cpp")


get_moduli <- function(seurat, membership = NULL, gene.clusters = NULL,
                       assay = DefaultAssay(seurat), npcs = 30, points = NULL,
                       pow = 1, n.cores = 1, downsample = 1, seed = 123, 
                       save.intermediary = F){
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
      comb <- comb[runif(length(comb)) <= downsample]
      points <- append(points, comb)
    }
  }
  
  is.large.enougth <- vapply(
      points,
      function(set) length(unique(unlist(gene.clusters[set]))) >= npcs,
      FUN.VALUE = rep_len(T, length(comb))
  )
  points <- points[is.large.enougth]
  n.pts <- length(points)
  
  embeddings <- filebacked.big.matrix(
    nrow = length(Cells(seurat))*npcs,
    ncol = length(points),
    backingfile = "embeddings.bin",
    descriptorfile = "embeddings.desc"
  )
  
  emb.descr <- describe(embeddings)
  
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  clusterEvalQ( cl = cl,{
    library(bigmemory)
    library(Seurat)
  })
  
  message(paste(format(Sys.time(), "%a %b %d %X"),"Building embeddings"))
  
  foreach(i = seq_along(points))%dopar%{
    
    emb <- attach.big.matrix(emb.descr)
    feature.subset <- unique(unlist(gene.clusters[points[[i]]]))
    approx  <- (length(feature.subset) >= npcs*3)
    
    emb[,i] <- c(
      Embeddings(RunPCA(seurat, assay = assay, features = feature.subset, npcs = npcs,
                          verbose = F, approx = approx))
    )
  }
  message(paste(format(Sys.time(), "%a %b %d %X"),"Finished building embeddings"))
  
  
  message(paste(format(Sys.time(), "%a %b %d %X"), "Computing distances"))
  
  answ <- NULL
  for(i in seq_along(points)){
    emb <- attach.big.matrix(emb.descr)
    dd <- dist_Cpp(emb@address, i, npcs, pow)
    answ <- cbind( answ, dd)
  }
  message("asdadsasdasd")
  
  metric.matrix <- foreach(i = seq_along(points),
                           .noexport = c("seurat", "gene.clusters", "points"),
                           .export = "dist_Cpp",
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
    npcs = npcs
  )
  names(out$gene.clusters) <- cluster.labels
  return(out)
}


visualize_moduli <- function(moduli, clusters = NULL, title = NULL, n_neighbors = 15){
  set.seed(123)
  umap.coords <- umap(as.dist(moduli$metric.matrix), n_neighbors = n_neighbors)
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

partition_moduli <- function(moduli, k = 5, verbose = T){
  knn_moduli <- create_knn(moduli$metric.matrix, 5)
  jaccard <- similarity(knn_moduli, method = "jaccard", loops = T)
  jaccard.graph <- graph_from_adjacency_matrix(jaccard, mode = "undirected", weighted = T, diag = F)

  clusters <- igraph::cluster_louvain(jaccard.graph)
  
  if(verbose){
    print(modularity(jaccard.graph, clusters$membership, weights = E(jaccard.graph)$weight))
    print(table(clusters$membership))
  }
  return(clusters$membership)
}


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




