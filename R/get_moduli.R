#' @export
get_moduli <- function(seuratObject, membership = NULL, gene.clusters = NULL,
                       assay = DefaultAssay(seuratObject), npcs = 30, weight.by.var = T,
                       points = NULL, n.cores = 1, n.pts = 0, n.cells = 0, stride = 10000,
                       seed = 123, filebacked = T, persist = F, verbose = T){
  # gene.clusters : clusters of genes as a list of vectors of gene names
  if(sum(c(is.null(membership), is.null(gene.clusters))) != 1){
    stop("Error: provide exatcly one of the arguments membership or gene.clusters")
  }
  
  if(!is.null(membership)){
    gene.clusters <- list()
    cluster.ids <- sort(unique(membership))
    gene.clusters <- list()
    for(k in cluster.ids){
      gene.clusters[[k]] <- names(membership[membership == k])
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



