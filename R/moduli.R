

#' @export
update_moduli <- function(moduli){
  out <- moduli
  out$gene.clusters <- data.frame(id = seq_along(moduli$gene.clusters))
  out$gene.clusters$genes <- moduli$gene.clusters
  
  out$points <- data.frame(id = seq_along(moduli$points))
  out$points$clusters <- moduli$points
  return(out)
}



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
            in_complement = !any(gene.clusters %in% pt.clusters)
    )
  }
  return(moduli$points$id[sapply(moduli$points$clusters, f)])
}

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

#' @export
find_genes <- function(moduli, genes, ignore.case = T){
  out <- list()
  for(i in seq_along(moduli$gene.clusters$genes)){
    if(ignore.case){
      found <- tolower(moduli$gene.clusters$genes[[i]]) %in% tolower(genes)
    } else {
      found <- moduli$gene.clusters$genes[[i]] %in% genes
    }
    out[[i]] <- moduli$gene.clusters$genes[[i]][found]
  }
  names(out) <- moduli$gene.clusters$id
  out <- out[sapply(out, function(x) length(x) > 0)]
  return(out)
}


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




