

#' Clusters moduli space using the Louvain or the Leiden algorithm
#' 
#' Create clusters of analysis (i.e. clusters of points in the moduli) by applying 
#' \link[igraph]{cluster_louvain} or \link[leiden]{leiden} to the snn graph of the moduli.
#' 
#' @param moduli A moduli object with a graph in the \code{snn.graph} slot
#' @param algorithm Which algorithm to use, the options are "Louvain" and "Leiden", default value is
#' "Louvain"
#' @param ... Additional arguments passed to \link[leiden]{leiden}
#' @param  enrich Whether to enrich the analysis clusters with differentially expressed gene clusters
#' by calling \link{enrich_analyis_clusters}. Default value is TRUE
#' @param thld Significance threshold for gene clusters, only used if \code{enrich} is TRUE
#' @param seed Random seed, default value is 123
#' 
#' @return A moduli object with analysis clusters saved in the \code{analysis.clusters} slot.
#' 
#' @examples
#' data("pbmc_small_moduli")
#' pbmc_small_moduli <- get_snn(pbmc_small_moduli, 4)
#' pbmc_small_moduli <- cluster_moduli_space(pbmc_small_moduli)
#' 
#' @export
cluster_moduli_space <- function(moduli, algoritm = c("Louvain", "Leiden"), ...,
                                 enrich = T, thld = 0.05, seed = 123){
  if(is.null(moduli$snn.graph)){
    stop("Error: snn graph required, run get_snn first")
  }
  algoritm <- match.arg(algoritm)
  if(algoritm == "Louvain"){
    set.seed(seed)
    communities <- igraph::cluster_louvain(moduli$snn.graph)
    membership <- communities$membership
  } else if(algoritm == "Leiden"){
    if(!requireNamespace("leiden", quietly = TRUE)){
      stop("Error: package TomKellyGenetics/leiden required")
    }
    membership <- leiden::leiden( moduli$snn.graph, seed = seed, ...)
  }
  
  moduli <- moduli
  moduli$analysis.clusters <- data.frame(id = sort(unique(membership)))
  moduli$analysis.clusters$points <- lapply(moduli$analysis.clusters$id,
                                           function(p) moduli$points$id[membership == p])
  
  if(enrich) moduli <- enrich_analysis_clusters(moduli, thld)
  
  return(moduli)
}

#' Enriches analysis clusters with differentially expressed gene clusters
#' 
#' Applies Fisher's exact test to find differentially expressed gene clusters in
#' each cluster of analysis.
#' 
#' @param moduli A moduli object with analysis clusters save to the \code{analysis.clusters} slot
#' @param thld Significance threshold
#' @return A moduli object with diffentially expressed gene clusters and associated p-values
#' saved in the \code{analysis.clusters$exp.gene.clusters} and \code{analysis.clusters$enrichment.p.vals}
#' respectively, ordered by of significance.
#' 
#' @examples
#' data("pbmc_small_moduli")
#' pbmc_small_moduli <- get_snn(pbmc_small_moduli, 4)
#' pbmc_small_moduli <- cluster_moduli_space(pbmc_small_moduli, enrich = F)
#' pbmc_small_moduli <- enrich_analysis_clusters(pbmc_small_moduli)
#' 
#' 
#' @export
enrich_analysis_clusters <- function(moduli, thld = 0.05){
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
  moduli$analysis.clusters$exp.gene.clusters <- exp.gene.clusters
  moduli$analysis.clusters$enrichment.p.vals <- enrichment.p.vals
  return(moduli)
}

#' Retrieved analysis cluster metadata
#' 
#' @param moduli A moduli object with analysis clusters save to the \code{analysis.clusters} slot
#' @returns A data.frame with ids, sizes, and, if available, expressed gene clusters and associated
#' p-values of each analysis cluster.
#' 
#' @examples
#' data("pbmc_small_moduli")
#' pbmc_small_moduli <- get_snn(pbmc_small_moduli, 4)
#' pbmc_small_moduli <- cluster_moduli_space(pbmc_small_moduli)
#' analysis_cluster_metadata(pbmc_small_moduli)
#' 
#' @export
analysis_cluster_metadata <- function(moduli){
  out <- data.frame(id = moduli$analysis.clusters$id)
  out$size <- sapply(moduli$analysis.clusters$points, length)
  out$exp.gene.clusters <- moduli$analysis.clusters$exp.gene.clusters
  out$enrichment.p.vals <- moduli$analysis.clusters$enrichment.p.vals
  return(out)
}




