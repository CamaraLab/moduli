

#' Clusters moduli space using the Leiden algorithm
#' 
#' Create clusters of analysis by applying \link[leiden]{leiden} to the snn graph
#' of the moduli
#' 
#' @param moduli A moduli object with a graph in the \code{snn.graph} slot
#' @param resolution_paramenter Resolution passed to the Leiden algorithm (see \link[leiden]{leiden}), default
#' value is 1
#' @param partition_type Type of partition to use (see \link[leiden]{leiden}), default value is
#' "RBConfigurationVertexPartition"
#' @param  enrich Whether to enrich the analysis clusters with differentially expressed gene clusters
#' by calling \link{enrich_analyis_clusters}. Default value is TRUE
#' @param thld Significance threshold for gene clusters, only used if \code{enrich} is TRUE
#' @param seed Random seed, default value is 123
#' @param ... Additional arguments passed to \link[leiden]{leiden}
#' 
#' 
#' @return A moduli object with analysis clusters saved in the \code{analysis.clusters} slot.
#' 
#' @export
cluster_moduli_space <- function(moduli, resolution_parameter = 1,
                                 partition_type = "RBConfigurationVertexPartition",
                                 enrich = T, thld = 0.05, seed = 123, ...){
  if(is.null(moduli$snn.graph)){
    stop("Error: snn graph required, run get_snn first")
  }
  
  membership <- leiden::leiden(
    moduli$snn.graph,
    partition_type = partition_type,
    resolution_parameter = resolution_parameter,
    seed = seed,
    ...
  )
  out <- moduli
  out$analysis.clusters <- data.frame(id = sort(unique(membership)))
  out$analysis.clusters$points <- lapply(out$analysis.clusters$id,
                                           function(p) out$points$id[membership == p])
  
  if(enrich) out <- enrich_analyis_clusters(out, thld)
  
  return(out)
}

#' Enriches analysis clusters with differentially expressed gene clusters
#' 
#' Applies Fisher's exact test to find differentially expressed gene clusters in
#' each analysis cluster.
#' 
#' @param moduli A moduli object with analysis clusters save to the \code{analysis.clusters} slot
#' @param thld Significance threshold
#' @return A moduli object with diffentially expressed gene clusters and associated p-values
#' saved in the \code{analysis.clusters$exp.gene.clusters} and \code{analysis.clusters$enrichment.p.vals}
#' respectively, ordered by of significance.
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
  out <- moduli
  out$analysis.clusters$exp.gene.clusters <- exp.gene.clusters
  out$enrichment.p.vals <- enrichment.p.vals
  return(out)
}


