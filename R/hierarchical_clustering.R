
#' Organizes the analysis clusters 
#' 
#' Clusters the moduli space hierarchically, starting by the analysis clusters
#' obtained by cluster_moduli_space.
#' 
#' @param moduli A moduli object with analysis clusters
#' @param method Agglomeration method used (see \link[stats]{hclust}). The default value
#' is \code{"averege"}.
#' @param thld Significance threshold for gene cluster enrichment
#' 
#' @return A moduli object with a dendrogram representation of the clusters in the
#' \code{dendro} slot and associated metadata in the \code{cluster.group.metadata}
#' slot (see \link{cluster_group_metadata})
#' 
#' @examples 
#' data("pbmc_small_moduli")
#' 
#' # clustering moduli space
#' pbmc_small_moduli <- get_snn(pbmc_small_moduli, 4)
#' pbmc_small_moduli <- cluster_moduli_space(pbmc_small_moduli)
#' 
#' pbmc_small_moduli <- organize_analysis_clusters(pbmc_small_moduli)
#' 
#' 
#' @export
organize_analysis_clusters <- function(moduli, method = "average", thld = 0.05){
  
  n.ac <- nrow(moduli$analysis.clusters)
  mean.metric <- matrix(nrow = n.ac, ncol = n.ac)
  metric.matrix <- as.matrix(moduli$metric)
  for(i in 2:n.ac){
    s1 <- moduli$points$id %in% moduli$analysis.clusters$points[[i]]
    for(j in 1:(i-1)){
      s2 <- moduli$points$id %in% moduli$analysis.clusters$points[[j]]
      mean.metric[i, j] <- mean(metric.matrix[s1,s2])
    }
  }
  rownames(mean.metric) <- moduli$analysis.clusters$id
  colnames(mean.metric) <- moduli$analysis.clusters$id
  hcls <- hclust(as.dist(mean.metric), method = method,
                 members = sapply(moduli$analysis.clusters$points, length))
  moduli$dendro <- as.dendrogram(hcls)
  
  moduli$cluster.group.metadata <- cluster_group_metadata(moduli, thld)
    
  return(moduli)
}


#' Retrieves analysis cluster groups
#' 
#' @param moduli A moduli object with a dendrogram representation
#' 
#' @return A data.frame with ids of analysis cluster groups, their heights, an
#' indicator if they are leafs (i.e. correspond to a single analysis cluster) and
#' analysis clusters associated to each group. Cluster groups are ordered
#' by depth first search.
#' 
#' @examples 
#' data("pbmc_small_moduli")
#' 
#' # clustering moduli space
#' pbmc_small_moduli <- get_snn(pbmc_small_moduli, 4)
#' pbmc_small_moduli <- cluster_moduli_space(pbmc_small_moduli)
#' 
#' pbmc_small_moduli <- organize_analysis_clusters(pbmc_small_moduli)
#' get_analysis_cluster_groups(pbmc_small_moduli)
#' 
#' @export
get_analysis_cluster_groups <- function(moduli){
  i <- 1
  n.nodes <- 2*nrow(moduli$analysis.clusters) - 1
  height <- numeric(n.nodes)
  leaf <- logical(n.nodes)
  clusters <- NULL
  
  dfs <- function(node){
    pos <- i
    i <<- i + 1
    height[pos] <<- attr(node, "height")
    if(is.null(attr(node, "leaf"))){
      leaf[pos] <<- F
      pl <- dfs(node[[1]])
      pr <- dfs(node[[2]])
      clusters[[pos]] <<- append(clusters[[pl]], clusters[[pr]])
    } else {
      leaf[pos] <<- T
      clusters[[pos]] <<- as.integer(attr(node, "label"))
    }
    return(pos)
  }
  dfs(moduli$dendro)
  out <- data.frame(id = 1:n.nodes, height = height, leaf = leaf, clusters = I(clusters))
  return(out)
}

#' Retrieves metadata associated to analysis cluster groups
#' 
#' @param moduli A moduli object with a dendrogram representation
#' @param thld Significance threshold
#' 
#' @return A data.frame with ids of analysis cluster groups, their heights, an
#' indicator if they are leafs (i.e. correspond to a single analysis cluster), number
#' of points in the moduli present in the analysis cluster group, differentially 
#' expressed and suppressed gene clusters with associated p-values (Fisher's exact test).
#' Cluster groups are ordered by depth first search.
#' 
#' @examples 
#' data("pbmc_small_moduli")
#' 
#' # clustering moduli space
#' pbmc_small_moduli <- get_snn(pbmc_small_moduli, 4)
#' pbmc_small_moduli <- cluster_moduli_space(pbmc_small_moduli)
#' 
#' pbmc_small_moduli <- organize_analysis_clusters(pbmc_small_moduli)
#' cluster_group_metadata(pbmc_small_moduli)
#' 
#' @export
cluster_group_metadata <- function(moduli, thld = 0.05){
  
  if(!is.null(moduli$cluster.group.metadata)) return(moduli$cluster.group.metadata)
  
  clst.grps <- get_analysis_cluster_groups(moduli)
  points <- lapply(
    clst.grps$clusters,
    function(cg) unlist(moduli$analysis.clusters$points[moduli$analysis.clusters$id %in% cg])
  )
  exp.data <- .find_expr(moduli, points, thld)
  out <- data.frame(id = clst.grps$id, height = clst.grps$height, leaf = clst.grps$leaf)
  out$n.pts <- sapply(points, length)
  out <- cbind(out, exp.data)
  return(out)
}


#' Plot a dendrogram of the analysis cluster groups
#' 
#' Create an interactive visualization of analysis cluster groups.
#' 
#' @param moduli A moduli object with a dendrogram representation
#' @param hang Whether to hang dendrogram leafs, the default value is FALSE. Setting to
#' TRUE can help distinguish leafs.
#' 
#' @return A plotly visualization. Hover text contains about analysis cluster groups. Each leaf
#' of the dendrogram represents an analysis cluster and is labeled with the respective id.
#' 
#' @examples
#' data("pbmc_small_moduli")
#' 
#' # clustering moduli space
#' pbmc_small_moduli <- get_snn(pbmc_small_moduli, 4)
#' pbmc_small_moduli <- cluster_moduli_space(pbmc_small_moduli)
#' 
#' pbmc_small_moduli <- organize_analysis_clusters(pbmc_small_moduli)
#' visualize_analysis_dendrogram(pbmc_small_moduli)
#' 
#' @export
visualize_analysis_dendrogram <- function(moduli, hang = F){
  metadata <- moduli$cluster.group.metadata
  
  txt <- paste0(
    "id: ", metadata$id,
    "<br>height: ", signif(metadata$height, 3),
    "<br>n.pts: ", metadata$n.pts,
    "<br>exp.gene.clst: ", sapply(metadata$exp.gene.clusters, function(x) paste0(x, collapse = ", ")), 
    "<br>sup.gene.clst: ", sapply(metadata$sup.gene.clusters, function(x) paste0(x, collapse = ", "))
  )
  
  if(hang){
    ggd <- dendextend::as.ggdend(dendextend::hang.dendrogram(moduli$dendro))
  } else {
    ggd <- dendextend::as.ggdend(moduli$dendro)
  }
  plt <- plotly::add_segments(
    plotly::plot_ly(),
    data = ggd$segments, 
    x = ~x,
    y = ~y,
    xend = ~xend,
    yend = ~yend, 
    hoverinfo = "none",
    color = I("black"),
    showlegend = F
  )
  plt <- plotly::add_markers(
    plt,
    data = ggd$nodes,
    x = ~x,
    y = ~y,
    color = I("black"),
    text = txt,
    hoverinfo = "text"
  )
  plt <- plotly::add_annotations(
    plt,
    x = ggd$labels$x,
    y = ggd$labels$y - max(metadata$height)/100,
    text = ggd$labels$label,
    yanchor = "top",
    showarrow = F,
    textangle = 90
  )
  plt <- plotly::layout(
    plt,
    title = "Dendrogram of analysis",
    xaxis = list(visible = FALSE),
    yaxis = list(title = "height", zeroline = FALSE)
  )
  return(plt)
}
