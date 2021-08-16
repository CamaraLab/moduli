#' @export
visualize_moduli_space <- function(moduli, n_neighbors = 15, mark.points = NULL,
                                   color.clusters = NULL, title = "UMAP plot of moduli space",
                                   seed = 123){
  set.seed(seed)
  if(!is.null(mark.points)){
    mark.legend <- c("marked", "unmarked")
    mark.factors <- ifelse(moduli$points$id %in% mark.points, mark.legend[1], mark.legend[2])
    mark.factors <- factor(mark.factors, levels = mark.legend)
  }
  
  if(!is.null(moduli$analysis.clusters) &&  is.null(color.clusters)){
    color.clusters <- sort(moduli$analysis.clusters$id)
    if(length(color.clusters) > 12) color.clusters <- color.clusters[1:12]
  }
 
  if(!is.null(color.clusters) && length(color.clusters) > 12){
    warning("Only the first 12 elements of color.clusters will be colored")
    color.clusters <- color.clusters[1:12]
  }
  
  
  umap.coords <- uwot::umap(moduli$metric, n_neighbors = n_neighbors)
  data = data.frame(
    UMAP_1 = umap.coords[,1],
    UMAP_2 = umap.coords[,2]
  )
  if(!is.null(mark.points)){
    data$mark.factors <- mark.factors
  }
  if(!is.null(moduli$analysis.clusters)){
    partition.factors <- character(nrow(moduli$points))
    point.memberships <- integer(nrow(moduli$points))
    for(i in 1:nrow(moduli$analysis.clusters)){
      point.memberships[moduli$points$id %in% moduli$analysis.clusters$points[[i]]] <- moduli$analysis.clusters$id[i] 
    }
    partition.factors[!(point.memberships %in% color.clusters)] <- "others"
    for(p in color.clusters){
      partition.factors[point.memberships == p] <- paste("analysis cluster", p)
    }
    data$partition <- factor(partition.factors, levels = c(paste("analysis cluster", color.clusters), "others") )
  }
  
  # compiling hoover text  
  gene.cluster.info <- sapply(moduli$points$clusters, function(pt) paste(sort(pt), collapse = " "))
  # marking differerentially expreessed clusters with *
  if(!is.null(moduli$analysis.clusters$exp.gene.clusters)){
    for(i in 1:nrow(moduli$analysis.clusters)){
      members <- moduli$points$id %in% moduli$analysis.clusters$points[[i]]
      splt.labels <- strsplit(gene.cluster.info[members], split = " ")
      for(gcl in moduli$analysis.clusters$exp.gene.clusters[[i]]){
        exp <- paste0("^", gcl, "$")
        rep <- paste0(gcl,"*")
        splt.labels <- lapply(splt.labels, function(x) sub(exp, rep, x))
      }
      gene.cluster.info[members] <- sapply(splt.labels, paste, collapse = " ")
    }
  }
 
  txt <- paste("point id:", moduli$points$id, "<br>gene clusters:", gene.cluster.info)
  if(!is.null(moduli$analysis.clusters)) txt <- paste0(txt,"<br>analysis cluster: ", point.memberships)
  
  # colors
  if(length(color.clusters) >= 3){
    cols <- c(RColorBrewer::brewer.pal(n = max(length(color.clusters)), name = "Paired"), "black")
  }
  if(length(color.clusters) == 2){
    cols <- c("red", "blue" , "black")
  }
  if(length(color.clusters) == 1){
    cols <- c("red", "black")
  }
  if(length(color.clusters) == 0){
    cols <- I("black")
  }
  
  
  if(is.null(mark.points) && is.null(moduli$analysis.clusters)){
    plt <- plotly::plot_ly(
      data = data,
      x = ~UMAP_1,
      y = ~UMAP_2,
      type = "scatter",
      mode = "markers",
      color = I("black"),
      text = ~txt,
      hoverinfo = "text"
    )
  }
  
  if(!is.null(mark.points) && is.null(moduli$analysis.clusters)){
    plt <- plotly::plot_ly(
      data = data,
      x = ~UMAP_1,
      y = ~UMAP_2,
      type = "scatter",
      mode = "markers",
      color = I("black"),
      symbol = ~mark.factors,
      symbols = c("x", "o"),
      text = ~txt,
      hoverinfo = "text"
    )
  }
  
  if(is.null(mark.points) && !is.null(moduli$analysis.clusters)){
    plt <- plotly::plot_ly(
      data = data,
      x = ~UMAP_1,
      y = ~UMAP_2,
      type = "scatter",
      mode = "markers",
      color = ~partition,
      colors = cols,
      text = ~txt,
      hoverinfo = "text"
    )
  }
  
  if(!is.null(mark.points) && !is.null(moduli$analysis.clusters)){
    plt <- plotly::plot_ly(
      data = data,
      x = ~UMAP_1,
      y = ~UMAP_2,
      type = "scatter",
      mode = "markers",
      color = ~partition,
      colors = cols,
      symbol = ~mark.factors,
      symbols = c("x", "o"),
      text = ~txt,
      hoverinfo = "text"
    )
  }
 
  plt <- plotly::layout(plt, title = title)
  return(plt)
}

#' @export
visualize_gene_space <- function(moduli, n_neighbors = 15, color.clusters = NULL, mark.genes = NULL,
                                 ignore.case = T,
                                 metric = NULL, title = "UMAP plot of gene space",
                                 seed = 123){

  if(is.null(color.clusters)){
    color.clusters <- sort(moduli$gene.clusters$id)
    if(length(color.clusters) > 12) color.clusters <- color.clusters[1:12]
  }
  if(!is.null(color.clusters) && length(color.clusters) > 12){
    warning("Only the first 12 elements of color.clusters will be colored")
    color.clusters <- color.clusters[1:12]
  }
  
  if(is.null(metric)){
    gene.exp <- t(FetchData(moduli$seurat[[moduli$assay]],
                            vars = unique(unlist(moduli$gene.clusters$genes)),
                            slot = "scale.data"))
    gene.names <- rownames(gene.exp)
    umap.coords <- uwot::umap(gene.exp, n_neighbors = n_neighbors, metric = "correlation")
  }else{
    gene.names <- rownames(as.matrix(gene.exp))
    umap.coords <- uwot::umap(metric, n_neighbors = n_neighbors)
  }
  
  set.seed(seed)
  data = data.frame(
    UMAP_1 = umap.coords[,1],
    UMAP_2 = umap.coords[,2]
  )
  
  if(!is.null(mark.genes)){
    mark.legend <- c("marked", "unmarked")
    if(ignore.case){
      mark.factors <- ifelse(tolower(gene.names) %in% tolower(mark.genes),
                             mark.legend[1], mark.legend[2])
    } else {
      mark.factors <- ifelse(gene.names %in% mark.genes, mark.legend[1], mark.legend[2])
    }
    data$mark.factors <- factor(mark.factors, levels = mark.legend)
  }
  
  # compiling gene cluster information
  membership <- integer(length(gene.names))
  for(i in 1:nrow(moduli$gene.clusters)){
    membership[gene.names %in% moduli$gene.clusters$genes[[i]]] <- moduli$gene.clusters$id[i]
  }
  partition.factors <- character(length(gene.names))
  partition.factors[!(membership %in% color.clusters)] <- "others"
  for(gcl in color.clusters){
    partition.factors[membership == gcl] <- paste("gene cluster", gcl)
  }
  data$partition <- factor(partition.factors, levels = c(paste("gene cluster", color.clusters), "others"))
  
  
  txt <- paste0(gene.names, "<br>gene cluster:", membership)
  
  # colors
  if(length(color.clusters) >= 3){
    cols <- c(RColorBrewer::brewer.pal(n = max(length(color.clusters)), name = "Paired"), "black")
  }
  if(length(color.clusters) == 2){
    cols <- c("red", "blue" , "black")
  }
  if(length(color.clusters) == 1){
    cols <- c("red", "black")
  }
  if(length(color.clusters) == 0){
    cols <- I("black")
  }
  
  
  if(is.null(mark.genes)){
    plt <- plotly::plot_ly(
      data = data,
      x = ~UMAP_1,
      y = ~UMAP_2,
      type = "scatter",
      mode = "markers",
      color = ~partition,
      colors = cols,
      text = ~txt,
      hoverinfo = "text"
    )
  } else {
    plt <- plotly::plot_ly(
      data = data,
      x = ~UMAP_1,
      y = ~UMAP_2,
      type = "scatter",
      mode = "markers",
      color = ~partition,
      colors = cols,
      symbol = ~mark.factors,
      symbols = c("x", "o"),
      text = ~txt,
      hoverinfo = "text"
    )
  }
  
  plt <- plotly::layout(plt, title = title)
  return(plt)
  
}

