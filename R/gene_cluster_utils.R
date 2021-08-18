

#' @export
gene_medioid_clustering <- function(seuratObject, n.clusters, slot = "scale.data", 
                                    features = VariableFeatures(seuratObject)){
  data <- FetchData(seuratObject, vars = features, slot = slot)
  dist.mat <- 1 - cor(data)
  part <- cluster::pam(dist.mat, n.clusters, diss = T)
  return(part)
}

# Delete?
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

#' Annotates genes clusters using g:OSt
#' 
#' Enriches gene clusters with most significant terms according to g:OSt by
#' using \link[gprofiler2]{gost}.
#' 
#' @param moduli A moduli object
#' @param organism Organism name, constructed by concatenating the fist letter of the
#' name with the family name
#' @param n.tems Number of terms to save
#' @param user_threshold p-value threshold for significance
#' @param ... Additional arguments passed to \link[gprofiler2]{gost}
#' 
#' @return A moduli object with the significant term names and associated p-values
#' saved to moduli$gene.clusters$term.names and moduli$gene.clusters$term.p.values
#' respectively, ordered by significance.
#' 
#' @export
annotate_gene_clusters <- function(moduli, organism, n.terms = 3, user_threshold = 0.01, ...){
  result <- gprofiler2::gost(moduli$gene.clusters$genes, organism = organism,
                             user_threshold = user_threshold, ...)$result
  p.vals <- NULL
  term.names <- NULL
  for(i in 1:nrow(moduli$gene.clusters)){
    name <- paste0("query_", i)
    p <- result$p_value[result$query == name]
    tn <- result$term_name[result$query == name]
    pos <- order(p)[1:min(n.terms, length(p))]
    p.vals[[i]] <- p[pos]
    term.names[[i]] <- tn[pos]
  }
  out <- moduli
  out$gene.clusters$term.names <- term.names
  out$gene.clusters$term.p.vals <- p.vals
  return(out)
}


#' Enriches gene clusters with laplacian scores
#' 
#' Computes for each gene cluster the graph laplacian score of the function that indicates the
#' presence of the gene cluster in the point, using the snn graph of the moduli.
#' 
#' @param moduli A moduli object with a graph in the \code{snn.graph} slot
#' @return A moduli object with laplacian scores and the rank of the gene cluster 
#' laplacian score saved in \code{gene.clusters$laplacian.score} and \code{gene.clusters$rank}
#' slots, respectively. The rank goes from most localized (smallest laplacian score) to least
#' localized (largest laplacian score).
#' 
#' @export
score_gene_clusters <- function(moduli){
  if(!requireNamespace("RayleighSelection", quietly = TRUE)){
    stop("Error: package CamaraLab/RayleighSelection required")
  }
  if(is.null(moduli$snn.graph)){
    stop("Error: snn graph required, run get_snn first")
  }
  f <- matrix(data = 0, nrow = nrow(moduli$gene.clusters), ncol = nrow(moduli$points))
  for(i in 1:nrow(moduli$points)){
    f[moduli$gene.clusters$id %in% moduli$points$clusters[[i]], i] <- 1
  }
  
  cplx <- RayleighSelection::graph_to_complex(
    igraph::as_adjacency_matrix(moduli$snn.graph, attr="weight"),
    clique = F
  )
  scores <- RayleighSelection::rayleigh_selection(cplx, f, num_perms = 1, one_forms = F)$R0
  out <- moduli
  out$gene.clusters$laplacian.score <- scores
  out$gene.clusters$rank <- rank(scores)
  return(out)
}


#' Retrieves gene cluster metadata
#' 
#' @param moduli A moduli object
#' @return A data.frame
#' 
#' @export
gene_cluster_metadata <- function(moduli){
  out <- data.frame(id = moduli$gene.clusters$id)
  out$size <- sapply(moduli$gene.clusters$genes, length)
  out$term.names <- moduli$gene.clusters$term.names
  out$term.p.vals <- moduli$gene.clusters$term.p.vals
  out$laplacian.score <- moduli$gene.clusters$laplacian.score
  out$rank <- moduli$gene.clusters$rank
  return(out)
}

