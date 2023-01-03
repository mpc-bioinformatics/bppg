

#' Generate bipartite peptide-protein graphs from a list of digested proteins via an edgelist
#'
#' @param edgelist Output from generate_edgelist (edgelist)
#'
#' @return List of subgraphs as igraph objects.
#' @export
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "2020_01_31_proteome_S_cerevisae.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' digested_proteins <- digest_fasta(fasta)
#' edgelist <- generate_edgelist(digested_proteins)
#' res <- generate_graphs_from_edgelist(edgelist)
#'
generate_graphs_from_edgelist <- function(edgelist) {

  #generate graph from edge matrix
  G <- igraph::graph_from_edgelist(edgelist, directed = FALSE)

  #assign vertex types to proteins and peptides for the graph to be bipartite
  igraph::V(G)[igraph::V(G)$name %in% edgelist[,1]]$type <- TRUE
  igraph::V(G)[igraph::V(G)$name %in% edgelist[,2]]$type <- FALSE
  ### TODO: export G

  #decompose graph into connected components
  subgraphs <- igraph::decompose(G)
  return(subgraphs)

}




