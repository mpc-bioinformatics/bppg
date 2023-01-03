


#' Generate graphs from a FASTA file
#'
#' @param fasta fasta file, already read into R by seqinr::read.fasta
#' @param collapse_protein_nodes collapse protein nodes?
#' @param ... additional arguments to bppg::digest_fasta()

#'
#' @return subgraphs (i.e. connected components) from the graph generated from the FASTA file.
#' @export
#'
#' @examples
#' #' library(seqinr)
#' file <- system.file("extdata", "2020_01_31_proteome_S_cerevisae.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' graphs <- bppg::generate_graphs_from_FASTA(fasta)
#'
generate_graphs_from_FASTA <- function(fasta, collapse_protein_nodes = TRUE, ...) {


  digested_proteins <- bppg::digest_fasta(fasta)#, ...)
  edgelist <- bppg::generate_edgelist(digested_proteins)
  graphs <- bppg::generate_graphs_from_edgelist(edgelist)

  if (collapse_protein_nodes) {
    graphs <- bppg::collapse_protein_nodes(graphs, sparse = TRUE, fast = FALSE, matrix = FALSE)
  }

  return(graphs)

}
