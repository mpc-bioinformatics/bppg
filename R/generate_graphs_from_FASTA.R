### TODO: Matrix package wird auf jeden Fall benötigt

#' Generate graphs from a FASTA file
#'
#' @param fasta fasta file, already read into R by seqinr::read.fasta
#' @param collapse_protein_nodes collapse protein nodes?
#' @param collapse_peptide_nodes collapse peptide nodes?
#' @param result_path path whereresults are saved. If NULL, results are not saved
#' @param suffix suffix for saving results
#' @param save_intermediate Save intermediate results?
#' @param ... additional arguments to bppg::digest_fasta()

#'
#' @return subgraphs (i.e. connected components) from the graph generated from the FASTA file.
#' @export
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' graphs <- bppg::generate_graphs_from_FASTA(fasta)
#'
generate_graphs_from_FASTA <- function(fasta, collapse_protein_nodes = TRUE,
                                       collapse_peptide_nodes = TRUE,
                                       result_path = NULL,
                                       suffix = NULL, save_intermediate = FALSE,
                                       prot_origin = NULL,
                                       ...) {

  message("Digesting FASTA file...")
  digested_proteins <- bppg::digest_fasta(fasta, ...)#, ...)
  message("Generating edgelist ...")
  edgelist <- bppg::generate_edgelist(digested_proteins, prot_origin = prot_origin)
  if (save_intermediate) {
    message("Saving edgelist ...")
    utils::write.table(edgelist, sep = "\t", row.names = FALSE,
                                           file = paste0(result_path, "edgelist_", suffix, ".txt"))
  }


  if (collapse_protein_nodes | collapse_peptide_nodes) {
    message("Collapsing nodes ...")
    edgelist_coll <- bppg::collapse_edgelist(edgelist,
                                             collapse_protein_nodes = collapse_protein_nodes,
                                             collapse_peptide_nodes = collapse_peptide_nodes)
  }

  if(collapse_protein_nodes & collapse_peptide_nodes) suffix2 <- "collprotpept_"
  if(collapse_peptide_nodes & !collapse_protein_nodes) suffix2 <- "collpept_"
  if(collapse_protein_nodes & !collapse_peptide_nodes) suffix2 <- "collprot_"
  if(!collapse_protein_nodes & !collapse_peptide_nodes) suffix2 <- NULL

  if(save_intermediate & (collapse_protein_nodes | collapse_peptide_nodes)) {
    utils::write.table(edgelist_coll, sep = "\t", row.names = FALSE,
                       file = paste0(result_path, "edgelist_", suffix2, suffix, ".txt"))
  }

  if (collapse_protein_nodes | collapse_peptide_nodes) {
  message("Generating graphs ...")
  graphs <- bppg::generate_graphs_from_edgelist(edgelist_coll)
  } else {
    message("Generating graphs ...")
    graphs <- bppg::generate_graphs_from_edgelist(edgelist)
  }


  if (save_intermediate) saveRDS(graphs,
                                 file = paste0(result_path, "subgraphs_", suffix2, suffix, ".rds"))

  return(graphs)
}



