#' Generate graphs from a FASTA file
#'
#' @param fasta                    \strong{list of vector of chars} \cr
#'                                 A fasta file, already read into R by seqinr::read.fasta().
#' @param collapse_protein_nodes   \strong{logical} \cr
#'                                 If \code{TRUE}, the protein nodes will be collapsed.
#' @param collapse_peptide_nodes   \strong{logical} \cr
#'                                 If \code{TRUE}, the peptide nodes will be collapsed.
#' @param result_path              \strong{character} \cr
#'                                 The path where results are saved. If \code{NULL}, results are not saved.
#' @param suffix                   \strong{character} \cr
#'                                 The suffix for saving results.
#' @param save_intermediate        \strong{logical} \cr
#'                                 If \code{TRUE}, the intermediate results will also be saved.
#' @param prot_origin              \strong{character vector} \cr
#'                                 The origin of protein, e.g. organism etc.
#' @param ...                      Additional arguments to bppg::digest_fasta()
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

generate_graphs_from_FASTA <- function(fasta,
                                       collapse_protein_nodes = TRUE,
                                       collapse_peptide_nodes = TRUE,
                                       result_path = NULL,
                                       suffix = NULL,
                                       save_intermediate = FALSE,
                                       prot_origin = NULL,
                                       ...) {

  message("Digesting FASTA file ...")
  digested_proteins <- bppg::digest_fasta(fasta, ...)
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



