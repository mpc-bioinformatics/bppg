
### TODO: cat funktioniert nicht einwandfrei
### TODO; Matrix package wird auf jeden Fall benötigt

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
#' #' library(seqinr)
#' file <- system.file("extdata", "2020_01_31_proteome_S_cerevisae.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' graphs <- bppg::generate_graphs_from_FASTA(fasta)
#'
generate_graphs_from_FASTA <- function(fasta, collapse_protein_nodes = TRUE,
                                       collapse_peptide_nodes = TRUE,
                                       result_path = NULL,
                                       suffix = NA, save_intermediate = FALSE,
                                       ...) {

  cat("Digesting FASTA file...")
  digested_proteins <- bppg::digest_fasta(fasta)#, ...)
  cat("Generating graphs ...")
  edgelist <- bppg::generate_edgelist(digested_proteins)
  if(save_intermediate) utils::write.table(edgelist, sep = "\t", row.names = FALSE,
                                    file = paste0(result_path, "edgelist_", suffix, ".txt"))

  graphs <- bppg::generate_graphs_from_edgelist(edgelist)
  if(save_intermediate) saveRDS(graphs,
                                    file = paste0(result_path, "subgraphs_", suffix, ".rds"))

  if (collapse_protein_nodes) {
    cat("Collapsing protein nodes ...")
    graphs <- bppg::collapse_protein_nodes(graphs, sparse = TRUE, fast = TRUE, fc = FALSE)
    if(save_intermediate) saveRDS(graphs,
                                  file = paste0(result_path, "subgraphs_coll_prot_", suffix, ".rds"))
  }

  if (collapse_peptide_nodes) {
    cat("Collapsing peptide nodes ...")
    graphs <- bppg::collapse_peptide_nodes(graphs, sparse = TRUE, fast = TRUE, fc = FALSE)
    if(save_intermediate) saveRDS(graphs,
                                  file = paste0(result_path, "subgraphs_coll_prot_pep_", suffix, ".rds"))
  }

  return(graphs)
}
