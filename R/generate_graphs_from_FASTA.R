#' Generate graphs from a FASTA file
#'
#' @param fasta                    \strong{list of vector of chars} \cr
#'                                 A fasta file, already read into R by 
#'                                 seqinr::read.fasta().
#' @param collProtNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the protein nodes will
#'                                 be collapsed.
#' @param collPeptNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the peptide nodes will
#'                                 be collapsed.
#' @param result_path              \strong{character} \cr
#'                                 The path where results are saved. If 
#'                                 \code{NULL}, results are not saved.
#' @param suffix                   \strong{character} \cr
#'                                 The suffix for saving results.
#' @param save_intermediate        \strong{logical} \cr
#'                                 If \code{TRUE}, the intermediate results
#'                                 will also be saved.
#' @param prot_origin              \strong{character vector} \cr
#'                                 The origin of protein, e.g. organism etc.
#' @param ...                      Additional arguments to bppg::digestFASTA()
#'
#' @return subgraphs (i.e. connected components) from the graph generated from
#'         the FASTA file.
#' @export
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' graphs <- bppg::generateGraphsFromFASTA(fasta)
#'

generateGraphsFromFASTA <- function(fasta,
    collProtNodes = TRUE,
    collPeptNodes = TRUE,
    result_path = NULL,
    suffix = NULL,
    save_intermediate = FALSE,
    prot_origin = NULL,
    ...) {
    message("Digesting FASTA file ...")
    digested_proteins <- bppg::digestFASTA(fasta, ...)
    message("Generating edgelist ...")
    edgelist <- bppg::generateEdgelist(digested_proteins,
        prot_origin = prot_origin)
    if (save_intermediate) {
        message("Saving edgelist ...")
        utils::write.table(edgelist, sep = "\t", row.names = FALSE,
            file = paste0(result_path, "edgelist_",
                suffix, ".txt"))
    }

    if (collProtNodes || collPeptNodes) {
        message("Collapsing nodes ...")
        edgelist_coll <- .collapseEdgelist(edgelist,
            collProtNodes = collProtNodes, collPeptNodes = collPeptNodes)
    }

    if (collProtNodes && collPeptNodes) suffix2 <- "collprotpept_"
    if (collPeptNodes && !collProtNodes) suffix2 <- "collpept_"
    if (collProtNodes && !collPeptNodes) suffix2 <- "collprot_"
    if (!collProtNodes && !collPeptNodes) suffix2 <- NULL

    if (save_intermediate && (collProtNodes || collPeptNodes)) {
        utils::write.table(edgelist_coll, sep = "\t", row.names = FALSE,
            file = paste0(result_path, "edgelist_", suffix2, suffix, ".txt"))
    }

    if (collProtNodes || collPeptNodes) {
        message("Generating graphs ...")
        graphs <- .generateGraphsFromEdgelist(edgelist_coll)
    } else {
        message("Generating graphs ...")
        graphs <- .generateGraphsFromEdgelist(edgelist)
    }


    if (save_intermediate) {
        saveRDS(graphs, file = paste0(result_path, "subgraphs_",
                suffix2, suffix, ".rds"))
    }
    return(graphs)
}
