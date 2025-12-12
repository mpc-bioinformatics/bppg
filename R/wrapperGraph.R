#' Functions in this file:
#' generateGraphsFromFASTA
#' generateGraphsFromQuantData

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
    edgelist <- bppg::digestFASTA(fasta, protOrigin = prot_origin, ...)
    if (save_intermediate) {
        message("Saving edgelist ...")
        utils::write.table(edgelist, sep = "\t", row.names = FALSE,
            file = file.path(result_path, paste0("edgelist_", suffix, ".txt")))
    }

    graphs <- generateGraphsFromEdgelist(edgelist, collProtNodes, collPeptNodes)

    if (collProtNodes && collPeptNodes) suffix2 <- "collprotpept_"
    if (collPeptNodes && !collProtNodes) suffix2 <- "collpept_"
    if (collProtNodes && !collPeptNodes) suffix2 <- "collprot_"
    if (!collProtNodes && !collPeptNodes) suffix2 <- NULL

    if (save_intermediate) {
        saveRDS(graphs, file = file.path(result_path, paste0("subgraphs_", suffix2, suffix, ".rds")))
    }
    return(graphs)
}


#' Generate graphs from quantitative peptide-level data
#'
#' @param D                        \strong{data.frame} \cr
#'                                 A data set with peptide sequence as first
#'                                 column and peptide intensities in subsequent
#'                                 columns, e.g. created with
#'                                 [bppg::readMqPeptideTable()].
#' @param fasta                    \strong{list of vector of characters} \cr
#'                                 A fasta file used for identification of
#'                                 peptides in already read into R by
#'                                 [seqinr::read.fasta()].
#' @param outpath                  \strong{character} \cr
#'                                 The output path for the results.
#' @param missed_cleavages         \strong{integer} \cr
#'                                 The number of allowed missed cleavages
#'                                 in a peptide.
#' @param min_aa                   \strong{integer} \cr
#'                                 The minimum number of amino acids
#'                                 in a peptide.
#' @param max_aa                   \strong{integer} \cr
#'                                 The maximum number of amino acids
#'                                 in a peptide.
#' @param id_columns               \strong{integer vector} \cr
#'                                 The columns of D that contain ID information
#'                                 (the rest should contain only peptide
#'                                 intensities, properly normalized).
#' @param seq_column               \strong{character} \cr
#'                                 The column name of the column with the
#'                                 peptide sequences.
#' @param collProtNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the protein nodes
#'                                 will be collapsed.
#' @param collPeptNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the peptide nodes
#'                                 will be collapsed.
#' @param suffix                   \strong{character} \cr
#'                                 The suffix for output files.
#' @param ...                      currently not in use
#'
#' @return A list of list of graphs
#' @export
#'
#' @seealso [bppg::readMqPeptideTable()], [seqinr::read.fasta()],
#'          [generateQuantGraphs()], [bppg::generateGraphsFromFASTA()]
#'
#' @examples
#' TODO!!!

generateGraphsFromQuantData <- function(D,
    fasta,
    outpath = NULL,
    #normalize = FALSE,
    missed_cleavages = 2,
    min_aa = 6,
    max_aa = 50,
    id_columns = 1,
    seq_column = "Sequence",
    collProtNodes = TRUE,
    collPeptNodes = FALSE,
    suffix = "",
    ...) {
    message("Digesting FASTA file...")
    edgelist <- bppg::digestFASTA(fasta,
        missed_cleavages = missed_cleavages,
        min_aa = min_aa, max_aa = max_aa)

    if (!is.null(outpath)) {
        openxlsx::write.xlsx(edgelist, file = paste0(outpath,
                "edgelist_fasta_", suffix, ".xlsx"),
            overwrite = TRUE, keepNA = TRUE)
    }

    ## aggregate replicates by calculating the mean
    group <- factor(limma::strsplit2(colnames(D), split = "_")[, 1])
    D_aggr <- bppg::aggregateReplicates(D, method = "mean", missing.limit = 0.4,
        group = group, id_col = id_columns)

    if (!is.null(outpath)) {
        openxlsx::write.xlsx(SummarizedExperiment::assays(D_aggr)$intensities, 
            file = paste0(outpath, "aggr_peptides_", suffix, ".xlsx"),
            overwrite = TRUE, keepNA = TRUE)
    }

    ## calculate the peptide ratio table
    groups  <- levels(group)
    peptide_ratios <- bppg::calculatePeptideRatios(D = D_aggr,
        group_levels = groups)
    if (!is.null(outpath)) {
        openxlsx::write.xlsx(
            SummarizedExperiment::assays(peptide_ratios)$logRatios, 
            file = paste0(outpath,"peptide_ratios_", suffix, ".xlsx"),
            overwrite = TRUE, keepNA = TRUE)
    }

    ## Generierung der Graphen (man braucht peptide_ratios und fast_edgelist!)
    graphs <- generateQuantGraphs(exp_peptide_ratios = peptide_ratios,
        fasta_edgelist = edgelist,
        outpath = outpath, seq_column = seq_column,
        collProtNodes = collProtNodes,
        collPeptNodes = collPeptNodes,
        suffix = suffix)
    return(graphs)

}
