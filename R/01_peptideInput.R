#' Functions in this file:
#' readMqPeptideTable
#' 

#' Import of MaxQuant's peptide.txt-table.
#'
#' @param path                      \strong{character} \cr
#'                                  The path to the peptides.txt table
#' @param LFQ                       \strong{logical} \cr
#'                                  If \code{TRUE}, LFQ intensities are used,
#'                                  if FALSE, raw (unnormalized) intensities
#' @param remove_contaminants       \strong{logical} \cr
#'                                  If \code{TRUE}, peptide sequences from
#'                                  potential contaminants are removed
#' @param rename_columns            \strong{logical} \cr
#'                                  If \code{TRUE}, "Intensity." or 
#'                                  "LFQ.intensity." are removed
#' @param zeroToNA                  \strong{logical} \cr
#'                                  If \code{TRUE}, zeros are converted to NAs.
#' @param remove_empty_rows         \strong{logical} \cr
#'                                  If \code{TRUE}, rows with only NAs are 
#'                                  removed.
#' @param further_columns_to_keep   \strong{integer vector} \cr
#'                                  Indices of additional columns to keep, 
#'                                  except peptide sequence and intensities
#'
#' @return A data frame with sequences and intensities.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "peptides.txt", package = "bppg")
#' D <- readMqPeptideTable(path = file, LFQ = TRUE, remove_contaminants = FALSE)

readMqPeptideTable <- function(path, LFQ = FALSE, remove_contaminants = FALSE,
    rename_columns = TRUE, zeroToNA = TRUE,
    remove_empty_rows = TRUE,
    further_columns_to_keep = NULL) {
    checkmate::assertFileExists(path, access = "", extension = NULL)
    checkmate::assertFlag(LFQ)
    checkmate::assertFlag(remove_contaminants)
    checkmate::assertFlag(rename_columns)
    checkmate::assertFlag(zeroToNA)
    checkmate::assertFlag(remove_empty_rows)
    checkmate::assertVector(further_columns_to_keep, null.ok = TRUE)

    D <- utils::read.table(path, sep = "\t", header = TRUE)

    ## remove decoy entries:
    ind_decoy <- D$Reverse == "+"
    D <- D[!ind_decoy, ]
    print(paste0("Removed ", sum(ind_decoy), " decoy sequences."))

    ind_cont <- D$Potential.contaminant == "+"
    if (remove_contaminants) {
        D <- D[!ind_cont, ]
        print(paste0("Removed ", sum(ind_cont), " contaminant sequences."))
    }

    ## search for intensity columns or LFQ values
    if (LFQ) {
        intensities <- D[, grep("LFQ", colnames(D))]
        if (rename_columns) {
            colnames(intensities) <- stringr::str_replace(colnames(intensities),
                "LFQ.intensity.", "")
        }
    } else {
        intensities <- D[, grep("Intensity.", colnames(D))]
        if (rename_columns) {
            colnames(intensities) <- stringr::str_replace(colnames(intensities),
                "Intensity.", "")
        }
    }

    if (zeroToNA) {
        intensities[intensities == 0] <- NA
        if (remove_empty_rows) {
            validvalues <- rowSums(!is.na(intensities))
            D <- D[validvalues >= 1, ]
            intensities <- intensities[validvalues >= 1, ]
        }
    }

    if (is.null(further_columns_to_keep)) {
        RES <- data.frame(Sequence = D$Sequence, intensities)
    } else {
        further_columns <- D[, further_columns_to_keep, drop = FALSE]
        colnames(further_columns) <- further_columns_to_keep
        RES <- data.frame(Sequence = D$Sequence, further_columns, intensities)
    }
    return(RES)
}
