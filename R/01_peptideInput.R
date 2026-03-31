# Functions in this file:
# .extractIntensities
# readMqPeptideTable
# readSpecPeptideTable
#


#' Helper function that extracts the itensitie columns and columns of interest 
#' from a given dataframe.
#'
#' @param D                         \strong{data.frame} \cr
#'                                  Data frame of peptides.txt from MaxQuant
#' @param col_pattern               \strong{character} \cr
#'                                  Pattern to recognize intensity columns by.
#'                                  Should be "Intensity." or "LFQ.intensity."
#' @param rename_columns            \strong{logical} \cr
#'                                  If \code{TRUE}, "Intensity." or
#'                                  "LFQ.intensity." are removed
#' @return returns intensity dataframe
#'
.extractIntensities <- function(D, col_pattern, rename_columns){
    intensities <- D[, grep(col_pattern, colnames(D))]
    if (rename_columns) {
        colnames(intensities) <- stringr::str_replace(colnames(intensities),
            col_pattern, "")
    }
    return(intensities)
}

#' Import of MaxQuant's peptide.txt-table.
#'
#' @param path                      \strong{character} \cr
#'                                  The path to the peptides.txt table.
#' @param group                     \strong{character} \cr
#'                                  List or vector of group names corresponding
#'                                  to the order of samples.
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
#' @param verbose                   \strong{logical} \cr
#'                                  If \code{TRUE}, additional information on
#'                                  each iteration of the optimization is 
#'                                  printed
#' @return A SummarizedExperiment with intensities, sequences, and optional data
#'         for the rowData dataframe.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "peptides.txt", package = "bppg")
#' D <- readMqPeptideTable(path = file, LFQ = TRUE, remove_contaminants = FALSE)

readMqPeptideTable <- function(path, group = NULL, LFQ = FALSE,
    remove_contaminants = FALSE,
    rename_columns = TRUE, zeroToNA = TRUE,
    remove_empty_rows = TRUE,
    further_columns_to_keep = NULL,
    verbose = FALSE) {
    checkmate::assertFileExists(path, access = "", extension = NULL)
    checkmate::assertFlag(LFQ)
    checkmate::assertFlag(remove_contaminants)
    checkmate::assertFlag(rename_columns)
    checkmate::assertFlag(zeroToNA)
    checkmate::assertFlag(remove_empty_rows)
    checkmate::assertVector(further_columns_to_keep, null.ok = TRUE)
    checkmate::assertFlag(verbose)

    D <- utils::read.table(path, sep = "\t", header = TRUE)
    rownames(D) <- D$Sequence

    ## remove decoy entries:
    ind_decoy <- D$Reverse == "+"
    D <- D[!ind_decoy, ]
    if (verbose) print(paste0("Removed ", sum(ind_decoy), " decoy sequences."))

    ind_cont <- D$Potential.contaminant == "+"
    if (remove_contaminants) {
        D <- D[!ind_cont, ]
        if (verbose) print(paste0("Removed ", sum(ind_cont),
                " contaminant sequences."))
    }

    if (LFQ) {
        intensities <- .extractIntensities(D, "LFQ.intensity.", rename_columns)
    } else {
        intensities <- .extractIntensities(D, "Intensity.", rename_columns)
    }

    if (zeroToNA) {
        intensities[intensities == 0] <- NA
        if (remove_empty_rows) {
            validvalues <- rowSums(!is.na(intensities))
            D <- D[validvalues >= 1, ]
            intensities <- intensities[validvalues >= 1, ]
        }
    }

    if(is.null(group)){
        colDF <- data.frame(sample = colnames(intensities))
    } else {
        colDF <- data.frame(sample = colnames(intensities), group = group)
    }

    if (is.null(further_columns_to_keep)) {
        rowDF <- data.frame(Sequence = D$Sequence)
    } else {
        further_columns <- D[, further_columns_to_keep, drop = FALSE]
        colnames(further_columns) <- further_columns_to_keep
        rowDF <- data.frame(Sequence = D$Sequence, further_columns)
    }
    return(SummarizedExperiment::SummarizedExperiment(
        assays = list(intensities=intensities), 
        colData = colDF, rowData = rowDF))
}


#' Import of Spectronauts's peptide_quant-table, with raw.PEP.Quantity values.
#'
#' @param path                      \strong{character} \cr
#'                                  The path to the peptides.txt table
#' @param group                     \strong{character} \cr
#'                                  List or vector of group names corresponding
#'                                  to the order of samples.
#' @param remove_contaminants       \strong{logical} \cr
#'                                  If \code{TRUE}, peptide sequences from 
#'                                  potential contaminants are removed
#' @param remove_decoys             \strong{logical} \cr
#'                                  If \code{TRUE}, decoy peptides are removed
#' @param rename_columns            \strong{logical} \cr
#'                                  If \code{TRUE}, "raw.PEP.Quantity" are removed
#' @param cut_off                   \strong{integer} \cr
#'                                  Values below this threshold will be set to zero
#' @param zeroToNA                  \strong{logical} \cr
#'                                  If \code{TRUE}, zeros are converted to NAs. Should be NA for downstream analysis
#' @param remove_empty_rows         \strong{logical} \cr
#'                                  If \code{TRUE}, rows with only NAs are removed.
#' @param further_columns_to_keep   \strong{integer vector} \cr
#'                                  Indices of additional columns to keep, except peptide sequence and intensities
#'
#' @return A data frame with sequences and intensities.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "spec_peptides.tsv", package = "bppg") # TODO
#' D <- readSpecPeptideTable(path = file, remove_contaminants = FALSE)

readSpecPeptideTable <- function(path, group = NULL, remove_contaminants = FALSE,
    remove_decoys = TRUE, rename_columns = TRUE,
    cut_off = 1000, zeroToNA = TRUE,
    remove_empty_rows = TRUE, 
    further_columns_to_keep = NULL,
    verbose = FALSE) {
    checkmate::assertFileExists(path, access = "", extension = NULL)
    checkmate::assertFlag(remove_contaminants)
    checkmate::assertFlag(rename_columns)
    checkmate::assertNumeric(cut_off, lower = 0)
    checkmate::assertFlag(zeroToNA)
    checkmate::assertFlag(remove_empty_rows)
    checkmate::assertVector(further_columns_to_keep, null.ok = TRUE)
    checkmate::assertFlag(verbose)
    
    D <- utils::read.table(path, sep = "\t", header = TRUE)

    # all columns in Spectronaut are optional
    # need to check which ones are there/ communicate which ones have to be
    # Quantification Data Filtering removes decoys in Spectronaut
    # die spalte sorgt für die duplication, da EG statt PEP
    if (remove_decoys) {
        ind_decoy <- D$EG.IsDecoy == "True"
        D <- D[!ind_decoy, ]
        if (verbose) print(paste("Removed", sum(ind_decoy), "decoy sequences."))
    }

  # remove duplicates
    ind_dub <- duplicated(D)
    D <- D[!ind_dub, ]
    intensities <- D[, grep("raw.PEP.Quantity", colnames(D))]
    rownames(intensities) <- D$PEP.GroupingKey

    # structure: [1] C1_R1.raw.PEP.Quantity zu X.1..C1_R1.raw.PEP.Quantity, leave sample name
    if (rename_columns) colnames(intensities) <- lapply(colnames(intensities), 
        FUN = function(x) {stringr::str_split(x, "\\.")[[1]][4]})

    # valid intensity cut off, values too low tend to be false positive identifications
    low_intensity <- intensities < cut_off
    intensities[low_intensity] <- 0
    if (verbose) print(paste("Removed", sum(low_intensity, na.rm = TRUE), 
            "intensities below:", cut_off))

    # "filtered" also as an option?
    if (zeroToNA) {
        intensities[intensities == 0] <- NA
        if (remove_empty_rows) {
            validvalues <- rowSums(!is.na(intensities))
            ind_full <- validvalues >= 1
            D <- D[ind_full, ]
            intensities <- intensities[ind_full, ]
        if (verbose) print(paste("Removed", sum(!ind_full, na.rm = TRUE),
            "empty rows"))
        }
    }

    if(is.null(group)){
        colDF <- data.frame(sample = colnames(intensities))
    } else {
        colDF <- data.frame(sample = colnames(intensities), group = group)
    }
    if (is.null(further_columns_to_keep)) {
        rowDF <- data.frame(Sequence = rownames(intensities))
    } else {
        further_columns <- D[, further_columns_to_keep, drop = FALSE]
        colnames(further_columns) <- further_columns_to_keep
        rowDF <- data.frame(Sequence = rownames(intensities), further_columns)
    }

    RES <- SummarizedExperiment::SummarizedExperiment(
        assays = list(intensities = intensities),
        colData = colDF, rowData = rowDF)
    
    return(RES)
}