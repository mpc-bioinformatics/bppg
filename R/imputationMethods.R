#' Impute missing values with half min value per row/peptide. This approach is 
#' applicable for MNAR values in proteomics data. It considers group wise 
#' information before dataset wide information
#' 
#' @param D             \strong{data.frame} \cr
#'                      Data frame of the columns of one group, includes 
#'                      missing values needed to be imputed. Used for imputation
#'                      value first.
#' @param intensities   \strong{data.frame} \cr
#'                      Data frame of the complete dataset. Is used as a second 
#'                      reference for the imputation, if n groupspecific
#'                      information is available.
#' @return A vector with imputation values for each peptide for one group. It
#'  can be used instead of an missing value in the aggregation step of bppg.
#' 
#' @examples 
#' file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
#' D <- bppg::readMqPeptideTable(file)
#' D_norm <- bppg::normalizePeptideIntensities(D)
#'
#' df <- SummarizedExperiment::assays(D_norm)$intensities_norm
#' group <- rep(1:9, each = 3)
#'
#' imputed <- list()
#' for (g in levels(factor(group))) {
#'    imputed[[g]] <- bppg:::.min2impute(D = df[g == group],intensities = df)
#' }
#' D_min <- data.frame(imputed)
#' 
#' @importFrom matrixStats rowMins

.min2impute<- function(D, intensities) {
    lod <- function(x) {
        # row wise min value halfed x[1], group specific
        if (all(is.na(x[-1]))) {
            imp_val <- x[1] / 2
            # row wise min value halfed, dataset specific
        } else {
            imp_val <- min(x[-1], na.rm = TRUE) / 2 # returns inf if all NA
        }
        return(imp_val)
    }
    min_row <- matrixStats::rowMins(as.matrix(intensities), na.rm = TRUE) 
    D_imp <- apply(cbind(min_row, D), 1, lod) # TODO vapply
    return(D_imp)
}


#' Impute missing values with the mean of column/sample specific value. This 
#' approach is applicable for MAR values in proteomics data.
#' 
#' @param intensities   \strong{data.frame} \cr
#'                      Data frame of the complete dataset.
#' @return A data.frame with missing values imputed with the mean of each column
#' /sample.
#' 
#' @examples 
#' file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
#' D <- bppg::readMqPeptideTable(file)
#' D_norm <- bppg::normalizePeptideIntensities(D)
#'
#' df <- SummarizedExperiment::assays(D_norm)$intensities_norm
#' D_mean <- bppg:::.colMeanImputation(intensities = df)
#' 


.colMeanImputation <- function(intensities){ # at this point we expect 
    # TODO vapply
    for (j in seq_len(ncol(intensities))) {

        col_mean <- mean(intensities[, j], na.rm = TRUE)

        intensities[is.na(intensities[, j]), j] <- col_mean
    }

    return(intensities)
}


# missForest 
#' Performs missing value imputation using the \pkg{missForest} package,
#' a non-parametric imputation method based on random forests. Missing
#' values are predicted from observed values using an iterative random
#' forest model.
#'
#' This method can be applied to proteomics datasets containing missing
#' values and does not assume a specific data distribution.
#' 
#' @param intensities     \strong{data.frame} \cr
#'                        Data frame of the complete dataset.
#' 
#' @return A complete numeric matrix with imputed values replacing missing 
#' entries.
#' 
#' @examples 
#' file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
#' D <- bppg::readMqPeptideTable(file)
#' D_norm <- bppg::normalizePeptideIntensities(D)
#'
#' df <- SummarizedExperiment::assays(D_norm)$intensities_norm
#' D_mean <- bppg:::.missForest(intensities = df)
#' 
#' @importFrom missForest missForest

.missForest <- function(intensities) {
    # missForest imputation
    D_imp <- missForest::missForest(as.matrix(intensities), 
        verbose = FALSE)$ximp
    return(D_imp)
}

#' Impute missing values using quantile regression
#' imputation of left-censored data (QRILC) implemented in \pkg{imputeLCMD}
#'
#' @param intensities   \strong{data.frame} \cr
#'                      Numeric matrix containing peptide intensity values.
#'
#' @param ...           Additional arguments passed through to QRILIC.
#'
#' @return A numeric matrix with imputed intensity values.
#'
#' @examples 
#' file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
#' D <- bppg::readMqPeptideTable(file)
#' D_norm <- bppg::normalizePeptideIntensities(D)
#'
#' df <- SummarizedExperiment::assays(D_norm)$intensities_norm
#' D_mean <- bppg:::.QRILC(intensities = df)
#' 
#' @importFrom imputeLCMD impute.QRILC

.QRILC <- function(intensities, ...) {
    D_log <- log2(as.matrix(intensities))

    D_imp_log <- imputeLCMD::impute.QRILC(D_log, ...)[[1]]

    D_imp <- 2^D_imp_log

    return(D_imp)
}

# Impute missing values using Bayesian Principal Component
#' Analysis (BPCA) from the \pkg{pcaMethods} package.
#' BPCA is particularly useful for high-dimensional datasets with
#' correlated features, such as proteomics intensity matrices.
#'
#' @param intensities   \strong{data.frame} \cr
#'                      Numeric matrix containing peptide intensity values.
#'
#' @param ...           Additional arguments passed through to QRILIC.
#'
#' @return A numeric matrix with imputed intensity values.
#' 
#' @examples 
#' file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
#' D <- bppg::readMqPeptideTable(file)
#' D_norm <- bppg::normalizePeptideIntensities(D)
#'
#' df <- SummarizedExperiment::assays(D_norm)$intensities_norm
#' D_mean <- bppg:::.BPCA(intensities = df)
#' 
#' @importFrom pcaMethods pca completeObs

.BPCA <- function(intensities, ...) {
    D_log <- log2(as.matrix(intensities))

    fit <- pcaMethods::pca(
        D_log,
        method = "bpca",
        nPcs = 2, ...
    )

    D_imp_log <- pcaMethods::completeObs(fit)
    D_imp <- 2^D_imp_log

    return(D_imp)
}
