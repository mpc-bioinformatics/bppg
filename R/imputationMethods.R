#' Impute missing values with half min value per row
#' 
#' Missing values are imputated using half of the minimum observed intensity
#' per peptide row. This method is commonly used for MNAR values in proteomics data, 
#' where missingness is often associated with low abundance peptides.
#' By imputing missing values with half of the minimum observed intensity, we can provide 
#' a conservative estimate for the missing values while preserving the overall distribution of the data.
#' 
#' @param D           Numeric matrix of intensity values with missing values (NA) to be imputed.
#'                   
#' @param min_row     Numeric vector containing the minimum observed intensity for each row (peptide) in the data matrix D. 
#'                    This vector is used to calculate the imputed values for the missing entries in D.              
#'
#' @return            Numeric matrix with imputed values for missing entries, where each missing value 
#'                    is replaced by half of the minimum observed intensity for the corresponding row (peptide) in the original data matrix D.

.min2impute<- function(D, intensities) {
    lod <- function(x) {
        # row wise min value halfed, group specific
        # this is currently not the case
        if (all(is.na(x[-1]))) {
            imp_val <- x[1] / 2
            # row wise min value halfed, dataset specific
        } else {
            imp_val <- min(x[-1], na.rm = TRUE) / 2 # returns inf if all NA
        }
        return(imp_val)
    }
    min_row <- matrixStats::rowMins(as.matrix(intensities), na.rm = TRUE) 
    D_imp <- t(apply(cbind(min_row, D), 1, lod)) # TODO vapply
    return(D_imp)
}


#' Imputation method using half of the minimum observed intensity per row (peptide) to impute missing values in proteomics data. 
#' Helper function that extracts intensity columns and applies the half minimum imputation method to the data. 
#' 
#' @param D data.frame of peptides.txt from MaxQuant
#'          
#' @param intensities data.frame containing only intensity columns (e.g. "Intensity." or "LFQ.intensity.")
#'                    extracted from the original data frame D using the .extractIntensities function.
#' 
#' @return D with imputed values for missing entries, where each missing value is replaced by
#'           half of the minimum observed intensity for the corresponding row (peptide) in the original data matrix D.

.colImputation <- function(D, intensities){

    for (j in seq_len(ncol(D))) {

        col_mean <- mean(D[, j], na.rm = TRUE)

        D[is.na(D[, j]), j] <- col_mean
    }

    return(D)
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
#' @param D data.frame of peptides.txt from MaxQuant
#' 

#' @return A numeric matrix with imputed intensity values replacing missing entries.

.missForest <- function(intensities) {
    # missForest imputation
    D_imp <- missForest::missForest(as.matrix(intensities), verbose = FALSE)$ximp
    return(D_imp)
}

#' Impute missing values using MLE and QRILC
#'
#' Performs missing value imputation using the
#' package {imputeLCMD}. 
#'
#' Missing at random (MAR) values are imputed using maximum
#' likelihood estimation (MLE), while missing not at random
#' (MNAR) values are imputed using quantile regression
#' imputation of left-censored data (QRILC).
#'
#' A model selector matrix is used to assign missing values to
#' either the MAR or MNAR model. In the current implementation,
#' the selector is initialized with only 1s, meaning all missing
#' values are treated as MAR.
#'
#' @param D [matrix]
#' Numeric matrix containing peptide intensity values.
#'
#' @param ... Additional arguments passed through the pipeline.
#'
#' @return A numeric matrix with imputed intensity values.
#'
#' @references
#' Based on the {imputeLCMD} package:
#' {https://bioconductor.org/packages/imputeLCMD}

.QRILC <- function(D, ...) {
    D_log <- log2(as.matrix(D))

    model.selector <- matrix(
        1,
        nrow = nrow(D_log),
        ncol = ncol(D_log)
    )

    D_imp_log <- imputeLCMD::impute.MAR.MNAR(
        dataSet.mvs = D_log,
        model.selector = model.selector,
        method.MAR = "MLE",
        method.MNAR = "QRILC",
        ...
    )

    D_imp <- 2^D_imp_log

    return(D_imp)
}

# Impute missing values using Bayesian PCA
#'
#' Performs missing value imputation using Bayesian Principal Component
#' Analysis (BPCA) from the \pkg{pcaMethods} package. Missing values are
#' estimated by modeling the latent structure of the data using principal
#' components and reconstructing incomplete observations.
#'
#' BPCA is particularly useful for high-dimensional datasets with
#' correlated features, such as proteomics intensity matrices.
#'
#' @param D [matrix]
#' Numeric matrix containing peptide intensity values.
#'
#' @param intensities [matrix]
#' Full intensity matrix passed through aggregateReplicates.
#' Currently unused in this imputation method.
#'
#' @return A numeric matrix with missing values imputed.
#'
#' @references
#' Based on the pcaMethods package:
#' https://bioconductor.org/packages/pcaMethods
#' 

.BPCA <- function(D, intensities) {
    D_log <- log2(as.matrix(D))

    fit <- pcaMethods::pca(
        D_log,
        method = "bpca",
        nPcs = 2
    )

    D_imp_log <- pcaMethods::completeObs(fit)
    D_imp <- 2^D_imp_log

    return(D_imp)
}