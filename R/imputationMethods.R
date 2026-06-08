#' Impute missing values with half min value per row
#' @param D           \strong{data.matrix} \cr
#'                    Matrix to impute
#' @param min_row     \strong{float vector} \cr
#'                    possible to set unique min value
#'
#' @return            vector with imputation values for each row (peptide)
#' @export

min_2_impute <- function(D, intensities) {
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


#' Extract intensity columns from MaxQuant dataframe
#'
#' Helper function that extracts intensity columns and optionally
#' renames them by removing a given pattern.
#'
#' @param D data.frame of peptides.txt from MaxQuant
#' @param col_mean calculates mean of each column 
#' (e.g. "Intensity." or "LFQ.intensity.")
#' @return data.frame with only intensity columns

col_imputation <- function(D, intensities){

    for (j in seq_len(ncol(D))) {

        col_mean <- mean(D[, j], na.rm = TRUE)

        D[is.na(D[, j]), j] <- col_mean
    }

    return(D)
}


# missForest 
#' Imputation method using random forest to predict missing values based on observed data.
#' @param D data.frame of peptides.txt from MaxQuant
#' @param D_imp applies missForest imputation to the data.frame D, 
#'which contains only intensity columns 
#' (e.g. "Intensity." or "LFQ.intensity.")
#' @return data.frame with only intensity columns

missForest <- function(D, intensities) {
    # missForest imputation
    D_imp <- missForest::missForest(D, verbose = FALSE)$ximp
    return(D_imp)
}

# QRILC_MAR_MNAR
#' Imputation method using QRILC for MNAR values and KNN for MAR values
#' @param D data.frame of peptides.txt from MaxQuant
#' @return data.frame with only intensity columns
#' QRILC_MAR_MNAR applies imputation using QRILC for MNAR values and KNN for MAR values.
#' This method is designed to handle both types of missingness in proteomics data, where MAR values are imputed using KNN and MNAR values are imputed using QRILC.
#' The function takes a data.frame D containing only intensity columns and returns a data.frame with imputed values for both MAR and MNAR missingness.
#' The model.selector matrix is used to specify which imputation method to apply to each missing value, with 1 indicating that KNN should be used for MAR values and 0 indicating that QRILC should be used for MNAR values. The imputeLCMD::impute.MAR.MNAR function is then called to perform the imputation based on the specified methods.
#' The resulting imputed data.frame is returned as the output of the function.

QRILC_MAR_MNAR <- function(D, ...) {

    model.selector <- matrix(
        1,
        nrow = nrow(D),
        ncol = ncol(D)
    )

    D_imp <- imputeLCMD::impute.MAR.MNAR(
        dataSet.mvs = as.matrix(D),
        model.selector = model.selector,
        method.MAR = "MLE",
        method.MNAR = "QRILC"
    )

    return(D_imp)
}

# BPCA
#' Imputation method using Bayesian Principal Component Analysis (BPCA) to predict missing values based on observed data.
#' @param D data.frame of peptides.txt from MaxQuant
#' @param D_imp applies BPCA imputation to the data.frame D, which contains only intensity columns 
#' (e.g. "Intensity." or "LFQ.intensity.")
#' @return data.frame with only intensity columns
#' 

BPCA <- function(D, intensities) {
    D_imp <- pcaMethods::pca(
        as.matrix(D), 
        method = "bpca", 
        nPcs = 2
        )

    D_imp <- pcaMethods::completeObs(D_imp)
    
    return(D_imp)
}