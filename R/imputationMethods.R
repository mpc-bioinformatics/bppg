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