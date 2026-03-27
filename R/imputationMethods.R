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
        imp_val <- min(x[-1], na.rm = TRUE) / 2
        if (is.na(imp_val)) {
            imp_val <- min_row / 2
            # row wise min value halfed, dataset specific
        }
        return(imp_val)
    }
    min_row <- matrixStats::rowMins(as.matrix(intensities), na.rm = TRUE) 
    D_imp <- t(apply(cbind(min_row, D), 1, lod)) # TODO vapply
    return(D_imp)
}