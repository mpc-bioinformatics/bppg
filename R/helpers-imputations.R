#' Impute missing values with half min value per row
#' @param D           \strong{data.matrix} \cr
#'                    Matrix to impute
#' @param min_row     \strong{float vector} \cr
#'                    possible to set unique min value
#'
#' @return            vector with imputation values for each row (peptide)
#' @export

min_2_impute <- function(D, min_row){
  lod <- function(x) {
    imp_val <- min(x, na.rm = TRUE) / 2 # row wise min value halfed, group specific
    if(is.na(imp_val)){
      imp_val <- min_row / 2   # row wise min value halfed, dataset specific
    }
    return(imp_val)
  }

  D_imp <- t(apply(cbind(D, min_row), 1, lod))
  return(D_imp)
}