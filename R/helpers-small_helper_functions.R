
#' Geometric mean
#'
#' @param x vector with numbers
#'
#' @return geometric mean of the provided data points
#' @export
#'
#' @examples # TODO
geom_mean <- function(x) {
  n <- length(x)
  return(prod(x)^(1/n))
}



