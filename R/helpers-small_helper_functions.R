
#' Geometric mean
#'
#' @param x vector with numbers
#'
#' @return geometric mean of the provided data points
#' @export
#'
#' @examples # TODO
geom_mean <- function(x,useprod = FALSE) {
  n <- length(x)

  if (useprod) {
  return(prod(x)^(1/n))
  } else {
    return(exp(mean(log(x))/n))
  }
}



