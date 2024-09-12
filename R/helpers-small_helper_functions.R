#' Geometric mean
#'
#' @param x vector with numbers
#' @param useprod if TRUE, prod(x)^(1/n) will be calculated, otherwise exp(mean(log(x)))
#'
#' @return geometric mean of the provided data points
#' @export
#'
#' @examples # TODO
geom_mean <- function(x, useprod = FALSE) {
  n <- length(x)

  if (useprod) {
  return(prod(x)^(1/n))
  } else {
    return(exp(mean(log(x))))
  }
}



