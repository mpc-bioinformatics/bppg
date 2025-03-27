#' Calculate the geometric mean.
#'
#' @param x         \strong{numeric vector} \cr
#'                  Input data.
#' @param useprod   \strong{logical} \cr
#'                  If \code{TRUE}, prod(x)^(1/n) will be calculated, otherwise exp(mean(log(x))).
#'
#' @return The geometric mean of the provided data points.
#' @export
#'
#' @examples
#' data <- c(1,6,3.5)
#' result <- geom_mean(data, useprod = FALSE)

geom_mean <- function(x, useprod = FALSE) {
  n <- length(x)

  if (useprod) {
  return(prod(x)^(1/n))
  } else {
    return(exp(mean(log(x))))
  }
}



