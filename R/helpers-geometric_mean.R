
#' Geometric mean
#'
#' @param x
#' @param na.rm
#' @param log_base
#'
#' @return
#' @export
#'
#' @examples
geom_mean <- function(x, na.rm = FALSE, log_base = 2) {
  res <- log_base^mean(log(x, base = log_base), na.rm = na.rm)
  return(res)
}
