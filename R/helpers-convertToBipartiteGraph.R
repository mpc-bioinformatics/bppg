#' Conversion of submatrices to subgraphs.
#'
#' @param x   \strong{matrix} \cr
#'            An element of a submatrix list.
#'
#' @return A graph as igraph object.
#' @export
#'
#' @examples
#' M <- matrix(c(1,0,1,1), nrow = 2, byrow = TRUE)
#' bppg::.convertToBipartiteGraph(M)

.convertToBipartiteGraph <- function(x) {
  if ("list" %in% class(x)) { # class list if it contains peptide ratios
    S <- x$X
  } else   {
    S <- x
  }

  G <- igraph::graph_from_incidence_matrix(S)
  return(G)
}
