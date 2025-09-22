#' Conversion of submatrices to subgraphs.
#'
#' @param x   \strong{matrix} \cr
#'            An element of a submatrix list.
#'
#' @return A graph as igraph object.
#'
#'
#' @examples
#' M <- matrix(c(1,0,1,1), nrow = 2, byrow = TRUE)
#' bppg:::.convertToBipartiteGraph(M)

.convertToBipartiteGraph <- function(x) {
    ## class list if it contains peptide ratios
    if ("list" %in% class(x)) { 
        S <- x$X
    } else   {
        S <- x
    }

    G <- igraph::graph_from_biadjacency_matrix(S)
    return(G)
}
