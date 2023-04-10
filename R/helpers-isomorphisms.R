#' Enchanced version of the igraph::isomorphic function that also considers the
#' node type in bipartite graphs, e.g. that W- and M-shaped graphs are NOT isomorphic
#'
#' @param graph1 First graph.
#' @param graph2 Second graph.
#' @param ... currently unused
#'
#' @return TRUE if graphs are isomorphic, FALSE if not.
#' @export
#'
#' @examples # TODO
isomorphic_bipartite <- function(graph1, graph2, ...) {

  iso <- igraph::isomorphic(graph1, graph2)

  if(iso) {
    ### TODO: ist es nötig, die kanonische Permutation zu berechnen?
    cG1 <- igraph::canonical_permutation(graph1)
    cG1 <- igraph::permute(graph1, cG1$labeling)
    cG2 <- igraph::canonical_permutation(graph2)
    cG2 <- igraph::permute(graph2, cG2$labeling)

    iso <- all(igraph::V(cG1)$type == igraph::V(cG2)$type)

  }
  return(iso)
}
