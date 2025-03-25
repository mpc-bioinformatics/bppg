#' Enchanced version of the igraph::isomorphic function that also considers the
#' node type in bipartite graphs, e.g. that W- and M-shaped graphs are NOT isomorphic
#'
#' @param graph1   \strong{graph (igraph)} \cr
#'                 First graph.
#' @param graph2   \strong{graph (igraph)} \cr
#'                 Second graph.
#' @param ...      currently unused
#'
#' @return TRUE if graphs are isomorphic, FALSE if not.
#' @export
#'
#' @seealso [generate_graphs_from_edgelist()]
#'
#' @examples
#'

isomorphic_bipartite <- function(graph1, graph2, ...) {

  ## direct graphs if they are not directed yet
  if (!igraph::is_directed(graph1))   graph1 <- bppg::direct_bipartite_graph(graph1)
  if (!igraph::is_directed(graph2))   graph2 <- bppg::direct_bipartite_graph(graph2)

  iso <- igraph::isomorphic(graph1, graph2, method = "vf2")

  return(iso)
}




#' Transform a bipartite graph into a directed graph.
#'
#' @param bip_graph   \strong{graph (igraph)} \cr
#'                    A bipartite graph.
#' @param from_type   \strong{logical} \cr
#'                    If \code{TRUE}, the edges will go out from the vertices with the type \code{TRUE} from the bipartite graph.
#'
#' @return A bipartite graph that is know directed.
#' @export
#'
#' @examples
#'
#' @importFrom igraph %->%
#'

direct_bipartite_graph <- function(bip_graph, from_type = FALSE){


  # turn undirected into directed edges
  bip_graph <- igraph::as.directed(bip_graph, mode = "arbitrary")

  from_vertices <- igraph::V(bip_graph)[igraph::V(bip_graph)$type == from_type]
  to_vertices <- igraph::V(bip_graph)[igraph::V(bip_graph)$type == !from_type]

  # reverse edges going from the "to-group" to the "from-group"
  bip_graph <- igraph::reverse_edges(bip_graph, igraph::E(bip_graph)[to_vertices %->% from_vertices])

  return(bip_graph)
}


