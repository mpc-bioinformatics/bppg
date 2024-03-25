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
#'
#'

#graph1 <- G4
#graph2 <- G5


isomorphic_bipartite <- function(graph1, graph2, ...) {

  #iso <- igraph::isomorphic(graph1, graph2)

  ## direct graphs if they are not directed yet
  if (!igraph::is_directed(graph1))   graph1 <- bppg::direct_bipartite_graph(graph1)
  if (!igraph::is_directed(graph2))   graph2 <- bppg::direct_bipartite_graph(graph2)
  #graph2 <- bppg::direct_bipartite_graph(graph2)


 # igraph::V(graph1)$color <- c("black", "white")[igraph::V(graph1)$type+1]
#  igraph::V(graph2)$color <- c("black", "white")[igraph::V(graph2)$type+1]


  iso <- igraph::isomorphic(graph1, graph2, method = "vf2")#, vertex.color1 = igraph::V(graph1)$color, vertex.color2 = igraph::V(graph2)$color)


  # if(iso) {
  #   ## list all attributes except "type" to remove them before comparing the graphs
  #
  #
  #   cG1 <- igraph::canonical_permutation(graph1)#, colors  = igraph::V(graph1)$type)
  #   cG1 <- igraph::permute(graph1, cG1$labeling)
  #
  #   cG2 <- igraph::canonical_permutation(graph2)#, colors  = igraph::V(graph2)$type)
  #   cG2 <- igraph::permute(graph2, cG2$labeling)
  #
  #
  #   attribute_list <- unique(c(igraph::vertex_attr_names(graph1), igraph::vertex_attr_names(graph2)))
  #   attribute_list <- attribute_list[!(attribute_list %in% c("type", "color"))]
  #
  #
  #   # if there are any attributes other than "type", they will be removed
  #   if (length(attribute_list) > 0) {
  #     for (i in 1:length(attribute_list)) {
  #       cG1 <- try({delete_vertex_attr(cG1, name = attribute_list[[i]])})
  #       cG2 <- try({delete_vertex_attr(cG2, name = attribute_list[[i]])})
  #     }
  #   }
  #
  #   iso <- igraph::identical_graphs(cG1, cG2, attrs = FALSE)#all(igraph::V(cG1)$type == igraph::V(cG2)$type)
  #
  # }
  return(iso)
}




#' Title
#'
#' @param bip_graph
#' @param from_type TODO
#'
#' @return a bipartite graph that is know directed
#' @export
#'
#' @examples
#' # TODO
direct_bipartite_graph <- function(bip_graph, from_type = FALSE){


  # turn undirected into directed edges
  bip_graph <- as.directed(bip_graph, mode = "arbitrary")

  from_vertices <- igraph::V(bip_graph)[igraph::V(bip_graph)$type == from_type]
  to_vertices <- igraph::V(bip_graph)[igraph::V(bip_graph)$type == !from_type]

  # reverse edges going from the "to-group" to the "from-group"
  bip_graph <- reverse_edges(bip_graph, igraph::E(bip_graph)[to_vertices %->% from_vertices])

  return(bip_graph)
}


