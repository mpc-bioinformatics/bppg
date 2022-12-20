#' Plotting of bipartite peptide-protein graphs.
#'
#' @param G A bipartite graph (igraph object).
#' @param vertex.label.dist Distance of the label from center of the vertex (0 = centered in vertex).
#' @param legend Add legend?
#' @param vertex.color Colours for the different vertex types.
#' @param vertex.size Size of vertices.
#' @param vertex.label.cex Size of vertex labels.
#' @param edge.width Width of the edges.
#' @param vertex.size2 Vertex size 2.
#' @param useCanonicalPermutation Convert the graph into the canonical permutation before plotting?
#' @param three_shapes Use a separate shape for the unique peptides?
#' @param node_names "letters&numbers" uses letters for protein nodes and numbers for peptide nodes,
#'                   "empty" uses no node names, "keep" keeps node names of G
#' @param ... Additional arguments for plot.igraph.
#'
#' @return Plot of one bipartite graph.
#' @export
#'
#' @examples
#' biadjacency_matrix <- matrix(c(1,1,1,0), nrow = 2)
#' G <- igraph::graph_from_incidence_matrix(biadjacency_matrix)
#' bppg::plotBipartiteGraph(G, three_shapes = TRUE, useCanonicalPermutation = TRUE)
plotBipartiteGraph <- function(G, vertex.label.dist = 0, legend = TRUE,
                               vertex.color = c("mediumseagreen", "cadetblue2", "coral1"),
                               vertex.size = 15, vertex.label.cex = 1, edge.width = 1, vertex.size2=15,
                               useCanonicalPermutation = FALSE, three_shapes = FALSE,
                               node_names = "letters&numbers",
                               ...) {

  igraph::V(G)$type <- !igraph::V(G)$type # switch node types so that proteins are at the top
  # 0 = proteins, 1 = peptides

  # calculate canonical permutation of the graph
  if (useCanonicalPermutation) {
    cG <- igraph::canonical_permutation(G)
    G <- igraph::permute(G, cG$labeling)
  }

  Layout <- igraph::layout.bipartite(G)
  names_G <- character(length(igraph::V(G)))
  pos_proteins <- Layout[,1][Layout[,2] == 1]
  pos_peptides <- Layout[,1][Layout[,2] == 0]

  ### proteins = Letters, peptides = numbers
  if (node_names == "letters&numbers") {
  names_G[Layout[,2] == 1] <- LETTERS[rank(pos_proteins)]
  names_peptides <- 1:sum(Layout[,2] == 0)
  names_G[Layout[,2] == 0] <- names_peptides[rank(pos_peptides)]
    G <- igraph::set_vertex_attr(G, name = "name", value = names_G)
  }
  ### empty node names
  if (node_names == "empty") {
    G <- igraph::set_vertex_attr(G, name = "name", value = rep("", igraph::gorder(G)))
  }



  if (three_shapes) {

    ## introduce 3rd node type for the unique peptides
    type <- integer(length(igraph::V(G)))
    type[!igraph::V(G)$type] <- 1                  # "protein"
    type[igraph::V(G)$type] <- 2                   # "shared peptide"
    type[igraph::V(G)$type & igraph::degree(G) == 1] <- 3  # "unique peptide"

    ### add diamond as a possible shape
    mydiamond <- function(coords, v=NULL, params) {
      vertex.color <- params("vertex", "color")
      if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
      }
      vertex.size <- 1/200 * params("vertex", "size")
      if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
      }

      graphics::symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
              stars=1.2*cbind(vertex.size, vertex.size, vertex.size, vertex.size),
              add=TRUE, inches=FALSE)
    }
    igraph::add_shape("diamond", clip = igraph::shape_noclip, plot=mydiamond)


    vertex.shapes = c("circle", "crectangle", "diamond")[type]
  } else {
    type <- igraph::V(G)$type+1
    vertex.shapes = c("circle", "crectangle")[type]
  }



  if (legend) graphics::par(mar = c(10, 4, 4, 2) + 0.1)
  plot(G, layout = igraph::layout_as_bipartite, vertex.color=vertex.color[type],
       vertex.shape = vertex.shapes,
       vertex.label.degree = c(-pi/2, pi/2)[igraph::V(G)$type+1],
       vertex.label.dist = vertex.label.dist,
       vertex.size = vertex.size, vertex.label.cex = vertex.label.cex,
       edge.width = edge.width, vertex.size2=vertex.size2, ...)

  if (legend) {
    if (three_shapes) {
      pch <- c(19, 15, 18)
      legend_text <- c("protein", "shared peptide", "unique peptide")
    } else {
      pch <-  c(19, 15)
      legend_text <- c("protein", "peptide")
    }


    legend(x = -1, y = -1.3, legend = legend_text, col = vertex.color, pch = pch)

    graphics::par(mar = c(5, 4, 4, 2) + 0.1)
  }
}
