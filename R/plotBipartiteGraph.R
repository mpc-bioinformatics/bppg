
### TODO: labelling of the nodes (letters/numbers or keep or )

#### TODO: Farbskala für die Peptid-Knoten einbauen, um die Peptid-Ratios darzustellen (Studienprojekt)

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
#' @param ... Additional arguments for plot.igraph.
#' @param node_labels_proteins "letters" or "acessions"
#' @param node_labels_peptides "numbers" or "pep_ratios" or "pep_ratio_aggr"
#' @param round_digits Number of digits to round the peptide ratios to.
#' @param use_edge_attributes Use edge attributes for plotting (e.g. deleted edges will be dashed)?
#'
#' @return Plot of one bipartite graph.
#' @export
#'
#' @examples
#' biadjacency_matrix <- matrix(c(1,1,1,0), nrow = 2)
#' G <- igraph::graph_from_incidence_matrix(biadjacency_matrix)
#' plotBipartiteGraph(G, three_shapes = TRUE, useCanonicalPermutation = TRUE)
plotBipartiteGraph <- function(G, vertex.label.dist = 0, legend = TRUE,
                               vertex.color = c("mediumseagreen", "cadetblue2", "coral1"),
                               vertex.size = 15, vertex.label.cex = 1, edge.width = 1, vertex.size2=15,
                               useCanonicalPermutation = FALSE, three_shapes = FALSE,
                               node_labels_proteins = "letters",
                               node_labels_peptides = "numbers",
                               round_digits = 2, use_edge_attributes = FALSE,
                               legend.x = NULL, legend.y = NULL,
                               ...) {

  igraph::V(G)$type <- !igraph::V(G)$type           # switch node types so that proteins are at the top
  # 0 = proteins, 1 = peptides

  if (useCanonicalPermutation) {
    cG <- igraph::canonical_permutation(G)
    G <- igraph::permute(G, cG$labeling)
  }
  Layout <- igraph::layout.bipartite(G)
  names_G <- character(length(igraph::V(G)))

  pos_proteins <- Layout[,1][Layout[,2] == 1]
  pos_peptides <- Layout[,1][Layout[,2] == 0]

  if (node_labels_proteins == "letters") {
    #### TODO: was ist, wenn es mehr als 26 Proteine gibt?
    names_G[Layout[,2] == 1] <- LETTERS[rank(pos_proteins)]
  }
  if (node_labels_proteins == "accessions") {
    names_G[Layout[,2] == 1] <- limma::strsplit2(V(G)$name[Layout[,2] == 1], ";")[,1]
  }
  # nicht geordnete Zahlen
  if (node_labels_proteins == "numbers_noord") {
    names_G[Layout[,2] == 1] <- 1:length(pos_proteins)
  }



  if (node_labels_peptides == "numbers") {
    names_peptides <- 1:sum(Layout[,2] == 0)
    names_G[Layout[,2] == 0] <- names_peptides[rank(pos_peptides)]
  }
  if (node_labels_peptides == "pep_ratios") {
    pep_ratios <- V(G)$pep_ratio
    names_G[Layout[,2] == 0] <- round(pep_ratios[Layout[,2] == 0],round_digits)
  }
  if (node_labels_peptides == "pep_ratio_aggr") {
    pep_ratios <- V(G)$pep_ratio_aggr
    names_G[Layout[,2] == 0] <- round(pep_ratios[Layout[,2] == 0],round_digits)
  }
  if (node_labels_peptides == "") {
    names_G[Layout[,2] == 0] <- NA
  }


  G <- igraph::set_vertex_attr(G, name = "name", value = names_G)

  #################################


  # if (node_labels == "letters+numbers") {
  #   names_G[Layout[,2] == 1] <- LETTERS[rank(pos_proteins)]
  #   names_peptides <- 1:sum(Layout[,2] == 0)
  #   names_G[Layout[,2] == 0] <- names_peptides[rank(pos_peptides)]
  #
  #   G <- igraph::set_vertex_attr(G, name = "name", value = names_G)
  # }
  # if (node_labels == "peptide_ratios") {
  #   pep_ratios <- V(G)$pep_ratio
  #   names_G[Layout[,2] == 1] <- limma::strsplit2(V(G)$name[Layout[,2] == 1], ";")[,1]
  #
  #  # names_peptides <- 1:sum(Layout[,2] == 0)
  #   names_G[Layout[,2] == 0] <- round(pep_ratios[Layout[,2] == 0],2)
  #   G <- igraph::set_vertex_attr(G, name = "name", value = names_G)
  # }

  type <- integer(length(igraph::V(G)))
  type[!igraph::V(G)$type] <- 1                  # "protein"
  type[igraph::V(G)$type] <- 2                   # "shared peptide"
  type[igraph::V(G)$type & igraph::degree(G) == 1] <- 3  # "unique peptide"

  if (three_shapes) {
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
    igraph::add_shape("diamond", clip= igraph::shape_noclip,
              plot=mydiamond)
    vertex.shapes = c("circle", "crectangle", "diamond")[type]
  } else {
    vertex.shapes = c("circle", "crectangle")[igraph::V(G)$type+1]
  }

  #if (legend) graphics::par(mar = c(10, 4, 4, 2) + 0.1)

  if (use_edge_attributes) {
    edge.lty <- E(G)$deleted + 1
  } else {
    edge.lty <- 1
  }

  plot(G, layout = igraph::layout_as_bipartite, vertex.color=vertex.color[type],
       vertex.shape = vertex.shapes,
       vertex.label.degree = c(-pi/2, pi/2)[igraph::V(G)$type+1],
       vertex.label.dist = vertex.label.dist,
       vertex.size = vertex.size, vertex.label.cex = vertex.label.cex,
       edge.width = edge.width, vertex.size2=vertex.size2, edge.lty = edge.lty, ...)

  if (legend & three_shapes) {
    legend(x = legend.x, y = legend.y, legend = c("protein", "shared peptide", "unique peptide"),
           col = vertex.color, pch = c(19, 15, 18))
  }
  if (legend & !three_shapes) {
    legend(x = legend.x, y = legend.y, legend = c("protein", "shared peptide", "unique peptide"),
           col = vertex.color, pch = c(19, 15, 15))
  }
}
