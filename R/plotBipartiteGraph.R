
### TODO: labelling of the nodes (letters/numbers or keep or )

#### TODO: Farbskala für die Peptid-Knoten einbauen, um die Peptid-Ratios darzustellen (Studienprojekt)

#' Plotting of bipartite peptide-protein graphs.
#'
#' @param G                         \strong{igraph graph object} \cr
#'                                  A bipartite peptide-protein graph.
#' @param vertex.label.dist         \strong{numeric} \cr
#'                                  The distance of the label from center of the vertex (0 = centered in vertex).
#' @param legend                    \strong{logical} \cr
#'                                  If \code{TRUE}, a legend will be added.
#' @param vertex.color              \strong{character vector} \cr
#'                                  The colours for the different vertex types.
#' @param vertex.size               \strong{numeric} \cr
#'                                  The size of vertices.
#' @param vertex.label.cex          \strong{numeric} \cr
#'                                  The size of vertex labels.
#' @param edge.width                \strong{numeric} \cr
#'                                  The width of the edges.
#' @param vertex.size2              \strong{numeric} \cr
#'                                  The vertex size 2.
#' @param useCanonicalPermutation   \strong{logical} \cr
#'                                  If \code{TRUE}, the graph will be converted into the canonical permutation before plotting.
#' @param three_shapes              \strong{logical} \cr
#'                                  If \code{TRUE}, a separate shape will be used for the unique peptides.
#' @param imputed_encoding          \strong{logical} \cr
#'                                  If \code{TRUE}, imputed values will be colored in the third color
#' @param node_labels_proteins      \strong{character} \cr
#'                                  The type of labels for the proteins. Options are "letters" or "accessions".
#' @param node_labels_peptides      \strong{character} \cr
#'                                  The type of labels for the peptides. Options are "numbers" or "pep_ratios" or "log_pep_ratios" or "pep_ratio_aggr".
#' @param round_digits              \strong{integer} \cr
#'                                  The number of digits to round the peptide ratios to.
#' @param use_edge_attributes       \strong{logical} \cr
#'                                  If \code{TRUE}, edge attributes will be used for plotting (e.g. deleted edges will be dashed)
#' @param legend.x                  \strong{numeric or character} \cr
#'                                  The x-coordinate of the legend or a keyword for the position. See [graphics::legend()] for details.
#' @param legend.y                  \strong{numeric or character} \cr
#'                                  The y-coordinate of the legend or a keyword for the position. See [graphics::legend()] for details.
#' @param ...                       Additional arguments for plot.igraph.
#'
#' @return Plot of one bipartite graph.
#' @export
#'
#' @examples
#' biadjacency_matrix <- matrix(c(1,1,1,0), nrow = 2)
#' G <- igraph::graph_from_biadjacency_matrix(biadjacency_matrix)
#' plotBipartiteGraph(G, three_shapes = TRUE, useCanonicalPermutation = TRUE)

plotBipartiteGraph <- function(G, vertex.label.dist = 0, legend = TRUE,
                               vertex.color = c("mediumseagreen", "cadetblue2", "coral1"),
                               vertex.size = 15, vertex.label.cex = 1, edge.width = 1, vertex.size2=15,
                               useCanonicalPermutation = FALSE, three_shapes = FALSE,
                               imputed_encoding = FALSE,
                               node_labels_proteins = "letters",
                               node_labels_peptides = "numbers",
                               round_digits = 2, use_edge_attributes = FALSE,
                               legend.x = "bottom", legend.y = NULL,
                               ...) {

  igraph::V(G)$type <- !igraph::V(G)$type           # switch node types so that proteins are at the top
  # 0 = proteins, 1 = peptides

  if (useCanonicalPermutation) {
    cG <- igraph::canonical_permutation(G)
    G <- igraph::permute(G, cG$labeling)
  }
  Layout <- igraph::layout_as_bipartite(G)
  names_G <- character(length(igraph::V(G)))

  pos_proteins <- Layout[, 1][Layout[, 2] == 1]
  pos_peptides <- Layout[, 1][Layout[, 2] == 0]

  if (node_labels_proteins == "letters") {
    #### TODO: was ist, wenn es mehr als 26 Proteine gibt?
    names_G[Layout[, 2] == 1] <- LETTERS[rank(pos_proteins)]
  }
  if (node_labels_proteins == "accessions") {
    names_G[Layout[, 2] == 1] <- limma::strsplit2(igraph::V(G)$name[Layout[, 2] == 1], ";")[, 1]
  }
  # nicht geordnete Zahlen
  if (node_labels_proteins == "numbers_noord") {
    names_G[Layout[, 2] == 1] <- 1:length(pos_proteins)
  }



  if (node_labels_peptides == "numbers") {
    names_peptides <- 1:sum(Layout[, 2] == 0)
    names_G[Layout[, 2] == 0] <- names_peptides[rank(pos_peptides)]
  }
  if (node_labels_peptides == "pep_ratios") {
    pep_ratios <- igraph::V(G)$pep_ratio
    names_G[Layout[, 2] == 0] <- round(pep_ratios[Layout[, 2] == 0], round_digits)
  }
  if (node_labels_peptides == "log_pep_ratios") {
    pep_ratios <- log2(igraph::V(G)$pep_ratio)
    names_G[Layout[, 2] == 0] <- round(pep_ratios[Layout[, 2] == 0], round_digits)
  }
  if (node_labels_peptides == "pep_ratio_aggr") {
    pep_ratios <- igraph::V(G)$pep_ratio_aggr
    names_G[Layout[, 2] == 0] <- round(pep_ratios[Layout[, 2] == 0], round_digits)
  }
  if (node_labels_peptides == "") {
    names_G[Layout[, 2] == 0] <- NA
  }


  G <- igraph::set_vertex_attr(G, name = "name", value = names_G)

  #################################

  type <- integer(length(igraph::V(G)))
  type[!igraph::V(G)$type] <- 1                          # "protein"
  type[igraph::V(G)$type] <- 2                           # "shared peptide"
  type[igraph::V(G)$type & igraph::degree(G) == 1] <- 3  # "unique peptide"

  if (imputed_encoding){
    type2 <- integer(length(igraph::V(G)))
    type2[!igraph::V(G)$type] <- 1                           # "protein"
    type2[igraph::V(G)$type] <- 2                            # "peptide"
    type2[igraph::V(G)$imputed] <- 3   # "imputed ratio"
  } else {
    type2 <- type
  }

  if (three_shapes) {
    mydiamond <- function(coords, v = NULL, params) {
      vertex.color <- params("vertex", "color")
      if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
      }
      vertex.size <- 1 / 200 * params("vertex", "size")
      if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
      }

      graphics::symbols(x = coords[, 1], y = coords[, 2], bg = vertex.color,
              star = 1.2 * cbind(vertex.size, vertex.size, vertex.size, vertex.size),
              add = TRUE, inches = FALSE)
    }
    igraph::add_shape("diamond", clip = igraph::shape_noclip,
                      plot = mydiamond)
    vertex.shapes = c("circle", "crectangle", "diamond")[type]
  } else {
    vertex.shapes = c("circle", "crectangle")[igraph::V(G)$type + 1]
  }

  #if (legend) graphics::par(mar = c(10, 4, 4, 2) + 0.1)

  if (use_edge_attributes) {
    edge.lty <- igraph::E(G)$deleted + 1
  } else {
    edge.lty <- 1
  }

  plot(G, layout = igraph::layout_as_bipartite, vertex.color = vertex.color[type2],
       vertex.shape = vertex.shapes,
       vertex.label.degree = c(-pi / 2, pi / 2)[igraph::V(G)$type + 1],
       vertex.label.dist = vertex.label.dist,
       vertex.size = vertex.size, vertex.label.cex = vertex.label.cex,
       edge.width = edge.width, vertex.size2=vertex.size2, edge.lty = edge.lty, ...)

  ## TODO add imputation to legend
  if (legend && three_shapes) {
    if (imputed_encoding){
      legend(x = legend.x, y = legend.y, legend = c("protein", "shared peptide", "unique peptide", "not imputed peptide", "imputed peptide", "imputed protein"),
           col = c(vertex.color[1], "black", "black", vertex.color[2:3], vertex.color[3]), pch = c(19, 0, 5, 20, 20, 19))
    } else {
      legend(x = legend.x, y = legend.y, legend = c("protein", "shared peptide", "unique peptide"),
           col = vertex.color, pch = c(19, 15, 18))
    }
  }
  if (legend && !three_shapes) {
    if (imputed_encoding){
      legend(x = legend.x, y = legend.y, legend = c("protein", "shared peptide", "unique peptide", "not imputed peptide", "imputed peptide", "imputed protein"),
           col = c(vertex.color[1], "black", "black", vertex.color[2:3], vertex.color[3]), pch = c(19, 0, 0, 20, 20, 19))
    } else {
    legend(x = legend.x, y = legend.y, legend = c("protein", "shared peptide", "unique peptide"),
           col = c(vertex.color[1], "black", "black", vertex.color[2:3]), pch = c(19, 15, 15))
    }
  }
}
