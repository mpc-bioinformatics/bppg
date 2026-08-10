# Functions in this file:
# setNodeLabels
# plotBipartiteGraph



#' Set names for plotting with plotBipartiteGraph.
#'
#' @param G                         \strong{igraph graph object} \cr
#'                                  A bipartite peptide-protein graph.
#' @param node_labels_proteins      \strong{character} \cr
#'                                  The type of labels for the proteins. Options
#'                                  are "letters" or "acessions".
#' @param node_labels_peptides      \strong{character} \cr
#'                                  The type of labels for the peptides. Options
#'                                  are"numbers" or "pep_ratios" or
#'                                  "pep_ratio_aggr".
#' @param round_digits              \strong{integer} \cr
#'                                  The number of digits to round the peptide
#'                                  ratios to.
#' @return Graph with updated names
#'
#' @importFrom igraph layout_as_bipartite set_vertex_attr V
#' @importFrom limma strsplit2
.setNodeLabels <- function(G, node_labels_peptides, node_labels_proteins,
    round_digits) {
    Layout <- igraph::layout_as_bipartite(G)
    names_G <- character(length(igraph::V(G)))

    pos_proteins <- Layout[, 1][Layout[, 2] == 1]
    pos_peptides <- Layout[, 1][Layout[, 2] == 0]

    if (node_labels_proteins == "letters") {
        names_G[Layout[, 2] == 1] <- LETTERS[rank(pos_proteins)]
    }
    if (node_labels_proteins == "accessions") {
        names_G[Layout[, 2] == 1] <- limma::strsplit2(
            igraph::V(G)$name[Layout[, 2] == 1], ";")[, 1]
    }
    ## nicht geordnete Zahlen
    if (node_labels_proteins == "numbers_noord") {
        names_G[Layout[, 2] == 1] <- seq_along(pos_proteins)
    }

    if (node_labels_peptides == "numbers") {
        names_peptides <- seq_len(sum(Layout[, 2] == 0))
        names_G[Layout[, 2] == 0] <- names_peptides[rank(pos_peptides)]
    }
    if (node_labels_peptides == "pep_logRatios") {
        pep_logRatios <- igraph::V(G)$pep_logRatio
        names_G[Layout[, 2] == 0] <- round(pep_logRatios[Layout[, 2] == 0],
            round_digits)
    }
    if (node_labels_peptides == "pep_ratios_mean") {
        pep_logRatios <- igraph::V(G)$pep_logRatio_mean
        names_G[Layout[, 2] == 0] <- round(pep_logRatios[Layout[, 2] == 0],
            round_digits)
    }
    if (node_labels_peptides == "") {
        names_G[Layout[, 2] == 0] <- ""
    }

    igraph::set_vertex_attr(G, name = "name", value = names_G)
}


#' Plotting of bipartite peptide-protein graphs.
#'
#' @param G                         \strong{igraph graph object} \cr
#'                                  A bipartite peptide-protein graph.
#' @param legend                    \strong{logical} \cr
#'                                  If \code{TRUE}, a legend will be added.
#' @param vertex.color              \strong{character vector} \cr
#'                                  The colours for the different vertex types.
#' @param vertex.size               \strong{numeric} \cr
#'                                  The size of vertices.
#' @param edge.width                \strong{numeric} \cr
#'                                  The width of the edges.
#' @param node_labels_proteins      \strong{character} \cr
#'                                  The type of labels for the proteins. Options
#'                                  are "letters" or "accessions".
#' @param node_labels_peptides      \strong{character} \cr
#'                                  The type of labels for the peptides. Options
#'                                  are "numbers" or "pep_logRatios" or
#'                                  "pep_logRatio_mean".
#' @param round_digits              \strong{integer} \cr
#'                                  The number of digits to round the peptide
#'                                  ratios to.
#' @param output_path               \strong{character} \cr
#'                                  file path for optional save of figure.
#' @param ...                       Additional arguments for ggsave.
#'
#' @return Plot of one bipartite graph.
#' @export
#'
#' @examples
#'
#' file <- system.file("extdata", "quantGraphs_collpept.rds", package = "bppg")
#' graphs <- readRDS(file)
#' G <- graphs$"1_2"[[2]]
#'
#' plotBipartiteGraph(G)
#'
#' @importFrom igraph add_shape canonical_permutation layout_as_bipartite 
#' @importFrom igraph permute V
#' @importFrom ggraph geom_edge_link geom_node_label geom_node_point ggraph
#' @importFrom ggraph theme_graph
#' @importFrom ggplot2 aes ggsave scale_fill_manual scale_shape_manual theme
#' 
plotBipartiteGraph <- function(G, legend = TRUE,
    vertex.color = c("mediumseagreen", "cadetblue2", "coral1"),
    vertex.size = 15, edge.width = 0.5,
    node_labels_proteins = "letters",
    node_labels_peptides = "numbers",
    round_digits = 2,
    output_path = NULL,
    ...) {
    name <- node_type <- NULL

    ## switch node types so that proteins are at the top
    ## 0 = proteins, 1 = peptides
    igraph::V(G)$type <- !igraph::V(G)$type

    G <- .setNodeLabels(G, node_labels_peptides, node_labels_proteins,
        round_digits)
    #################################

    igraph::V(G)$node_type[!igraph::V(G)$type] <- "Protein"
    igraph::V(G)$node_type[igraph::V(G)$type] <- "Shared Peptide"
    igraph::V(G)$node_type[igraph::V(G)$type
        & igraph::degree(G) == 1] <- "Unique Peptide"

    shape_values <- c(Protein = 21,
        `Shared Peptide` = 22,
        `Unique Peptide` = 23)
    color_values <- c(Protein = vertex.color[1], 
        `Shared Peptide` = vertex.color[2],
        `Unique Peptide` = vertex.color[3])

    p <- ggraph::ggraph(G, layout = "bipartite") +
        ggraph::geom_edge_link(linewidth = edge.width) + 
        ggraph::geom_node_point(size = vertex.size, 
            ggplot2::aes(fill = factor(node_type), shape = factor(node_type))) +
        ggraph::geom_node_label(ggplot2::aes(label = name, 
            fill = factor(node_type)), size = vertex.size/3, 
            show.legend = FALSE, family = "sans") + 
        ggplot2::scale_shape_manual(values = shape_values, name = "Node Type") +
        ggplot2::scale_fill_manual(values = color_values, name = "Node Type") +
        ggraph::theme_graph(base_family = "sans") +
        ggplot2::theme(legend.position = if (legend) "bottom" else "none")

    if (!is.null(output_path)) {
        ggplot2::ggsave(
            filename = output_path,
            plot = p,
            ...
        )
    }

    return(p)
}
