# Functions in this file:
# myDiamond
# setNodeLabels
# .plotBipartiteGraph

 library(igraph)
 library(ggraph)
# plotBPPG <- function(G) {
#     ## add Node labels for three categories
#     V(G)$type <- !V(G)$type
#     type <- integer(length(igraph::V(G)))
#     type[!igraph::V(G)$type] <- "Protein"
#     type[igraph::V(G)$type] <- "Shared Peptides"
#     type[igraph::V(G)$type & igraph::degree(G) == 1] <- "Unique Peptides"

#     #set edge attributes for edges to unique peptides
#     uniquePeptides <- V(G)[igraph::V(G)$type & igraph::degree(G) == 1]
#     E(G)$unique <- FALSE
#     E(G)[.from(uniquePeptides)]$unique <- TRUE

#     layout <- create_layout(G, layout = 'igraph', algorithm = 'bipartite')
#     layout$nodeType <- type
#     # layout$name[layout$type] <- round(V(G)$pep_logRatio [layout$type], 2)
#     layout$name[layout$type] <- "peptide"

#     layout <- layout[order(layout$x),]
#     layout$x[!layout$type] <- (seq(0, max(layout$x), length.out = sum(!layout$type)) + layout$x[!layout$type])/2
#     layout <- layout[order(layout$.ggraph.index),]

#     ggraph(G, layout = layout) + 
#         geom_edge_link0(aes(edge_alpha = factor(unique))) + 
#         geom_node_point(size = 15, aes(fill = factor(nodeType), shape = factor(nodeType))) +
#         geom_node_label(aes(label = name, fill = factor(nodeType)), size = 4, show.legend = FALSE)  +
#         scale_edge_alpha_manual(name = "Edge to unique peptide", values = c(0.3, 1)) + 
#         scale_shape_manual(values = c(21,22,23), name = "Node Type") + 
#         scale_fill_manual(values = c("mediumseagreen", "cadetblue2", "#FF8C00"), name = "Node Type") +
#         coord_cartesian(ylim=c(-0.1,1.1)) + theme_graph()

# }

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
        #### TODO: was ist, wenn es mehr als 26 Proteine gibt?
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
    if (node_labels_peptides == "pep_ratios") {
        pep_ratios <- igraph::V(G)$pep_ratio
        names_G[Layout[, 2] == 0] <- round(pep_ratios[Layout[, 2] == 0],
            round_digits)
    }
    if (node_labels_peptides == "pep_ratio_aggr") {
        pep_ratios <- igraph::V(G)$pep_ratio_aggr
        names_G[Layout[, 2] == 0] <- round(pep_ratios[Layout[, 2] == 0],
            round_digits)
    }
    if (node_labels_peptides == "") {
        names_G[Layout[, 2] == 0] <- NA
    }

    igraph::set_vertex_attr(G, name = "name", value = names_G)
}

#' Plot a bipartite graph of peptides and proteins.
#' 
#' @param G                         \strong{igraph graph object} \cr
#'                                  A bipartite peptide-protein graph.
#' @param legend                    \strong{logical} \cr
#'                                  Whether to display the legend.
#' @param vertex.color              \strong{character vector} \cr
#'                                  Colors for the protein, shared peptide, and unique peptide nodes, respectively.
#' @param vertex.size               \strong{numeric} \cr
#'                                  Size of the nodes.
#' @param vertex.label.cex          \strong{numeric} \cr
#'                                  Size of the node labels.        
#' @param edge.width                \strong{numeric} \cr        
#'                                 Width of the edges.  
#' @param useCanonicalPermutation   \strong{logical} \cr
#'                                 Whether to use the canonical permutation of the graph.
#' @param three_shapes               \strong{logical} \cr
#'                                Whether to use three different shapes for the nodes (protein, shared peptide, unique peptide).
#' @param node_labels_proteins      \strong{character} \cr
#'                                 The type of labels for the proteins. Options
#'                                are "letters" or "accessions".
#' @param node_labels_peptides      \strong{character} \cr
#'                                The type of labels for the peptides. Options
#'                               are "numbers", "pep_ratios", or "pep_ratio_aggr".
#'@param round_digits              \strong{integer} \cr
#'                                The number of digits to round the peptide ratios to.
#' @param save_path                 \strong{character} \cr
#'                                Path to save the plot as a JPEG file. If NULL, the plot is not saved.
#' @return A ggplot object representing the bipartite graph.     


.plotBipartiteGraph <- function(
    G,
    legend = TRUE,
    vertex.color = c("mediumseagreen", "cadetblue2", "coral1"),
    vertex.size = 15,
    vertex.label.cex = 1,
    edge.width = 1,
    useCanonicalPermutation = FALSE,
    three_shapes = FALSE,
    node_labels_proteins = "letters",
    node_labels_peptides = "numbers",
    round_digits = 2,
    save_path = NULL
) {
    igraph::V(G)$type <- !igraph::V(G)$type

    if (useCanonicalPermutation) {
        permutation <- igraph::canonical_permutation(G)
        G <- igraph::permute(G, permutation$labeling)
    }
    G <- .setNodeLabels(
        G,
        node_labels_peptides = node_labels_peptides,
        node_labels_proteins = node_labels_proteins,
        round_digits = round_digits
    )
    degree <- igraph::degree(G)

    igraph::V(G)$node_type <- ifelse(
        !igraph::V(G)$type,"protein",
        ifelse(
            degree == 1, "unique peptide", "shared peptide")
    )

    if (!three_shapes) {
        igraph::V(G)$node_shape_type <- ifelse(
            igraph::V(G)$node_type == "protein", "protein", "peptide")

        shape_values <- c(
            protein = 21, peptide = 22)
    } else {
        igraph::V(G)$node_shape_type <- igraph::V(G)$node_type
        shape_values <- c(protein = 21, `shared peptide` = 22,`unique peptide` = 23)
    }

    p <- ggraph::ggraph(G, layout = "bipartite") +
        ggraph::geom_edge_link(linewidth = edge.width) +
        ggraph::geom_node_point(
            ggplot2::aes(shape = node_shape_type, fill = node_type), size = vertex.size) +
        ggraph::geom_node_text(
            ggplot2::aes(label = name), size = vertex.label.cex, family = "sans") +
        ggplot2::scale_shape_manual(values = shape_values) +
        ggplot2::scale_fill_manual(values = c(protein = vertex.color[1], `shared peptide` = vertex.color[2],`unique peptide` = vertex.color[3])) +
        ggraph::theme_graph(base_family = "sans") +
        ggplot2::theme(legend.position = if (legend) "bottom" else "none")

    if (!is.null(save_path)) {
        ggplot2::ggsave(
            filename = save_path,
            plot = p,
            width = 10,
            height = 7,
            dpi = 300,
            device = "jpeg"
        )
    }

    return(p)
}