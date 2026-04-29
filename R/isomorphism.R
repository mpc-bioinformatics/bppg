# Functions in this file:
# .directBipartiteGraph
# .isomorphicBipartite
# .generatePrototypeList
#

#' Transform a bipartite graph into a directed graph.
#'
#' @param bip_graph   \strong{graph (igraph)} \cr
#'                    A bipartite graph.
#' @param from_type   \strong{logical} \cr
#'                    If \code{TRUE}, the edges will go out from the vertices
#'                    with the type \code{TRUE} from the bipartite graph.
#'
#' @return A bipartite graph that is know directed.
#'
#'
#' @importFrom igraph %->% as_directed E reverse_edges V
#'
.directBipartiteGraph <- function(bip_graph, from_type = FALSE) {

    ## turn undirected into directed edges
    bip_graph <- igraph::as_directed(bip_graph, mode = "arbitrary")

    from_vs <- igraph::V(bip_graph)[igraph::V(bip_graph)$type == from_type]
    to_vs <- igraph::V(bip_graph)[igraph::V(bip_graph)$type == !from_type]

    ## reverse edges going from the "to-group" to the "from-group"
    bip_graph <- igraph::reverse_edges(bip_graph,
        igraph::E(bip_graph)[to_vs %->% from_vs])

    return(bip_graph)
}

#' Enchanced version of the igraph::isomorphic function that also considers the
#' node type in bipartite graphs, e.g. that W- and M-shaped graphs are NOT
#' isomorphic
#'
#' @param graph1   \strong{graph (igraph)} \cr
#'                 First graph.
#' @param graph2   \strong{graph (igraph)} \cr
#'                 Second graph.
#' @param ...      currently unused
#'
#' @return TRUE if graphs are isomorphic, FALSE if not.
#'
#'
#' @seealso [generateGraphsFromEdgelist()]
#'
#' @examples
#'
#' M1 <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
#' G1 <- igraph::graph_from_biadjacency_matrix(M1)
#'
#' M2 <- matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE)
#' G2 <- igraph::graph_from_biadjacency_matrix(M2)
#'
#' .isomorphicBipartite(G1, G2)
#'
#' @importFrom igraph is_directed isomorphic

.isomorphicBipartite <- function(graph1, graph2, ...) {

    ## direct graphs if they are not directed yet
    if (!igraph::is_directed(graph1))   graph1 <- .directBipartiteGraph(graph1)
    if (!igraph::is_directed(graph2))   graph2 <- .directBipartiteGraph(graph2)

    igraph::isomorphic(graph1, graph2, method = "vf2")
}


#' Generates a list of graph prototypes for the different isomorphism classes
#' and their occurence.
#'
#'
#' @param G                  \strong{list of igraph graph objects} \cr
#'                           A list of graphs.
#' @param sort_by_nr_edges   \strong{logical} \cr
#'                           If \code{TRUE}, the list of prototypes is sorted by
#'                           number of edges.
#'
#' @return A list of prototype graphs plus their count.
#'
#'
#' @examples
#' file <- system.file("extdata", "quantGraphs_collpept.rds", package = "bppg")
#' graphs <- readRDS(file)
#' graphs <- unlist(graphs, recursive = FALSE)
#'
#' PL <- bppg:::.generatePrototypeList(graphs, sort_by_nr_edges = FALSE)
#' # prototype list contains 3 isomorphism classes (I, N and M shaped graphs)
#'
#' @importFrom pbapply setpb startpb closepb
#' @importFrom igraph gsize

.generatePrototypeList <- function(G, sort_by_nr_edges = FALSE) {
    counter <- integer(length(G))
    pb <- pbapply::startpb(min = 0, max = 1)
    i <- 1
    # While-loop over G to find isomorphism classes. All graphs belonging to the
    # class found in the current iteration are removed, so the length of G is
    # decreasing.
    while (i <= length(G)) {
        ## if end of list is reached:
        if (i == length(G)) {
            counter[i] <- 1
            i <- i + 1
            next
        }
        G_tmp <- G[[i]]
        ## Which graphs are isomorphic to G_tmp?
        x <- vapply(G[(i + 1):length(G)], function(x) {
            .isomorphicBipartite(x, G_tmp)
        }, logical(1))
        ind <- which(x)
        ## delete Graphs isomorphic to G_tmp graphs (-> list becomes smaller)
        ## G_tmp itself is a new isomorphism class.
        if (length(ind) > 0) {
            G <- G[-(ind + i)]
            counter[i] <- length(ind) + 1
            counter <- counter[-(ind + i)]
        } else {
            counter[i] <- 1
        }
        pbapply::setpb(pb, i / length(G))
        i <- i + 1
    }
    ## sort list of prototypes according to number of edges
    if (sort_by_nr_edges) {
        nr_edges <- vapply(G, igraph::gsize, FUN.VALUE = numeric(1))
        ord <- order(nr_edges)
        G <- G[ord]
        counter <- counter[ord]
    }
    pbapply::closepb(pb)
    return(list(graphs = G, counter = counter))
}
