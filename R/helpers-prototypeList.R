#' Generates a list of graph prototypes for the different isomorphism classes and their occurence.
#'
#'
#' @param G                  \strong{igraph graph object} \cr
#'                           A graph.
#' @param sort_by_nr_edges   \strong{logical} \cr
#'                           If \code{TRUE}, the list of prototypes is sorted by number of edges.
#'
#' @return A list of prototype graphs plus their count.
#' @export
#'
#' @examples # TODO
#'

generatePrototypeList <- function(G, sort_by_nr_edges = FALSE) {


  counter <- integer(length(G))
  k <- 1

  pb <- pbapply::startpb(min = 0, max = 1)


  i <- 1

  # go trough list of graphs
  while(i <= length(G)){

    ## if end of list is reached:
    if (i == length(G)) {
      counter[i] <- 1
      i <- i + 1
      next
    }


    G_tmp <- G[[i]]


    # Which graphs are isomorphic to G_tmp?
    x <- sapply(G[(i+1):length(G)], function(x) {
      isomorphic_bipartite(x, G_tmp)
    })
    ind <- which(x)


    # delete Graphs isomorphic to G_tmp graphs (-> list becomes smaller)
    # G_tmp itself is a new isomorphism class.
    if (length(ind) > 0) {
      G <- G[-(ind+i)]
      counter[i] <- length(ind) + 1
      counter <- counter[-(ind+i)]
    } else {
      counter[i] <- 1
    }

    pbapply::setpb(pb, i/length(G))
    i <- i + 1
  }


  # sort list of prototypes according to number of edges
  if (sort_by_nr_edges) {
    nr_edges <- sapply(G, igraph::gsize)
    ord <- order(nr_edges)

    G <- G[ord]
    counter <- counter[ord]
  }


  pbapply::closepb(pb)
  return(list(graphs = G, counter = counter))
}

