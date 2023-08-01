


#' Generates a list of graph prototypes for the different isomorphism classes and
#'
#'
#' @param G bla
#'
#' @return bla
#' @export
#'
#' @examples # TODO
#'
#'
generatePrototypeList <- function(G, sort_by_nr_edges = FALSE) {


  counter <- integer(length(G))
  k <- 1

  pb <- pbapply::startpb(min = 0, max = 1)

  # for (i in 1:length(G))

  i <- 1

  # Gehe Liste an Graphen durch
  while(i < length(G)){

    G_tmp <- G[[i]]


    # Welche Graphen sind isomorphic zu G_tmp?
    x <- sapply(G[(i+1):length(G)], function(x) {
      isomorphic_bipartite(x, G_tmp)
    })
    ind <- which(x)


    # delete Graphs isomorphic to G_tmp graphs (-> list becomes smallet)
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

# P <- generatePrototypeList(G)
#
# saveRDS(P, file = "Graph_prototypes_D4.rds")

### TODO: plotten und einzeln als graphml file abspeichern!

### WICHTIG: Laufzeit wird zu Beginn stark überschätzt, geht aber relativ schnell

