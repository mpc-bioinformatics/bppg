#
# G <- readRDS("../promotion_project_new/data/subgraphs_collprotpept_minAA7_mc2.rds")
#
# library(igraph)
#
#
# x <- sapply(G, gorder)
#
#
# G <- G[order(x)]






# X <- matrix(rep(FALSE, length(G)^2), nrow = length(G), ncol = length(G))
#
# for (i in 1:length(G)) {
#   print(i)
#   x <- sapply(G, function(x) {
#     isomorphic(x, G[[i]])
#   })
#
#   X[,i] <- x
# }
#








#' Generates a list of graph prototypes for the different isomorphism classes and
#'
#'
#' @param G bla
#'
#' @return bla
#' @export
#'
#' @examples # TODO
generatePrototypeList <- function(G, sort_by_nr_edges = FALSE) {





  #prototype_list <- G
  counter <- integer(length(G))
  k <- 1


  pb <- pbapply::startpb(min = 0, max = 1)

  # for (i in 1:length(G))

  i <- 1

  while(i < length(G)){
    #print(length(G))
    G_tmp <- G[[i]]
    #if (k == 1) {prototype_list[[k]] <- G_tmp; k <- k + 1; next}


    x <- sapply(G[(i+1):length(G)], function(x) {
      isomorphic_bipartite(x, G_tmp)
    })
    ind <- which(x)
    #ind <- ind[-1]  # entferne den G_tmp selbst (ist zu sich selbst isomorph)

    # delete duplicated graphs:
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

  #   for (j in 2:length(G)) {
  #
  #     #G1 <<- G_tmp
  #     #G2 <<- prototype_list[[j]]
  #     #print(prototype_list[[j]])
  #
  #     iso <- bppg::isomorphic_bipartite(G_tmp, prototype_list[[j]])
  #
  #     if (iso) {
  #       prototype_list[[j]] <- NULL
  #       k <- k + 1;
  #       next
  #     }
  #   }
  #
  # }


  if (sort_by_nr_edges) {
    nr_edges <- sapply(G, igraph::gsize)
    ord <- order(nr_edges)

    G <- G[ord]
    counter <- counter[ord]
  }


  pbapply::closepb(pb)
  return(list(G, counter))
}

# P <- generatePrototypeList(G)
#
# saveRDS(P, file = "Graph_prototypes_D4.rds")

### TODO: plotten und einzeln als graphml file abspeichern!

