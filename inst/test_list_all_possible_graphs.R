

library(igraph)
library(rje)

## m proteins, n peptides
## Spalten = Proteine, Zeilen = Peptide

################################################################################
### m = 1
M <- matrix(c(1), nrow = 1)
plot(graph_from_incidence_matrix(M), layout = igraph::layout_as_bipartite)

################################################################################
### m = 2

M <- matrix()


### Power set is set of possible nodes
x <- rje::powerSet(1:2)
x[[1]] <- NULL  ## delete empty set
x

poss_node_comb <- t(rje::powerSetMat(length(x)))  ## alle möglichen Kombinationen von Knoten

### es muss mindestens 2 Knoten geben, sonst wäre es m = 1

poss_node_comb <- poss_node_comb[,colSums(poss_node_comb) >=2]


G_list_tmp <- list()

for (i in 1:ncol(poss_node_comb)) {
  peptide_node_edges <- x[which(poss_node_comb[,i]==1)]

  ### TODO: Abfragen, dass alle Proteinknoten erwischt werden
  prot_nodes <- unique(unlist(peptide_node_edges))
  if(length(prot_nodes) < 2) next()

  M <- matrix(ncol = length(prot_nodes), nrow = length(peptide_node_edges))

  for(j in 1:length(peptide_node_edges)) {
    tmp <- prot_nodes %in% peptide_node_edges[[j]]
    M[j,] <- tmp
  }

  G <- igraph::graph_from_incidence_matrix(M)

  if (!igraph::is_connected(G)) next() ### nicht connected

  ### isomorphic zu einem der schon bestehenden Graphen?

  if(length(G_list_tmp)) {
    test_iso <- sapply(G_list_tmp, function(x) {
      return(isomorphic_bipartite(x, G))
    })
    ## falls Graph nicht isomorph zu irgendeinem vorherigen ist, ergänze Liste
    if(!any(test_iso)) {
      G_list_tmp <- rlist::list.append(G_list_tmp, G)
    }
  } else {
    G_list_tmp <- rlist::list.append(G_list_tmp, G)
  }
}


plotBipartiteGraph(G_list_tmp[[1]])
plotBipartiteGraph(G_list_tmp[[2]])


################################################################################
### m = 3

source("R/helpers-isomorphisms.R")
source("R/plotBipartiteGraph.R")

m = 3

library(devtools)
# install 'binaryLogic'
install_github("d4ndo/binaryLogic")
library(binaryLogic)

#library(R.utils)



### Power set is set of possible nodes
x <- rje::powerSet(1:m)
x[[1]] <- NULL  ## delete empty set
x

nr_poss_pep_nodes <- length(x)

nr_poss_node_comb <- 2^nr_poss_pep_nodes-1 # wieder Potenzmenge aller möglichen Knoten (außer leere Menge)
### leere Menge wird dann so codiert



#tmp <- binaryLogic::as.binary(2^2^6-1, n = 2^6)


#poss_node_comb <- t(rje::powerSetMat(length(x)))  ## alle möglichen Kombinationen von Knoten

### es muss mindestens log2(m+1) Knoten geben, sonst würden Proteine zusammengefasst

#poss_node_comb <- poss_node_comb[,colSums(poss_node_comb) >= log2(m+1)]


G_list_tmp <- list()

for (i in 1:nr_poss_node_comb) {
  print(i)

  ind <- binaryLogic::as.binary(i, n = nr_poss_pep_nodes)

  peptide_node_edges <- x[ind]

  ### TODO: Abfragen, dass alle Proteinknoten erwischt werden
  prot_nodes <- unique(unlist(peptide_node_edges))
  if(length(prot_nodes) < m) next()

  M <- matrix(ncol = length(prot_nodes), nrow = length(peptide_node_edges))

  for(j in 1:length(peptide_node_edges)) {
    tmp <- prot_nodes %in% peptide_node_edges[[j]]
    M[j,] <- tmp
  }

  G <- igraph::graph_from_incidence_matrix(M)

  if (!igraph::is_connected(G)) next() ### nicht connected

  if(any(duplicated(M, MARGIN = 2))) next()  ### doppelte Spalten, Protein-Knoten würden gemergt

  ### isomorphic zu einem der schon bestehenden Graphen?

  if(length(G_list_tmp)) {
    test_iso <- sapply(G_list_tmp, function(x) {
      return(isomorphic_bipartite(x, G))
    })
    ## falls Graph nicht isomorph zu irgendeinem vorherigen ist, ergänze Liste
    if(!any(test_iso)) {
      G_list_tmp <- rlist::list.append(G_list_tmp, G)
    }
  } else {
    G_list_tmp <- rlist::list.append(G_list_tmp, G)
  }
}





for (i in 1:length(G_list_tmp)) {
  plotBipartiteGraph(G_list_tmp[[i]])
}





