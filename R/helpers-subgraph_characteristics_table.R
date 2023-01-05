
### Calculate subgraph characteristics table

### TODO: irgendwas stimmt hier nicht? ZUmindest bei den Prototypen läuft hier was falsch!


calculate_subgraph_characteristics <- function(S, fastalevel = TRUE, comparison = NULL) {

  require(igraph)
  require(pbapply)

  Data <- NULL

  #add progress bar to loop
  number_of_iterations <- length(S)
  pb <- startpb(0, length(S))
  on.exit(closepb(pb))



  for (i in 1:length(S)) {

    G_tmp <- S[[i]]
    S_tmp <- igraph::as_incidence_matrix(G_tmp)


    nr_proteins <- sum(igraph::V(G_tmp)$type)
    nr_peptides <- sum(!igraph::V(G_tmp)$type)
    nr_edges <- igraph::gsize(G_tmp)

    nr_edges_per_pep_node <- degree(G_tmp)[!igraph::V(G_tmp)$type]
    nr_unique_peptides <- sum(nr_edges_per_pep_node == 1)
    nr_shared_peptides <- sum(nr_edges_per_pep_node > 1)


    D_tmp <- data.frame(ID = i,
                        nr_protein_nodes = nr_proteins,
                        nr_peptide_nodes = nr_peptides,
                        nr_unique_peptide_nodes = nr_unique_peptides,
                        nr_shared_peptide_nodes = nr_shared_peptides,
                        nr_edges = nr_edges
    )

    Data <- rbind(Data, D_tmp)

    setpb(pb, i)
  }

  #progress bar command
  invisible(NULL)

  return(Data)
}
