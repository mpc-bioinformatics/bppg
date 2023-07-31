
#' Generates a table with characteristics for each subgraph in a list.
#'
#' @param S list of subgraphs, where peptide and protein nodes are collapsed
#' @param S2 list of subgraphs, where only protein nodes are collapsed
#' @param S3 list of subgraphs, where nodes are not collapsed
#' @param fastalevel Are the subgraphs on fasta level?
#' @param comparison name of comparison, for quantitative level (not in use currently)
#' @param file where to save the table.
#'
#'
#' @return table
#' @export
#'
#' @examples
#' # TODO
calculate_subgraph_characteristics <- function(S, #S2, S3,
                                               fastalevel = TRUE,
                                               #comparison = NULL,
                                               file = NULL) {

  Data <- NULL

  if (fastalevel) {
    comparisons <- 1
  } else {
    comparisons <- names(S)
  }

  for (j in 1:length(comparisons)) {

    if (fastalevel) {
      S_tmp <- S
     # S2_tmp <- S2
     # S3_tmp <- S3
    } else {
      S_tmp <- S[[j]]
     # S2_tmp <- S2[[j]]
    #  S3_tmp <- S3[[j]]
    }

    print(comparisons[j])

    #add progress bar to loop
    number_of_iterations <- length(S_tmp)
    pb <- pbapply::startpb(0, length(S_tmp))
    on.exit(pbapply::closepb(pb))

    for (i in 1:length(S_tmp)) {

      G_tmp <- S_tmp[[i]]
     # G2_tmp <- S2_tmp[[i]]
     # G3_tmp <- S3_tmp[[i]]


    #S_tmp <- igraph::as_incidence_matrix(G_tmp)

    nr_protein_nodes <- sum(igraph::V(G_tmp)$type)
    nr_peptide_nodes <- sum(!igraph::V(G_tmp)$type)
    nr_edges <- igraph::gsize(G_tmp)

    nr_edges_per_pep_node <- igraph::degree(G_tmp)[!igraph::V(G_tmp)$type]
    nr_unique_peptides <- sum(nr_edges_per_pep_node == 1)
    nr_shared_peptides <- sum(nr_edges_per_pep_node > 1)


    protein_acc <- igraph::V(G_tmp)$name[igraph::V(G_tmp)$type]
    protein_acc <- strsplit(protein_acc, ";")
    nr_protein_accessions <- sum(sapply(protein_acc, length))

    peptide_seq <- igraph::V(G_tmp)$name[!igraph::V(G_tmp)$type]
    peptide_seq <- strsplit(peptide_seq, ";")
    nr_peptide_sequences <- sum(sapply(peptide_seq, length))

    # nr_edges_per_pep_node2 <- igraph::degree(G2_tmp)[!igraph::V(G2_tmp)$type]
    # nr_peptide_sequences_unique <- sum(nr_edges_per_pep_node2 == 1)
    # nr_peptide_sequences_shared <- sum(nr_edges_per_pep_node2 > 1)

### TODO: add nr of unique and shared peptide sequences
### TODO: add infor about graph type (isomorphism list!)


    D_tmp <- data.frame(graph_ID = i,
                        nr_protein_nodes = nr_protein_nodes,
                        nr_peptide_nodes = nr_peptide_nodes,
                        nr_unique_peptide_nodes = nr_unique_peptides,
                        nr_shared_peptide_nodes = nr_shared_peptides,
                        nr_edges = as.integer(nr_edges),
                        nr_protein_accessions = nr_protein_accessions,
                        nr_peptide_sequences = as.integer(nr_peptide_sequences),
                        # nr_peptide_sequences_unique = nr_peptide_sequences_unique,
                        # nr_peptide_sequences_shared = nr_peptide_sequences_shared,
                        comparison = comparisons[j]

    )

    Data <- rbind(Data, D_tmp)

    pbapply::setpb(pb, i)
  }
  #progress bar command
  invisible(NULL)
}
  if (!is.null(file)) openxlsx::write.xlsx(Data, file, overwrite = TRUE)

  return(Data)
}
