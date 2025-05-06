#' Generates a table with characteristics for each subgraph in a list.
#'
#' @param S            \strong{list of igraph graph objects} \cr
#'                     A list of subgraphs, where peptide and protein nodes are collapsed.
#' @param fastalevel   \strong{logical} \cr
#'                     If \code{TRUE}, the subgraphs should be on fasta level
#' @param prototype    \strong{logical} \cr
#'                     If \code{TRUE}, the subgraphs should be part of a prototype list
#' @param file         \strong{character} \cr
#'                     A file path where to save the table.
#'
#'
#' @return A table with the characteristics.
#' @export
#'
#' @examples
#'

calculate_subgraph_characteristics <- function(S, #S2, S3,
                                               fastalevel = TRUE,
                                               prototype = FALSE,
                                               #comparison = NULL,
                                               file = NULL) {

  if (prototype) {
    counter <- S$counter
    S <- S$graph
  }

  Data <- NULL


  ### TODO: das kann man auch anders lösen, indem man guckt ob es ne liste ist? Dann würde das Argument wegfallen
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



      unique_peptide_nodes <- igraph::V(G_tmp)[(igraph::degree(G_tmp) == 1 & !igraph::V(G_tmp)$type)]

      if (length(unique_peptide_nodes) == 0) { # Fall: keine uniquen Peptide im ganzen Graphen
        nr_prot_node_only_unique_pep <- 0
        nr_prot_node_unique_and_shared_pep <- 0
        nr_prot_node_only_shared_pep <-  nr_protein_nodes

      } else {
        if (length(unique_peptide_nodes) == 1 & nr_protein_nodes == 1) { # Fall: I-shaped graph
          nr_prot_node_only_unique_pep <- 1
          nr_prot_node_unique_and_shared_pep <- 0
          nr_prot_node_only_shared_pep <-  0
        } else {

          # neighborhood of the unique peptides (these are proteins with a unique peptide)
          NH_of_unique_peptides <- igraph::ego(G_tmp, order = 1, mindist = 1, nodes = unique_peptide_nodes)

          nr_prot_node_only_unique_pep <- 0
          nr_prot_node_unique_and_shared_pep <- length(NH_of_unique_peptides) # = Anzahl uniquer Peptide??
          nr_prot_node_only_shared_pep <-  nr_protein_nodes - nr_prot_node_unique_and_shared_pep#length(NH_of_unique_peptides)
        }
      }

      ### TODO: add nr of unique and shared peptide sequences
      ### TODO: add info about graph type (isomorphism list!)

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
                          nr_prot_node_only_unique_pep = nr_prot_node_only_unique_pep,
                          nr_prot_node_unique_and_shared_pep = nr_prot_node_unique_and_shared_pep,
                          nr_prot_node_only_shared_pep = nr_prot_node_only_shared_pep,
                          comparison = comparisons[j]

      )

      Data <- rbind(Data, D_tmp)

      pbapply::setpb(pb, i)
    }
    #progress bar command
    invisible(NULL)
  }

  if (prototype) {
    Data <- cbind(Data, counter = counter)
  }

  if (!is.null(file)) openxlsx::write.xlsx(Data, file, overwrite = TRUE)

  return(Data)
}
