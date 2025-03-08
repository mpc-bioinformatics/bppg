
#' Calculate number for table
#' @param subgraph_char_tab subgraphs characteristics tables
#' @param isomorph_list isomorph lists
#'
#' @return bla
#' @export
#'
#' @examples # TODO
calculate_summary_table <- function(subgraph_char_tab, isomorph_list) {

  D <- subgraph_char_tab

  ind_zero_pep <- which(D$nr_peptide_nodes == 0)
  if (length(ind_zero_pep) > 0) D <- D[-ind_zero_pep, ]

  ### largest graph (in terms of number of protein nodes)
  ind_largest <- which.max(D$nr_protein_nodes)

  ### 2nd largest graph (in terms of number of protein nodes)
  ind_largest2 <- which(D$nr_protein_nodes == sort(D$nr_protein_nodes, decreasing = TRUE)[2])

  c(
    sum(D$nr_protein_accessions),                 ### Nr of proteins
    sum(D$nr_protein_nodes),                ### Nr of protein nodes
    sum(D$nr_peptide_sequences),                  ### Nr of peptides
    sum(D$nr_peptide_nodes),                 ### Nr of peptide nodes
    sum(D$nr_edges),                     ### Nr of edges

    as.integer(nrow(D)),                ### Nr of graphs
    as.integer(sum(D$nr_protein_nodes == 1 & D$nr_peptide_nodes == 1)),                 ### Nr of graphs with only 1 protein node
    as.integer(length(isomorph_list$isomorph_list) - 1*(length(ind_zero_pep) > 0)), ### Nr of isomorphism classes

    D$Nr_prot_node[ind_largest],        ### largest system
    D$Nr_pep_node[ind_largest],
    D$Nr_edge[ind_largest],

    D$Nr_prot_node[ind_largest2],       ### 2nd largest system
    D$Nr_pep_node[ind_largest2],
    D$Nr_edge[ind_largest2]

  )

}





### TODO: Funktion, die dann die Vergleichstabelle ausrechnet
