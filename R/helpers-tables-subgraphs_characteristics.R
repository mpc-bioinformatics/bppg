#
# G <- readRDS("C:\\Users\\schorkka\\UNI\\Promotion\\promotion_project\\data\\D4_without_isoforms\\D4_quant\\preprocessed\\subgraphs_quant_collprot_minAA7_mc2.rds")
# G2 <- readRDS("C:\\Users\\schorkka\\UNI\\Promotion\\promotion_project\\data\\D4_without_isoforms\\D4_quant\\preprocessed\\subgraphs_quant_collprot_collpept_minAA7_mc2.rds")
#
# G_collprot <- G[[1]]
# G_collprot_collpept <- G2[[1]]
# i = 1
# library(igraph)
# calculate_subgraph_characteristics(G_collprot)




#' Calculate table with characteristics for each subgraph
#'
#' @param G_collprot list of graphs
#' @param comparison is added as separate column
#'
#' @return table
#' @export
#'
#' @examples # TODO
calculate_subgraph_characteristics <- function(G_collprot, #G_collprot_collpept,
                                               #fastalevel = TRUE,
                                               comparison = NULL) {

  Data <- NULL

  total <- length(G_collprot)
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  for (i in 1:length(G_collprot)) {

    setTxtProgressBar(pb, i)
    #print(i)

    G_tmp <- G_collprot[[i]]
    #G_collpept_tmp <- G_collprot_collpept[[i]]

    ## protein nodes:
    names(V(G_tmp))[V(G_tmp)$type]

    # peptide sequences
    peptides_seq <- V(G_tmp)[!V(G_tmp)$type]
    ego_peptides_seq <- ego_size(G_tmp, mindist = 1, nodes = V(G_tmp)[!V(G_tmp)$type])
    is_unique_peptide_seq <- (ego_peptides_seq == 1)
    is_shared_peptide_seq <- (ego_peptides_seq > 1)

    # # pepide nodes
    # peptides_node <- V(G_collpept_tmp)[!V(G_collpept_tmp)$type]
    # ego_peptides_node <- ego_size(G_collpept_tmp, mindist = 1, nodes = V(G_collpept_tmp)[!V(G_collpept_tmp)$type])
    # is_unique_peptide_node <- (ego_peptides_node == 1)
    # is_shared_peptide_node <- (ego_peptides_node > 1)

    #
    M <- as_incidence_matrix(G_tmp) # rows = peptides, columns = proteins
    is_unique_peptide <- rowSums(M) == 1
    nr_unique_pep_node_per_prot <- colSums(M[is_unique_peptide,, drop = FALSE])
    nr_shared_pep_node_per_prot <- colSums(M[!is_unique_peptide,, drop = FALSE])





    D_tmp <- data.frame(ID = i,
                        Nr_prot_acc = length(unlist(strsplit(names(V(G_tmp))[V(G_tmp)$type], ";"))),
                        Nr_prot_node = sum(V(G_tmp)$type),
                        Nr_prot_node_only_unique_pep = sum(nr_unique_pep_node_per_prot > 0 & nr_shared_pep_node_per_prot == 0),
                        Nr_prot_node_uniq_and_shared_pep = sum(nr_unique_pep_node_per_prot > 0 & nr_shared_pep_node_per_prot > 0),
                        Nr_prot_node_only_shared_pep = sum(nr_unique_pep_node_per_prot == 0 & nr_shared_pep_node_per_prot > 0),
                        Nr_prot_node_with_unique_pep = sum(nr_unique_pep_node_per_prot > 0),


                        #has_multiple_prot = (dim_submatrix_tmp[2] > 1),#
                        #has_prot_without_unique_pep = any(nr_unique_pep_node_per_prot == 0 & nr_shared_pep_node_per_prot > 0),#

                        Nr_pep_seq = sum(!V(G_tmp)$type),
                        Nr_pep_seq_unique = sum(is_unique_peptide_seq),
                        Nr_pep_seq_shared = sum(!is_shared_peptide_seq),
                        # Nr_pep_node = sum(!V(G_collpept_tmp)$type),
                        # Nr_pep_node_unique = sum(is_unique_peptide_node),
                        # Nr_pep_node_shared = sum(is_shared_peptide_node),

                        Nr_edge = gsize(G_tmp)
    )

    if (!is.null(comparison)) {
      D_tmp <- data.frame(comparison = comparison, D_tmp)
    }

    Data <- rbind(Data, D_tmp)
  }
  close(pb)
  return(Data)
}
