
#' Table with information on each protein node
#'
#' @param G   \strong{list of list of igraph objects} \cr
#'            The graphs with collapsed protein and peptide nodes.
#'
#' @return A data frame with information on number of unique/shared peptides.
#' @export
#'
#' @seealso For the generation of the list of lists of igraphs: [generate_graphs_from_quant_data()]
#'
#' @examples # TODO

calculate_proteinnode_info <- function(G) {

  G2 <- lapply(G, function(x) {
    pbapply::pblapply(x, bppg::add_uniqueness_attributes)
  })

  accessions <- NULL
  comparison <- NULL
  graphID <- NULL
  ind_within_graph <- NULL
  nr_peptides <- NULL
  nr_unique_peptides <- NULL
  nr_shared_peptides <- NULL


  for (i in 1:length(G2)){
    for (j in 1:length(G2[[i]])) {

      G_tmp <- G2[[i]][[j]]
      ind_proteins <- which(igraph::V(G_tmp)$type)

      for (k in 1:length(ind_proteins)) {
        ind <- ind_proteins[k]

        accessions_tmp <- igraph::V(G_tmp)$name[ind]
        accessions <- c(accessions, accessions_tmp)
        comparison <- c(comparison, names(G2)[i])
        graphID <- c(graphID, j)
        ind_within_graph <- c(ind_within_graph, k)
        nr_unique_peptides_tmp <- igraph::V(G_tmp)$nr_unique_peptides[ind]
        nr_shared_peptides_tmp <- igraph::V(G_tmp)$nr_shared_peptides[ind]
        nr_peptides <- c(nr_peptides, nr_unique_peptides_tmp + nr_shared_peptides_tmp)
        nr_unique_peptides <- c(nr_unique_peptides, nr_unique_peptides_tmp)
        nr_shared_peptides <- c(nr_shared_peptides, nr_shared_peptides_tmp)

      }
    }
  }

  D <- data.frame(accessions = accessions, comparison = comparison, graphID = graphID,
                  ind_within_graph = ind_within_graph, nr_peptides = nr_peptides,
                  nr_unique_peptides = nr_unique_peptides, nr_shared_peptides = nr_shared_peptides)

  return(D)
}












#D3 <- read.xlsx("data/D3_without_isoforms/D3_quant/table_subgraph_characteristics_D3_quant_LTS50.xlsx")


#### table with protein information ####
# library(pbapply)
#
#
#
# human_proteins <- names(fasta1)
# Ecoli_proteins <- names(fasta2)
#
#
# is_spikein <- NULL
# true_ratio <- NULL
# taxonomy <- NULL
#
#
# for (i in 1:length(G2)){
#   for (j in 1:length(G2[[i]])) {
#
#     G_tmp <- G2[[i]][[j]]
#     ind_proteins <- which(V(G_tmp)$type)
#
#     for (k in 1:length(ind_proteins)) {
#       ind <- ind_proteins[k]
#
#       accessions_tmp <- V(G_tmp)$name[ind]
#       accessions <- c(accessions, accessions_tmp)
#       comparison <- c(comparison, names(G2)[i])
#       graphID <- c(graphID, j)
#       ind_within_graph <- c(ind_within_graph, k)
#       nr_unique_peptides_tmp <- V(G_tmp)$nr_unique_peptides[ind]
#       nr_shared_peptides_tmp <- V(G_tmp)$nr_shared_peptides[ind]
#       nr_peptides <- c(nr_peptides, nr_unique_peptides_tmp + nr_shared_peptides_tmp)
#       nr_unique_peptides <- c(nr_unique_peptides, nr_unique_peptides_tmp)
#       nr_shared_peptides <- c(nr_shared_peptides, nr_shared_peptides_tmp)
#
#       ## is the protein a spike in (contains UPS in name)
#       accessions_tmp_split <- limma::strsplit2(accessions_tmp, ";")[1,]
#       taxonomy_tmp <- NULL
#       for (l in 1:length(accessions_tmp_split)) {
#         if (accessions_tmp_split[l] %in% human_proteins) {
#           taxonomy_tmp <- c(taxonomy_tmp, "human")
#         } else {
#           if (accessions_tmp_split[l] %in% Ecoli_proteins) {
#             taxonomy_tmp <- c(taxonomy_tmp, "Ecoli")
#           } else {
#             taxonomy_tmp <- c(taxonomy_tmp, "contaminant")
#           }
#         }
#       }
#       if (all(taxonomy_tmp == "human")) {
#         taxonomy <- c(taxonomy, "human")
#         is_spikein <- c(is_spikein, FALSE)
#         true_ratio <- c(true_ratio, 1)
#         next
#       }
#       if (all(taxonomy_tmp == "Ecoli")) {
#         taxonomy <- c(taxonomy, "Ecoli")
#         is_spikein <- c(is_spikein, TRUE)
#         true_ratio <- c(true_ratio, 1/3)
#         next
#       }
#       if (any(taxonomy_tmp == "contaminant")) {
#         taxonomy <- c(taxonomy, "contaminant")
#         is_spikein <- c(is_spikein, NA)
#         true_ratio <- c(true_ratio, NA)
#         next
#       }
#
#       if (any(taxonomy_tmp == "human") & any(taxonomy_tmp == "Ecoli")) {
#         taxonomy <- c(taxonomy, "mixed")
#         is_spikein <- c(is_spikein, NA)
#         true_ratio <- c(true_ratio, NA)
#         next
#       }
#     }
#   }
# }
#
# D <- data.frame(accessions = accessions, comparison = comparison, graphID = graphID,
#                 ind_within_graph = ind_within_graph, nr_peptides = nr_peptides,
#                 nr_unique_peptides = nr_unique_peptides, nr_shared_peptides = nr_shared_peptides,
#                 taxonomy = taxonomy,
#                 is_spikein = is_spikein,
#                 true_ratio = true_ratio)
