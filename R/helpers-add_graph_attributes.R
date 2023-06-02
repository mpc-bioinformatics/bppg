


#' Adds vertex attributes with uniqueness of peptides and number of unique peptides
#' for proteins
#'
#' @param G graph
#'
#' @return graph with 2 additional vertex attributes, uniqueness and nr_unique_peptides
#' @export
#'
#' @examples # TODO
add_uniqueness_attributes <- function(G) {

  ### FALSE = peptide, TRUE = protein
  igraph::V(G)$type

  uniqueness <- degree(G, V(G)) == 1
  uniqueness[V(G)$type] <- NA # attribute only for peptides

  G <- set_vertex_attr(G, "uniqueness", value = uniqueness)

  unique_peptide_nodes <- V(G)[V(G)$uniqueness & !is.na(V(G)$uniqueness)]
  shared_peptide_nodes <- V(G)[!V(G)$uniqueness & !is.na(V(G)$uniqueness)]

  neighborhood <-  ego(G, order = 1, mindist = 1, nodes = V(G))

  nr_unique_peptides <- sapply(neighborhood, function(x) sum(x%in% unique_peptide_nodes))
  nr_unique_peptides[!V(G)$type] <- NA ## attribute only for proteins
  G <- set_vertex_attr(G, "nr_unique_peptides", value = nr_unique_peptides)

  nr_shared_peptides <- sapply(neighborhood, function(x) sum(x%in% shared_peptide_nodes))
  nr_shared_peptides[!V(G)$type] <- NA ## attribute only for proteins
  G <- set_vertex_attr(G, "nr_shared_peptides", value = nr_shared_peptides)

}
