#' Adds vertex attributes with uniqueness of peptides and number of unique peptides for proteins.
#'
#' @param G \strong{igraph graph object} \cr
#'          A peptide-protein graph.
#'
#' @return A graph with 2 additional vertex attributes, uniqueness and nr_unique_peptides
#' @export
#'
#' @seealso [generate_graphs_from_FASTA()], [generate_quant_graphs()], [add_average_pep_ratio()]
#'
#' @examples

add_uniqueness_attributes <- function(G) {

  ### FALSE = peptide, TRUE = protein
  igraph::V(G)$type

  uniqueness <- igraph::degree(G, igraph::V(G)) == 1
  uniqueness[igraph::V(G)$type] <- NA # attribute only for peptides

  G <- igraph::set_vertex_attr(G, "uniqueness", value = uniqueness)

  unique_peptide_nodes <- igraph::V(G)[igraph::V(G)$uniqueness & !is.na(igraph::V(G)$uniqueness)]
  shared_peptide_nodes <- igraph::V(G)[!igraph::V(G)$uniqueness & !is.na(igraph::V(G)$uniqueness)]

  neighborhood <-  igraph::ego(G, order = 1, mindist = 1, nodes = igraph::V(G))

  nr_unique_peptides <- sapply(neighborhood, function(x) sum(x%in% unique_peptide_nodes))
  nr_unique_peptides[!igraph::V(G)$type] <- NA ## attribute only for proteins
  G <- igraph::set_vertex_attr(G, "nr_unique_peptides", value = nr_unique_peptides)

  nr_shared_peptides <- sapply(neighborhood, function(x) sum(x%in% shared_peptide_nodes))
  nr_shared_peptides[!igraph::V(G)$type] <- NA ## attribute only for proteins
  G <- igraph::set_vertex_attr(G, "nr_shared_peptides", value = nr_shared_peptides)

}






#' Adds average peptide ratios as a attribute to the graphs, if a list of peptide ratios is already present.
#'
#' @param G      \strong{igraph graph object} \cr
#'               A peptide-protein graph.
#' @param type   \strong{character} \cr
#'               !NOT USED AT THE MOMENT!
#'
#' @return A graph with added peptide ratio attributes.
#' @export
#'
#' @seealso [generate_graphs_from_FASTA()], [generate_quant_graphs()], [add_uniqueness_attributes()]
#'
#' @examples

add_average_pep_ratio <- function(G, type = "geom_mean") {

  pep_ratio <- igraph::V(G)$pep_ratio
  pep_ratio_split <- strsplit(pep_ratio, ";")

  pep_ratio_aggr <- sapply(pep_ratio_split, function(x){
    bppg::geom_mean(as.numeric(x))})


  nr_sequences <- sapply(pep_ratio_split, length)

  G <- igraph::set_vertex_attr(G, "pep_ratio_aggr", value = pep_ratio_aggr)
  G <- igraph::set_vertex_attr(G, "nr_sequences", value = nr_sequences)
  return(G)
}



