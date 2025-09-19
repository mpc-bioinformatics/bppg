
#' Table with information on each protein node
#'
#' @param G   \strong{list of list of igraph objects} \cr
#'            The graphs with collapsed protein and peptide nodes.
#'
#' @return A data frame with information on number of unique/shared peptides.
#'
#'
#' @seealso For the generation of the list of lists of igraphs: [generateGraphsFromQuantData()]
#'
#' @examples ## TODO

.calculateProteinNodeInfo <- function(G) {
    G2 <- lapply(G, function(x) {
        pbapply::pblapply(x, .addUniquenessAttributes)
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
                nr_unique_pept_tmp <- igraph::V(G_tmp)$nr_unique_peptides[ind]
                nr_shared_pept_tmp <- igraph::V(G_tmp)$nr_shared_peptides[ind]
                nr_peptides <- c(nr_peptides, 
                    nr_unique_pept_tmp + nr_shared_pept_tmp)
                nr_unique_peptides <- c(nr_unique_peptides, nr_unique_pept_tmp)
                nr_shared_peptides <- c(nr_shared_peptides, nr_shared_pept_tmp)
            }
        }
    }

    D <- data.frame(accessions = accessions, comparison = comparison, 
        graphID = graphID, ind_within_graph = ind_within_graph,
        nr_peptides = nr_peptides, nr_unique_peptides = nr_unique_peptides,
        nr_shared_peptides = nr_shared_peptides)
    return(D)
}
