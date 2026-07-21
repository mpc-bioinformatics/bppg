# Functions in this file:
# .calculateProteinNodeInfo
# .calculateSubgraphCharcateristics


#' Table with information on each protein node
#'
#' @param G   \strong{list of list of igraph objects} \cr
#'            The graphs with collapsed protein and peptide nodes.
#' @param verbose     \strong{logical} \cr
#'                    If \code{TRUE}, additional information on each iteration
#'                    of the optimization is printed
#'
#' @return A data frame with information on number of unique/shared peptides.
#'
#'
#' @seealso For the generation of the list of lists of igraphs:
#'  [generateGraphsFromQuantData()]
#'
#' @examples
#'
#' file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
#' graphs <- readRDS(file)
#'
#' ProtInfo <- bppg:::.calculateProteinNodeInfo(G = graphs, verbose = FALSE)
#'
#'
#' @importFrom pbapply pblapply pboptions
#' @importFrom igraph V

.calculateProteinNodeInfo <- function(G, verbose = FALSE) {
    if (!verbose) {
        pbo <- pbapply::pboptions(type = "none")
        on.exit(pbapply::pboptions(pbo), add = TRUE)
    }
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

    for (i in seq_along(G2)){
        for (j in seq_along(G2[[i]])) {
            G_tmp <- G2[[i]][[j]]
            ind_proteins <- which(igraph::V(G_tmp)$type)

            for (k in seq_along(ind_proteins)) {
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


#' Generates a table with characteristics for each subgraph in a list.
#'
#' @param S            \strong{list of igraph graph objects} \cr
#'                     A list of subgraphs, where peptide and protein nodes
#'                     are collapsed.
#' @param fastalevel   \strong{logical} \cr
#'                     If \code{TRUE}, the subgraphs should be on fasta level
#' @param prototype    \strong{logical} \cr
#'                     If \code{TRUE}, the subgraphs should be part of a
#'                     prototype list
#' @param file         \strong{character} \cr
#'                     A file path where to save the table.
#' @param verbose     \strong{logical} \cr
#'                    If \code{TRUE}, additional information on each iteration
#'                    of the optimization is printed
#'
#' @return A table with the characteristics.
#'
#'
#' @examples # TODO
#'
#' @importFrom pbapply closepb pboptions startpb setpb
#' @importFrom igraph gsize V
#' @importFrom openxlsx write.xlsx
## TODO was ist mit dem alten S2? ZU LANG
## TODO: enthält noch for-Schleifen
.calculateSubgraphCharacteristics <- function(S, #S2, S3,
    fastalevel = TRUE,
    prototype = FALSE,
    #comparison = NULL,
    file = NULL,
    verbose = FALSE) {
    if (!verbose) {
        pbo <- pbapply::pboptions(type = "none")
        on.exit(pbapply::pboptions(pbo), add = TRUE)
    }

    if (prototype) {
        counter <- S$counter
        S <- S$graph
    }

    Data <- NULL  ## TODO Allocation!

    ## TODO: das kann man auch anders lösen,
    ## indem man guckt ob es ne liste ist? Dann würde das Argument wegfallen
    if (fastalevel) {
        comparisons <- 1
    } else {
        comparisons <- names(S)
    }

    for (j in seq_along(comparisons)) {
        if (fastalevel) {
            S_tmp <- S
            ## S2_tmp <- S2
            ## S3_tmp <- S3
        } else {
            S_tmp <- S[[j]]
            ## S2_tmp <- S2[[j]]
            ## S3_tmp <- S3[[j]]
        }

        if (verbose) message(comparisons[j])

        ## add progress bar to loop
        pb <- pbapply::startpb(0, length(S_tmp))
        on.exit(pbapply::closepb(pb))

        for (i in seq_along(S_tmp)) {
            G_tmp <- S_tmp[[i]]

            nr_protein_nodes <- sum(igraph::V(G_tmp)$type)
            nr_peptide_nodes <- sum(!igraph::V(G_tmp)$type)
            nr_edges <- igraph::gsize(G_tmp)

            nr_edges_per_pep_node <- igraph::degree(G_tmp)[
                !igraph::V(G_tmp)$type]
            nr_unique_peptides <- sum(nr_edges_per_pep_node == 1)
            nr_shared_peptides <- sum(nr_edges_per_pep_node > 1)


            protein_acc <- igraph::V(G_tmp)$name[igraph::V(G_tmp)$type]
            protein_acc <- strsplit(protein_acc, ";")
            nr_protein_accessions <- sum(vapply(protein_acc,
                    length, FUN.VALUE = numeric(1)))

            peptide_seq <- igraph::V(G_tmp)$name[!igraph::V(G_tmp)$type]
            peptide_seq <- strsplit(peptide_seq, ";")
            nr_peptide_sequences <- sum(vapply(peptide_seq, length,
                    FUN.VALUE = numeric(1)))

            unique_pept_nodes <- igraph::V(G_tmp)[
                (igraph::degree(G_tmp) == 1 & !igraph::V(G_tmp)$type)]

            ## Fall: keine uniquen Peptide im ganzen Graphen
            if (length(unique_pept_nodes) == 0) {
                nr_prot_node_only_unique_pep <- 0
                nr_prot_node_unique_and_shared_pep <- 0
                nr_prot_node_only_shared_pep <-  nr_protein_nodes
            } else {
                ## Fall: I-shaped graph
                if (length(unique_pept_nodes) == 1 && nr_protein_nodes == 1) {
                    nr_prot_node_only_unique_pep <- 1
                    nr_prot_node_unique_and_shared_pep <- 0
                    nr_prot_node_only_shared_pep <-  0
                } else {
                    ## neighborhood of the unique peptides
                    ## (these are proteins with a unique peptide)
                    NH_of_unique_peptides <- igraph::ego(G_tmp, order = 1,
                        mindist = 1, nodes = unique_pept_nodes)

                    nr_prot_node_only_unique_pep <- 0
                    ## = Anzahl uniquer Peptide??
                    nr_prot_node_unique_and_shared_pep <- length(
                        NH_of_unique_peptides)
                    ## length(NH_of_unique_peptides)
                    nr_prot_node_only_shared_pep <-  nr_protein_nodes -
                        nr_prot_node_unique_and_shared_pep
                }
            }

            ## TODO: add nr of unique and shared peptide sequences
            ## TODO: add info about graph type (isomorphism list!)

            D_tmp <- data.frame(graph_ID = i,
                                nr_protein_nodes = nr_protein_nodes,
                                nr_peptide_nodes = nr_peptide_nodes,
                                nr_unique_peptide_nodes = nr_unique_peptides,
                                nr_shared_peptide_nodes = nr_shared_peptides,
                                nr_edges = as.integer(nr_edges),
                                nr_protein_accessions = nr_protein_accessions,
                                nr_peptide_sequences =
                                    as.integer(nr_peptide_sequences),
                                ## nr_peptide_sequences_unique =
                                    ## nr_peptide_sequences_unique,
                                ## nr_peptide_sequences_shared =
                                    ## nr_peptide_sequences_shared,
                                nr_prot_node_only_unique_pep =
                                    nr_prot_node_only_unique_pep,
                                nr_prot_node_unique_and_shared_pep =
                                    nr_prot_node_unique_and_shared_pep,
                                nr_prot_node_only_shared_pep =
                                    nr_prot_node_only_shared_pep,
                                comparison = comparisons[j])

            Data <- rbind(Data, D_tmp)
            pbapply::setpb(pb, i)
        }
        ## progress bar command
        invisible(NULL)
    }

    if (prototype) {
        Data <- cbind(Data, counter = counter)
    }

    if (!is.null(file)) openxlsx::write.xlsx(Data, file, overwrite = TRUE)

    return(Data)
}
