#' Functions in this file:
#' .addAveragePepRatio
#' .addUniquenessAttributes
#' .geomMean

#' Adds average peptide ratios as a attribute to the graphs, if a list of
#' peptide ratios is already present.
#'
#' @param G      \strong{igraph graph object} \cr
#'               A peptide-protein graph.
#' @param type   \strong{character} \cr
#'               !NOT USED AT THE MOMENT!
#'
#' @return A graph with added peptide ratio attributes.
#'
#'
#' @seealso [generateGraphsFromFASTA()], [generateQuantGraphs()],
#'          [.addUniquenessAttributes()]
#'
#' @examples
#' 
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' edgelist <- bppg::digestFASTA(fasta)
#'
#' file <- system.file("extdata", "peptides.txt", package = "bppg")
#' D <- bppg::readMqPeptideTable(path = file, LFQ = TRUE, remove_contaminants = FALSE)
#' group <- factor(rep(1:9, each = 3))
#' dAgg <- bppg::aggregateReplicates(D, group = group)
#' exp_peptide_ratios <- bppg::calculatePeptideRatios(dAgg)
#' graph <- bppg::generateQuantGraphs(exp_peptide_ratios, edgelist)
#'
#' res <- bppg:::.addAveragePepRatio(graph[[1]])

.addAveragePepRatio <- function(G, type = "geom_mean") {

    pep_ratio <- igraph::V(G)$pep_ratio
    pep_ratio_split <- strsplit(pep_ratio, ";")

    pep_ratio_aggr <- sapply(pep_ratio_split, function(x) {
        .geomMean(as.numeric(x))})

    nr_sequences <- sapply(pep_ratio_split, length)

    G <- igraph::set_vertex_attr(G, "pep_ratio_aggr", value = pep_ratio_aggr)
    G <- igraph::set_vertex_attr(G, "nr_sequences", value = nr_sequences)
    return(G)
}

#' Adds vertex attributes with uniqueness of peptides and number of unique
#' peptides for proteins.
#'
#' @param G \strong{igraph graph object} \cr
#'          A peptide-protein graph.
#'
#' @return A graph with 2 additional vertex attributes, uniqueness and
#'         nr_unique_peptides
#'
#'
#' @seealso [generateGraphsFromFASTA()], [generateQuantGraphs()],
#'          [.addAveragePepRatio()]
#'
#' @examples
#' 
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' edgelist <- digestFASTA(fasta)
#' graph <- bppg::generateGraphsFromEdgelist(edgelist)
#' 
#' res <- bppg:::.addUniquenessAttributes(graph[[1]])

.addUniquenessAttributes <- function(G) {
    ## FALSE = peptide, TRUE = protein
    igraph::V(G)$type

    uniqueness <- igraph::degree(G, igraph::V(G)) == 1
    ## attribute only for peptides
    uniqueness[igraph::V(G)$type] <- NA

    G <- igraph::set_vertex_attr(G, "uniqueness", value = uniqueness)

    unique_peptide_nodes <- igraph::V(G)[igraph::V(G)$uniqueness &
            !is.na(igraph::V(G)$uniqueness)]
    shared_peptide_nodes <- igraph::V(G)[!igraph::V(G)$uniqueness &
            !is.na(igraph::V(G)$uniqueness)]

    neighborhood <- igraph::ego(G, order = 1, mindist = 1, nodes = igraph::V(G))

    ## TODO VAPPLY
    nr_unique_peptides <- sapply(neighborhood, function(x) {
        sum(x %in% unique_peptide_nodes)
    })
    nr_unique_peptides[!igraph::V(G)$type] <- NA ## attribute only for proteins
    G <- igraph::set_vertex_attr(G, "nr_unique_peptides",
        value = nr_unique_peptides)

    nr_shared_peptides <- sapply(neighborhood, function(x) {
        sum(x%in% shared_peptide_nodes)
    })
    nr_shared_peptides[!igraph::V(G)$type] <- NA ## attribute only for proteins
    G <- igraph::set_vertex_attr(G, "nr_shared_peptides",
        value = nr_shared_peptides)
}

#' Calculate the geometric mean.
#'
#' @param x         \strong{numeric vector} \cr
#'                  Input data.
#' @param useprod   \strong{logical} \cr
#'                  If \code{TRUE}, prod(x)^(1/n) will be calculated, otherwise
#'                  exp(mean(log(x))).
#'
#' @return The geometric mean of the provided data points.
#'
#'
#' @examples
#' data <- c(1,6,3.5)
#' result <- bppg:::.geomMean(data, useprod = FALSE)

.geomMean <- function(x, useprod = FALSE, na.rm = FALSE) {
    n <- length(x)

    if (useprod) {
        return(prod(x, na.rm = na.rm)^(1 / n))
    } else {
        return(exp(mean(log(x), na.rm = na.rm)))
    }
}
