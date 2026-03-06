# Functions in this file:
# .addAveragePepRatio
# .addUniquenessAttributes
# .geomMean


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
#' @seealso [generateGraphsFromFASTA()], [generateQuantGraphs()]
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
    nr_unique_peptides <- vapply(neighborhood, function(x) {
        sum(x %in% unique_peptide_nodes)
    }, FUN.VALUE = numeric(1))
    nr_unique_peptides[!igraph::V(G)$type] <- NA ## attribute only for proteins
    G <- igraph::set_vertex_attr(G, "nr_unique_peptides",
        value = nr_unique_peptides)

    nr_shared_peptides <- vapply(neighborhood, function(x) {
        sum(x%in% shared_peptide_nodes)
    }, FUN.VALUE = numeric(1))
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
#' @param na.rm     \strong{logical} \cr
#'                  If \code{TRUE}, missing values are removed before the 
#'                  calculation. \code{FALSE} is default.
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
