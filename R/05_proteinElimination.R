#' Functions in this file:
#' proteinElimination



#' Tries to recursively remove protein nodes from graph while keeping the error
#' term in the optimization step low
#'
#' @param G                    \strong{igraph graph object} \cr
#'                             The graph to eliminate the proteins in.
#' @param threshold            \strong{numeric} \cr
#'                             The threshold for increase of error term.
#'                             The default 1.05 refers to 5% increase.
#' @param control              \strong{list} \cr
#'                              A list of control parameters for the
#'                              optimization step.
#'                              See \code{\link[Rsolnp]{solnp}} for details.
#' @param min_error_ref        \strong{numeric} \cr
#'                             The minimal error term on the whole graph using
#'                             all available protein nodes.
#' @param resDF                \strong{data.frame} \cr
#'                             Dataframe for collecting results from each iteration.
#'                             It contains information on the tested combinations
#'                             of protein nodes, the error terms and the currently
#'                             best combination.
#' @param protnodes_list       \strong{igraph node list} \cr
#'                             A current list of protein nodes.
#' @param res_best             \strong{list} \cr
#'                             A list containg the currently best solution including
#'                             the protein node combination, the error term and
#'                             the corresponding graph object.
#'
#' @return list containing the following elements:
#' \item{min_error_ref}{reference error term of the whole graph}
#' \item{protnodes_list}{list of all available protein nodes in the beginning}
#' \item{resDF}{dataframe with results of all iterations}
#' \item{res_best}{List of the overall best solution.}
#'
#'
#' @details
#' This function works in a recursive way. In the first iteration, the error
#' term on the whole graph is assessed. Then, one of the protein nodes is
#' deleted and the function is recursively applied. The results of all
#' iterations are collected in the resDF data.frame.
#'
#'  For starting the first iteration, only G, threshold and if necessary control
#'  have to be defined, everything else will be calculated during the first
#'  iteration for all future iterations.
#'
#'
#' @export
#'
#' @seealso [bppg::.minimizeSquaredError()]
#'
#' @examples ## TODO
#' file <- system.file("extdata", "quantGraphsForTesting.rds", package = "bppg")
#' graphs <- readRDS(file)
#' G <- graphs$sample1_sample2[[2]]
#' proteinElimination(G)
proteinElimination <- function(G,
    threshold = 1.05,
    control = list(), #list(trace = 0, delta = 1e-9),  # TODO! list()
    min_error_ref = NULL, # error from iteration 1
    resDF = NULL,
    protnodes_list = NULL, # based on the very first original graph in iter 1
    res_best = NULL) {
    checkmate::assertClass(G, classes = c("igraph"))
    checkmate::checkTRUE(igraph::is_bipartite(G))
    checkmate::assertNumeric(threshold)
    checkmate::assertList(control)
    checkmate::assertNumeric(min_error_ref, null.ok = TRUE)
    checkmate::assertDataFrame(resDF, null.ok = TRUE)
    checkmate::assertClass(protnodes_list, classes = "igraph.vs", null.ok = TRUE)
    checkmate::assertList(res_best, null.ok = TRUE)

    if (is.null(resDF)) { # first iteration
        G <- .addUniquenessAttributes(G)
        protnodes <- igraph::V(G)[igraph::V(G)$type]
        nr_unique_peptides <- igraph::V(G)$nr_unique_peptides[igraph::V(G)$type]
        min_error_ref <- .minimizeSquaredError(G, fixed.Ci = NULL,
            verbose = FALSE, control = control)$RES$res_squ_err
        resDF <- data.frame(comb = paste(protnodes, collapse = ","),
            n_proteins = length(protnodes), error = min_error_ref,
            current_best = TRUE)
        res_best <- list(G = list(G), comb = paste(protnodes, collapse = ","),
                         n_comb = length(protnodes), error = min_error_ref)
        protnodes_list <- protnodes
    }

    G <- .addUniquenessAttributes(G)
    nr_unique_peptides <- igraph::V(G)$nr_unique_peptides[igraph::V(G)$type]
    protnodes <- igraph::V(G)[igraph::V(G)$type]

    ## TODO: for loop into apply? -> too complicated for this recursive function???
    for (i in seq_along(protnodes)) {
        protnodes_tmp <- protnodes_list[-i]
        combination <- paste(protnodes_tmp, collapse = ",")
        if (combination %in% resDF$comb) next # skip already seen combinations
        res_tmp <- list(comb = combination,
                        n_proteins = length(protnodes_tmp),
                        error = NA, current_best = FALSE)
        if (nr_unique_peptides[i] > 0) { ## skip if deleted protein has unique peptides
            resDF <- rbind(resDF, res_tmp)
            next
        }
        G_tmp <- igraph::delete_vertices(G, protnodes[i])
        G_CC <- igraph::decompose(G_tmp)  #
        min_error_tmp <- 0

        errors <- vapply(G_CC, function(x) {
            .minimizeSquaredError(x, fixed.Ci = NULL, verbose = FALSE,
                control = control)$RES$res_squ_err}, FUN.VALUE = numeric(1))
        min_error_tmp <- sum(errors)
        res_tmp$error <- min_error_tmp
        resDF <- rbind(resDF, res_tmp)

        ## skip if error is NA or larger than the reference * threshold
        if (is.na(min_error_tmp) | min_error_tmp > min_error_ref * threshold) {
            next
        }
        if (res_tmp$n_proteins <= res_best$n_comb) {
            resDF$current_best[nrow(resDF)] <- TRUE
            res_best <- list(G = G_CC, comb = combination,
                             n_comb = res_tmp$n_proteins, error = res_tmp$error)
        }

        RES <- proteinElimination(G = G_tmp, threshold = threshold, control = control,
                                  min_error_ref = min_error_ref, resDF = resDF,
                                  protnodes_list = protnodes_tmp,
                                  res_best = res_best)
        resDF <- RES$resDF
        res_best <- RES$res_best
    }
    return(list(min_error_ref = min_error_ref, protnodes_list = protnodes,
        resDF = resDF, res_best = res_best))
}



