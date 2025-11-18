#' Extract protein ratio solutions from the result of iterateOverCi()
#'
#' @param S                             \strong{igraph graph object OR list} \cr
#'                                      An igraph graph of the bipartite
#'                                      peptide-protein graph with peptide
#'                                      ratios
#'                                      OR
#'                                      a list of the biadjacency matrix of the
#'                                      bipartite peptide-protein graph
#'                                      (named "X") and the measured peptide
#'                                      ratios (named "fc"). \cr
#'                                      Set \code{S_is_graph} depending on the
#'                                      input type.
#' @param res                           \strong{data.frame} \cr
#'                                      The data.frame resulting from the
#'                                      [bppg::iterateOverCi()] function.
#' @param use_results_from_other_proteins   \strong{logical} \cr
#'                                          If \code{TRUE}, the results from
#'                                          other proteins within the same graph
#'                                          will be used to calculate the
#'                                          optimal solution for each protein
#'                                          node.
#' @param verbose                       \strong{logical} \cr
#'                                      If \code{TRUE}, additional information
#'                                      will be printed.
#' @param job                           \strong{BatchExperiment job object} \cr
#'                                      Is used to print the job id and
#'                                      parameters in the output
#' @param S_is_graph                    \strong{logical} \cr
#'                                      If \code{TRUE}, S is an igraph object
#'                                      and if \code{FALSE} S is a list with the
#'                                      biadjacency matrix and fold changes.
#'
#' @return A data frame
#' @export
#'
#' @seealso [bppg::iterateOverCi()], [.minimizeSquaredError()]
#'
#' @examples ## TODO

automatedAnalysisIteratedCi_old <- function(S,
                                            res,
                                            use_results_from_other_proteins = FALSE,
                                            verbose = FALSE,
                                            job = NULL,
                                            S_is_graph = FALSE,
                                            error_tol = 1e-10,
                                            ratio_tol = 1e-4) {

    if (S_is_graph & !is.null(S)) {
        X <- igraph::as_biadjacency_matrix(S)
        fc <- stats::na.omit(igraph::V(S)$pep_ratio)
        S <- list(X = X, fc = fc)
    }

    if (!is.null(S)) {
        nr_proteins <- ncol(S$X)
        if (nr_proteins == 1) res <- NULL
    } else {
        nr_proteins <- 1
        res <- NULL
    }

    error_optimal <- NA
    Ri_optimal <- rep(NA, nr_proteins)
    Ci_optimal <- rep(NA, nr_proteins)

    if (nr_proteins == 1 & is.null(res)) {
        ## graph contains only one protein node and was skipped, so optimal
        ## solution has to be calculated here
        solution_optimal <- .minimizeSquaredError(S,
                                                  fixed.Ci = NULL,
                                                  verbose = FALSE,
                                                  # reciprocal = FALSE,
                                                  log_level = TRUE,
                                                  control = list(trace = 0))

        error_optimal <- solution_optimal$RES$res_squ_err
        Ri_optimal <- solution_optimal$Ri
        Ci_optimal <- solution_optimal$Ci

        D_tmp <- list(Accession = colnames(S$X),
                      comparison = NA,
                      graphID = NA,
                      proteinNr = 1,
                      error_optimal = error_optimal,
                      Ri_optimal = Ri_optimal[1],
                      Ci_optimal = Ci_optimal[1],
                      min_error = NA,
                      error_constant = NA,
                      Ri = Ri_optimal,
                      Ri_min = NA,
                      Ri_max = NA,
                      Ci = Ci_optimal,
                      Ci_min = NA,
                      Ci_max = NA,
                      case = 6)

        if(!is.null(job)) {
            D_tmp$graphID <- job$pars$prob.pars$k
            D_tmp$comparison <- job$prob.name
            D_tmp$job.id <- job$job.id
        }

        RES <- BBmisc::convertListOfRowsToDataFrame(list(D_tmp))
        return(RES)
    }

    D <- list()

    ############################################################################
    ###

    ## remove results where error is NA
    res <- res[!is.na(res$error),]


    for (i in 1:nr_proteins) {
        D_tmp <- list(Accession = colnames(S$X)[i],
                      comparison = NA,
                      graphID = NA,
                      proteinNr = i,
                      error_optimal = NA,
                      Ri_optimal = Ri_optimal[i],
                      Ci_optimal = Ci_optimal[i],
                      min_error = NA,
                      error_constant = NA,
                      Ri = NA,
                      Ri_min = NA,
                      Ri_max = NA,
                      Ci = NA,
                      Ci_min = NA,
                      Ci_max = NA,
                      case = 6)

        if (!is.null(job)) {
            D_tmp$graphID <- job$pars$prob.pars$k
            D_tmp$comparison <- job$prob.name
            D_tmp$job.id <- job$job.id
        }

        X_tmp <- res[res$protein == i, ] # results from current protein

        error <- X_tmp$error
        R <- X_tmp[, paste0("R", i)]
        C <- X_tmp[, paste0("C", i)]

        D_tmp$min_error <- min(error)
        # D_tmp$error_optimal <- error_optimal

        ## 1st check: is error constant?
        if (abs(diff(range(X_tmp$error))) <= error_tol) {
            ## constant error, i.e. single solution or interval with lower and
            ## upper border

            D_tmp$error_constant <- "yes"

            ind_min <- which.min(error)

            ## 2nd check: Is Ri constant too?
            if (abs(diff(range(R))) > ratio_tol) {

                D_tmp$Ri <- NA
                D_tmp$Ri_min <- min(R)
                D_tmp$Ri_max <- max(R)
                D_tmp$case <- 1
            } else {
                D_tmp$Ri <- .geomMean(R)
                D_tmp$Ri_min <- NA
                D_tmp$Ri_max <- NA
                D_tmp$case <- 2
            }

            D_tmp$Ci_min <- min(C)
            D_tmp$Ci_max <- max(C)

        } else {
            error <- X_tmp$error
            R <- X_tmp[, paste0("R", i)]
            C <- X_tmp[, paste0("C", i)]

            ## area in which the error is almost constant
            ind_min <- which.min(error)
            ind_min_tol <- which(abs(min(error) - error) <= 1e-10)

            if (length(ind_min_tol) <= 1) {

                D_tmp$error_constant <- "no"
                D_tmp$Ri <- R[ind_min_tol]
                D_tmp$Ri_min <- NA
                D_tmp$Ri_max <- NA
                D_tmp$Ci <- C[ind_min_tol]
                D_tmp$Ci_min <- NA
                D_tmp$Ci_max <- NA
                D_tmp$case <- 3


            } else {  ## multiple data points with nearly constant error
                D_tmp$error_constant <- "partially"


                ## Ri is constant but Ci is not
                if (all(R[ind_min_tol] == 0) |
                    abs(diff(range(log2(R[ind_min_tol])))) < ratio_tol) {
                    D_tmp$Ri <- .geomMean(R[ind_min_tol])
                    D_tmp$Ri_min <- NA
                    D_tmp$Ri_max <- NA
                    D_tmp$case <- 4

                } else {  ## Ri is not constant
                    D_tmp$Ri <- NA
                    D_tmp$Ri_min <- min(R[ind_min_tol])
                    D_tmp$Ri_max <- max(R[ind_min_tol])
                    D_tmp$case <- 5

                }

                D_tmp$Ci <- NA
                D_tmp$Ci_min <- min(C[ind_min_tol])
                D_tmp$Ci_max <- max(C[ind_min_tol])
            }
        }

        if (use_results_from_other_proteins) {
            ## use results from other proteins but only when Ci is not too extreme
            X_tmp3 <- res[res$protein != i,]
            # X_tmp3 <- X_tmp3[!is.na(X_tmp3$error),]
            C_tmp3 <- X_tmp3[, paste0("C", i)]
            X_tmp3 <- X_tmp3[X_tmp3$protein == i |
                                 (C_tmp3 >= 0.01 & C_tmp3 < 0.99),]
            R_3 <- X_tmp3[, paste0("R", i)]
            C_3 <- X_tmp3[, paste0("C", i)]

            ## see if solution can be enhanced by data from the other proteins
            ind <- which(abs(min(error) - X_tmp3$error) <= 1e-10)
            if (length(ind) > 0) {
                if (!is.na(D_tmp$Ri)) {
                    if (abs(diff(range(
                        log2(c(D_tmp$Ri, R_3[ind]))))) >= 1e-04) {
                        ## if there was 1 solution before there are multiple now
                        D_tmp$Ri_min <- min(c(D_tmp$Ri, R_3[ind]))
                        D_tmp$Ri_max <- max(c(D_tmp$Ri, R_3[ind]))
                        D_tmp$Ri <- NA
                    }
                } else {
                    ## if the range of solutions can be enhanced
                    D_tmp$Ri_min <- min(c(D_tmp$Ri_min, R_3[ind]))
                    D_tmp$Ri_max <- max(c(D_tmp$Ri_max, R_3[ind]))
                }
            }
        }
        D[[i]] <- D_tmp

    }

    RES <- BBmisc::convertListOfRowsToDataFrame(D)
    return(RES)
}
