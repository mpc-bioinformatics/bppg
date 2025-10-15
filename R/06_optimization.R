#' Functions in this file:
#' .errorEquation
#' .minimizeSquaredError
#' iterateOverCi
#' automatedAnalysisIteratedCi

#' Function to set up the error equations for the optimization problem
#'
#' @param Ri          \strong{numeric vector} \cr
#'                    Contains the (estimated) protein ratios.
#' @param Ci          \strong{numeric vector} \cr
#'                    Contains the protein weights (estimated, sum up to 1)
#' @param M           \strong{matrix} \cr
#'                    The biadjaceny matrix of the corresponding graphs.
#' @param rj          \strong{numeric vector} \cr
#'                    Contains the measured peptide ratios.
#' @param log_level   \strong{logical} \cr
#'                    If \code{TRUE}, the Ri are given on log2-level and need
#' to be back-transformed here
#' (this leads to a symmetric behavior during optimization)
#'
#' @return list containing the following elements:
#' \item{res_Mat}{matrix containing the estimated peptide
#'      ratios using the given Ri and Ci}
#' \item{res_equ}{vector of error terms for each peptide}
#' \item{res_squ_err}{sum of squared error terms}
#' \item{W}{internal weight matrix}
#'
#'
#'
#' @examples
#' Ri <- c(0.5, 1.3)
#' Ci <- c(0.3, 0.7)
#' M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
#' rj <- c(0.6, 1.2)
#' bppg:::.errorEquation(Ri, Ci, M, rj)
.errorEquation <- function(Ri,
    Ci,
    M,
    rj,
    log_level = FALSE) {
    m <- length(Ri) ## number of proteins
    n <- length(rj) ## number of peptides
    checkmate::assertNumeric(Ri)
    checkmate::assertNumeric(Ci, len = length(Ri))
    checkmate::assertNumeric(rj)
    checkmate::assertMatrix(M, ncols = length(Ri), nrows = length(rj),
                            mode = "integerish")
    stopifnot(all(M %in% c(0,1)))
    checkmate::assertFlag(log_level)

    ## backtransformation if necessary
    if (log_level) Ri <- 2^Ri
    ## multiply the delta values (biadjacency matrix M) with their weights Ci
    W <- sweep(M, MARGIN = 2, Ci, "*")
    ## sum of the weights per peptide
    W_sum <- rowSums(W)
    ## divide the weights by the sum of the weights per peptide
    W <- as.matrix(sweep(W, 1, W_sum, "/"))
    ## multiply Ri with the corresponding weight
    res_Mat <- sweep(W, MARGIN = 2, Ri, '*')

    ## error term per peptide (on log-scale)
    res_equ <- log(rj) - log(rowSums(res_Mat))

    ## sum of squared error terms
    res_squ_err <- sum(res_equ^2)

    return(list(res_Mat = res_Mat, res_equ = res_equ,
            res_squ_err = res_squ_err, W = W))
}



#' Calulate initial values for Ci (protein weights) for optimization
#'
#' @param fixed.Ci \strong{numeric vector} \cr
#'                    The fixed protein weights, variable weights set as NA.
#'                    Sum of fixed weights must not exceed 1.
#'                    If NULL, all Cis will be considered as variable.
#'                    This argument is needed to fix Ci on a grid point in the
#'                    iterated_Ci function.
#' @param m \strong{integer} \cr
#'                    number of proteins in the respective graph
#'
#' @returns Ci_start: vector with inital start values for Ci for the
#'                    optimization step
#'
#' @examples
.initializeCi <- function(fixed.Ci,
                          m) {
    #m <- length(fixed.Ci) # number of proteins
    is.Ci.fixed <- !is.null(fixed.Ci)
    which.Ci.fixed <- which(!is.na(fixed.Ci))

    ## Initialization of Ci:
    if (!is.Ci.fixed) {
        ## the algorithm starts with equal weights for each protein
        Ci_start <- rep(1 / m, m)
    } else {
        ## if at least one Ci is fixed, the algorithm distributes the remaining
        ## weight equally among the non-fixed proteins
        m2 <- m - length(which.Ci.fixed)
        ## sum of fixed Ci (as all Ci have to sum up tp 1)
        fixed.Ci.sum <- sum(fixed.Ci, na.rm = TRUE)
        Ci_start <- fixed.Ci
        ## starting values for the remaining Ci values
        Ci_start[is.na(Ci_start)] <- (1 - fixed.Ci.sum) / m2
    }
    return(Ci_start)
}




#' Calulate initial values for Ri (protein ratios) for optimization
#'
#' @param M \strong{matrix} \cr
#'                     Biadjacency matrix of the graph.
#' @param rj \strong{numeric vector} \cr
#'                   Vector of measured peptide ratios.
#' @param m \strong{integer(1)} \cr
#'                   number of proteins in the respective graph
#' @param log_level \strong{logical(1)} \cr
#'                  if TRUE, the Ri values are returned on log2-scale
#'
#' @returns Ri_start: vector with inital start values for Ri for the
#'                    optimization step
.initializeRi <- function(M, rj, m, log_level) {
    Ri_start <- rep(NA, m)
    for (j in 1:m) {
        tmp <- M * rj
        tmp[tmp == 0] <- NA    ## 0 -> peptide is not present in the protein
        uniquePep <- (rowSums(M) == 1) & (M[, j] == 1)
        if (any(uniquePep)) {
            Ri_start[j] <- .geomMean(tmp[uniquePep, j], na.rm = TRUE)
        } else {
            Ri_start[j] <- .geomMean(tmp[, j], na.rm = TRUE)
        }
    }
    if (log_level) Ri_start <- log2(Ri_start)
    return(Ri_start)
}




#' Calculate equality and inequality constraints for the optimization step
#'
#' @param fixed.Ci \strong{numeric vector} \cr
#'                    The fixed protein weights, variable weights set as NA.
#'                    Sum of fixed weights must not exceed 1.
#'                    If NULL, all Cis will be considered as variable.
#'                    This argument is needed to fix Ci on a grid point in the
#'                    iterated_Ci function.
#' @param log_level \strong{logical(1)} \cr
#'                  if TRUE, the Ri values are returned on log2-scale
#' @param m \strong{integer(1)} \cr
#'          number of proteins in the respective graph
#'
#' @returns list containing the following elements (see also
#'                  \code{\link[Rsolnp]{solnp}}):
#' \item{eqfun}{function for calculating inequality constraints}
#' \item{LB}{lower bound for eqfun}
#' \item{eqB}{vector of equality bounds for the variables}
.calcConstraints <- function(fixed.Ci, log_level, m) {
    if (!is.null(fixed.Ci)) {
        m2 <- m - sum(!is.na(fixed.Ci)) # number of free weights (not fixed)
        fixed.Ci.sum <- sum(fixed.Ci, na.rm = TRUE)
        eqfun <- function(x) sum(x[(m + 1):(m + m2)]) + fixed.Ci.sum - 1
        LB <- rep(0, m + m2)
        if (log_level) LB <- c(rep(-Inf, m), rep(0, m2))
        eqB <- 0
    } else {
        eqfun <- function(x) sum(x[(m + 1):(2 * m)]) - 1
        LB <- rep(0, 2 * m)
        if (log_level) LB <- c(rep(-Inf, m), rep(0, m))
        eqB <- 0
    }
    return(list(eqfun = eqfun, LB = LB, eqB = eqB))
}



#' Calculate objective function for the optimization step
#'
#' @param fixed.Ci    \strong{numeric vector} \cr
#'                    The fixed protein weights, variable weights set as NA.
#'                    Sum of fixed weights must not exceed 1.
#'                    If NULL, all Cis will be considered as variable.
#'                    This argument is needed to fix Ci on a grid point in the
#'                    iterated_Ci function.
#' @param log_level   \strong{logical(1)} \cr
#'                    if TRUE, the Ri values are returned on log2-scale
#' @param m           \strong{integer(1)} \cr
#'                    number of proteins in the respective graph
#' @param M           \strong{matrix} \cr
#'                    Biadjacency matrix of the graph.
#' @param rj          \strong{numeric vector} \cr
#'                    Contains the measured peptide ratios.
#'
#' @returns objective function that will be minimized
.calcObjectiveFunction <- function(fixed.Ci, log_level, m, M, rj) {
    if (!is.null(fixed.Ci)) {
        m2 <- m - sum(!is.na(fixed.Ci)) # number of free weights (not fixed)
        fun <- function(x) {
            Ri_tmp <- x[1:m]
            Ci_tmp <- fixed.Ci
            Ci_tmp[is.na(fixed.Ci)] <- x[(m + 1):(m + m2)]
            res <- .errorEquation(Ri = Ri_tmp, Ci = Ci_tmp, M = M, rj = rj,
                           log_level = log_level)$res_squ_err
            return(res)
        }
    } else {
        fun <- function(x) {
            Ri_tmp <- x[1:m]
            Ci_tmp <- x[(m + 1):(2 * m)]
            res <- .errorEquation(Ri = Ri_tmp, Ci = Ci_tmp, M = M, rj = rj,
                           log_level = log_level)$res_squ_err
            return(res)
        }
    }
    return(fun)
}



#' Function to set up the optimization problem and minimize the sum of squared
#' error terms
#'
#' @param S           \strong{list} \cr
#'                    A list of biadjacency matrix of the bipartite
#'                    peptide-protein graph (X) and
#'                    measured peptide ratios (fc).
#' @param fixed.Ci    \strong{numeric vector} \cr
#'                    The fixed protein weights, variable weights set as NA.
#'                    Sum of fixed weights must not exceed 1.
#'                    If NULL, all Cis will be considered as variable.
#'                    This argument is needed to fix Ci on a grid point in the
#'                    iterated_Ci function.
#' @param verbose     \strong{logical} \cr
#'                    If \code{TRUE}, additional information on each iteration
#'                    of the optimization is printed (see also rsolnp function
#'                    in package Rsolnp).
#' @param reciprocal  \strong{logical} \cr
#'                    If \code{TRUE}, the reciprocal of the peptide ratios is
#'                    used for the optimization.
#' @param log_level   \strong{logical} \cr
#'                    If \code{TRUE}, the Ri are log2-transformed before
#'                    optimization, allowing a symmetric consideration of
#'                    Ri < 0 and > 0.
#' @param control     \strong{list} \cr
#'                    The control parameters for solnp.
#'
#' @return list containing the following elements:
#' \item{Ri}{estimated protein ratios}
#' \item{Ci}{estimated protein weights}
#' \item{RES}{final result of .errorEquation(), which also contains the final,
#'  minimal error term}
#' \item{Tracking}{Tracking of Ri, Ci and error term for the
#'  different iterations}
#' \item{outer.iter}{Number of outer iterations needed for the optimization
#'  algorithm to converge or stop}
#' \item{convergence}{Indicates whether the solver has converged (0) or
#' not (1 or 2).}
#'
#' @examples
#' M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
#' rj <- c(0.6, 1.2)
#' S <- list(X = M, fc = rj)
#' bppg:::.minimizeSquaredError(S)
#'
## TODO SIMPLFY -> put options into differnt functions?
## TODO: parameter checks for this private function useful?
.minimizeSquaredError <- function(S,
    fixed.Ci = NULL,
    verbose = FALSE,
    log_level = TRUE,
    control = list(),
    ...) {
    M <- S$X  ## biadjacency matrix
    m <- ncol(S$X) ## number of proteins
    n <- nrow(S$X) ## number of peptides
    rj <- S$fc     ## given peptide ratios
    checkmate::assertList(S, types = c("matrix", "numeric"))
    checkmate::assertMatrix(M, mode = "integerish")
    stopifnot(all(S$X %in% c(0,1)))
    checkmate::assertNumeric(rj)
    checkmate::assertNumeric(fixed.Ci, len = length(rj), lower = 0, upper = 1,
                             null.ok = TRUE)
    stopifnot(sum(fixed.Ci, na.rm = TRUE) <= 1)
    checkmate::assertFlag(verbose)
    checkmate::assertFlag(log_level)
    checkmate::assertList(control)
    if (!verbose) control <- c(control, trace = 0)
    is.Ci.fixed <- !is.null(fixed.Ci)
    if (is.Ci.fixed) which.Ci.fixed <- which(!is.na(fixed.Ci))

    Ci_start <- .initializeCi(fixed.Ci, m)
    Ri_start <- .initializeRi(M, rj, m, log_level)
    if (is.Ci.fixed) {
        pars <- c(Ri_start, Ci_start[-which.Ci.fixed])
    } else {
        pars <- c(Ri_start, Ci_start)
    }
    ## initial error term
    RES <- .errorEquation(Ri = Ri_start, Ci = Ci_start, M = M, rj = rj,
        log_level = log_level)
    ### TODO: do we need tracking?
    track_colnames <- c("iter", "squ_err", paste0("R", 1:m), paste0("C", 1:m))
    Tracking <- matrix(c(0, RES$res_squ_err, Ri_start, Ci_start), nrow = 1)
    Tracking <- as.data.frame(Tracking)
    colnames(Tracking) <- track_colnames

    fun <- .calcObjectiveFunction(fixed.Ci, log_level, m, M, rj)
    constr <- .calcConstraints(fixed.Ci, log_level, m)

    res <- Rsolnp::solnp(pars = pars, fun = fun, LB = constr$LB,
                         eqfun = constr$eqfun, eqB = constr$eqB,
                         control = control)
    # extract optimal Ri and Ci values from optimization result
    Ri <- res$pars[1:m]
    if (is.Ci.fixed) {
        m2 <- m - sum(!is.na(fixed.Ci)) # number of free weights (not fixed)
        Ci_tmp <- res$pars[(m + 1):(m + m2)]
        Ci <- fixed.Ci
        Ci[is.na(Ci)] <- Ci_tmp
    } else {
        Ci <- res$pars[(m + 1):(2 * m)]
    }
    ## update RES
    RES <- .errorEquation(Ri = Ri, Ci = Ci, M = M, rj = rj,
        log_level = log_level)
    Tracking <- rbind(Tracking, c(1, RES$res_squ_err, Ri, Ci))
    if (log_level) Ri <- 2^Ri
    result <- list(Ri = Ri, Ci = Ci, RES = RES, Tracking = Tracking,
        outer.iter = res$outer.iter, convergence = res$convergence)
    return(result)
}




#' Iterate over possible Ci values
#'
#' @param S                        \strong{list} \cr
#'                                 A list of biadjacency matrix of the bipartite
#'                                 peptide-protein graph (named "X")
#'                                 and measured peptide ratios (named "fc").
#' @param grid.size                \strong{integer} \cr
#'                                 The number of grid points for the Cis.
#' @param omit_grid_borders        \strong{logical} \cr
#'                                 If \code{TRUE}, omit exact value of 1 and 0
#'                                 from the grid (recommended, as they may cause
#'                                 numerical issues).
#' @param grid.start               \strong{integer} \cr
#'                                 The start of the grid (default is 0).
#' @param grid.stop                \strong{integer} \cr
#'                                 The end of the grid (default is 1).
#' @param verbose                  \strong{logical} \cr
#'                                 If \code{TRUE}, print additional information
#'                                 (see solnp function).
#' @param verbose_opt              \strong{logical} \cr
#'                                 The \code{verbose} argument of the
#'                                 [.minimizeSquaredError()] function.
#' @param control                  \strong{list} \cr
#'                                 The \code{control} object to be passed to the
#'                                 [.minimizeSquaredError()] function.
#' @param extend_grid_at_borders   \strong{logical} \cr
#'                                 If \code{TRUE}, the grid will be extend close
#'                                 to the borders (0 and 1).
#'                                 While exact values of 0 and 1 may cause
#'                                 numerical problems, values close to those may
#'                                 be valuable to get a better estimate of the
#'                                 protein ratios.
#' @param log_level                \strong{logical} \cr
#'                                 The \code{log_level} argument will passed to
#'                                 the [.minimizeSquaredError()] function.
#'
#' @return
#' A dataframe containing the optimal Ci and Ri values together with the reached
#' minimal error term.
#'
#' @export
#'
#' @seealso [.minimizeSquaredError()],
#'          [bppg::automatedAnalysisIteratedCi()]
#'
#' @details
#' With .minimizeSquaredError() each protein node receives one estimate for the
#' protein ratio. However, in some cases, there are multiple possible values for
#' the protein ratios that lead to the same, minimal error term.
#' To get a better coverage of the optimal solutions, this function iterates
#' over a grid of possible weights Ci for each protein node.
#' For each protein node, the Ci value if fixed on a point on the grid, while
#' the other Ci values and all Ri values are optimized using
#' .minimizeSquaredError().
#' The resulting table can be used to assess a range of possible solutions for
#' the protein ratios.
#'
#'
#' @examples
#' M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
#' rj <- c(0.6, 1.2)
#' S <- list(X = M, fc = rj)
#' bppg:::.minimizeSquaredError(S)
#' ## example not complete? TODO
iterateOverCi <- function(S,
    grid.size = 1000,
    omit_grid_borders = TRUE,
    grid.start = 0,
    grid.stop = 1,
    verbose = FALSE,
    verbose_opt = FALSE,
    control = list(),
    extend_grid_at_borders = FALSE,
    log_level = TRUE) {

    n <- ncol(S$X) ## number of protein groups
    grid <- seq(grid.start, grid.stop, length.out = grid.size + 1)

    if (extend_grid_at_borders) {
        grid_min <- grid[2] ## 2. Element, da erstes 0
        grid_max <- grid[length(grid) - 1] ## 1. Element, da letztes 1

        grid_extend_min <- seq(grid.start, grid_min, length.out = 11)
        grid_extend_max <- seq(grid_max, grid.stop, length.out = 11)

        grid <- sort(unique(c(grid, grid_extend_min, grid_extend_max)))
    }

    if (omit_grid_borders) grid <- grid[-c(1, length(grid))]

    err_tmp <- rep(NA, length(grid))
    Ris_tmp <- matrix(nrow = length(grid), ncol = n)
    colnames(Ris_tmp) <- paste0("R", 1:n)
    Cis_tmp <- matrix(nrow = length(grid), ncol = n)
    colnames(Cis_tmp) <- paste0("C", 1:n)
    result <- NULL

    pb <- pbapply::startpb(0, n * length(grid))
    ## TODO VAPPLY?
    for (j in 1:n) {
        for (i in seq_along(grid)) {

            pbapply::setpb(pb, (j - 1) * length(grid) + i)

            if (verbose) print(paste0("j = ", j, " i = ", i))

            Ci_tmp <- rep(NA, n)
            Ci_tmp[j] <- grid[i]

            RES <- try({
                .minimizeSquaredError(S,
                    fixed.Ci = Ci_tmp,
                    verbose = verbose_opt,
                    #reciprocal = FALSE,
                    control = control,
                    log_level = log_level)
            })
            if ("try-error" %in% class(RES)) {
                if (grepl("reached elapsed time limit", RES)) {
                stop(paste0("protein ", j, " grid point ", i))
                } else {
                next
                }
            }
            err_tmp[i] <- RES$RES$res_squ_err
            Ris_tmp[i, ] <- RES$Ri
            Cis_tmp[i, ] <- RES$Ci
        }
        result_tmp <- data.frame(protein = rep(j, length(grid)),
            grid = grid, Ris_tmp, Cis_tmp, error = err_tmp)
        result <- rbind(result, result_tmp)
    }
    invisible(NULL)
    pbapply::closepb(pb)

    return(result)
}


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
#' @param res                           \strong{list} \cr
#'                                      The list resulting from the
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

automatedAnalysisIteratedCi <- function(S,
    res,
    use_results_from_other_proteins = FALSE,
    verbose = FALSE,
    job = NULL,
    S_is_graph = FALSE) {

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

        if (verbose) print(paste0("Protein ", i))

        X_tmp <- res[res$protein == i, ]
        X_tmp <- X_tmp[!is.na(X_tmp$error), ]
        C_tmp <- X_tmp[, paste0("C", i)]

        ## use results from other proteins but only when Ci is not too extreme
        X_tmp3 <- res[res$protein != i,]
        X_tmp3 <- X_tmp3[!is.na(X_tmp3$error),]
        C_tmp3 <- X_tmp3[, paste0("C", i)]
        X_tmp3 <- X_tmp3[X_tmp3$protein == i |
            (C_tmp3 >= 0.01 & C_tmp3 < 0.99),]
        R_3 <- X_tmp3[, paste0("R", i)]
        C_3 <- X_tmp3[, paste0("C", i)]

        error <- X_tmp$error
        R <- X_tmp[, paste0("R", i)]
        C <- X_tmp[, paste0("C", i)]

        D_tmp$min_error <- min(error_optimal, min(error), na.rm = TRUE)
        D_tmp$error_optimal <- error_optimal

        ## 1st check: is error constant?
        if (abs(diff(range(X_tmp$error))) <= 1e-10) {
            ## constant error, i.e. single solution or interval with lower and
            ## upper border

            D_tmp$error_constant <- "yes"

            if (verbose) print(paste0("Error is nearly constant (abs. diff. = ",
                    abs(diff(range(error))), ")."))
            if (verbose) print(paste0("Mean error is ",
                    mean(error, na.rm = TRUE), "."))
            if (verbose) print(paste0("Minimal error is ",
                    min(error, na.rm = TRUE) , "."))
            ind_min <- which.min(error)

            ## 2nd check: Is Ri constant too?
            if (abs(diff(range(R))) > 1e-4) {
                if (verbose) print(paste0("Range for R", i, ": ",
                        BBmisc::collapse(range(R), sep = " - "), "."))
                D_tmp$Ri <- NA
                D_tmp$Ri_min <- min(R)
                D_tmp$Ri_max <- max(R)
                D_tmp$case <- 1
            } else {
                if (verbose) print(paste0("Constant Solution for R", i, ": ",
                        .geomMean(R)))
                D_tmp$Ri <- .geomMean(R)
                D_tmp$Ri_min <- NA
                D_tmp$Ri_max <- NA
                D_tmp$case <- 2
            }
            if (verbose) print(paste0("Range for C", i, ": ",
                    BBmisc::collapse(range(C), sep = " - "), "."))
            D_tmp$Ci_min <- min(C)
            D_tmp$Ci_max <- max(C)

        } else {
            error <- X_tmp$error
            R <- X_tmp[, paste0("R", i)]
            C <- X_tmp[, paste0("C", i)]

            if (verbose & !is.na(error_optimal) & any(error < error_optimal)) {
                warning(paste0("Lower than optimal error detected.
                    Difference = ", abs(min(error) - error_optimal)))
            }

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

                if (verbose) print(paste0("Single optimal solution with error ",
                        min(error, na.rm = TRUE), "."))

                if (verbose) print(paste0("Single optimal solution is R",
                        i, "= ", Ri_optimal[i], " and C",
                        i, " = ", Ci_optimal[i], "."))
            } else {  ## multiple data points with nearly constant error
                D_tmp$ error_constant <- "partially"

                if (verbose) print(paste0("Multiple points with
                        optimal error."))

                ## Ri is constant but Ci is not
                if (all(R[ind_min_tol] == 0) |
                        abs(diff(range(log2(R[ind_min_tol])))) < 1e-4) {
                    D_tmp$Ri <- .geomMean(R[ind_min_tol])
                    D_tmp$Ri_min <- NA
                    D_tmp$Ri_max <- NA
                    D_tmp$case <- 4
                    if (verbose) print(paste0("Constant Solution for R",
                            i, ": ", .geomMean(R[ind_min_tol])))
                } else {  ## Ri is not constant
                    D_tmp$Ri <- NA
                    D_tmp$Ri_min <- min(R[ind_min_tol])
                    D_tmp$Ri_max <- max(R[ind_min_tol])
                    D_tmp$case <- 5
                    if (verbose) print(paste0("Range for R", i, ": ",
                            BBmisc::collapse(range(R[ind_min_tol]),
                                sep = " - "), "."))
                }
                if (verbose) print(paste0("Range for C", i, ": ",
                        BBmisc::collapse(range(C[ind_min_tol]),
                            sep = " - "), "."))
                D_tmp$Ci <- NA
                D_tmp$Ci_min <- min(C[ind_min_tol])
                D_tmp$Ci_max <- max(C[ind_min_tol])
            }
        }

        if (use_results_from_other_proteins) {
            ## see if solution can be enhanced by data form the other proteins
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

