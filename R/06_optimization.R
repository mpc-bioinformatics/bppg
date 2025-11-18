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
.errorEquation <- function(RiLog,
    Ci,
    M,
    rjLog) { #},
    #log_level = FALSE) {
    m <- length(RiLog) ## number of proteins
    n <- length(rjLog) ## number of peptides
    checkmate::assertNumeric(RiLog)
    checkmate::assertNumeric(Ci, len = length(RiLog))
    checkmate::assertNumeric(rjLog)
    checkmate::assertMatrix(M, ncols = length(RiLog), nrows = length(rjLog),
                            mode = "integerish")
    stopifnot(all(M %in% c(0,1)))
    #checkmate::assertFlag(log_level)

    ## backtransformation if necessary
    # if (log_level)
    Ri <- 2^RiLog
    ## multiply the delta values (biadjacency matrix M) with their weights Ci
    W <- sweep(M, MARGIN = 2, Ci, "*")
    ## sum of the weights per peptide
    W_sum <- rowSums(W)
    ## divide the weights by the sum of the weights per peptide
    W <- as.matrix(sweep(W, 1, W_sum, "/"))
    ## multiply Ri with the corresponding weight
    res_Mat <- sweep(W, MARGIN = 2, Ri, '*')

    ## error term per peptide (on log-scale)
    res_equ <- rjLog - log2(rowSums(res_Mat))

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
.initializeRi <- function(M, rjLog, m) {
    RiLog_start <- rep(NA, m)
    for (j in 1:m) {
        tmp <- M * rjLog
        tmp[tmp == 0] <- NA    ## 0 -> peptide is not present in the protein
        uniquePep <- (rowSums(M) == 1) & (M[, j] == 1)
        if (any(uniquePep)) {
            RiLog_start[j] <- mean(tmp[uniquePep, j], na.rm = TRUE) # .geomMean
        } else {
            RiLog_start[j] <- mean(tmp[, j], na.rm = TRUE) # .geomMean
        }
    }
    #if (log_level) Ri_start <- log2(Ri_start)
    return(RiLog_start)
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
.calcConstraints <- function(fixed.Ci, m) {
    if (!is.null(fixed.Ci)) {
        m2 <- m - sum(!is.na(fixed.Ci)) # number of free weights (not fixed)
        fixed.Ci.sum <- sum(fixed.Ci, na.rm = TRUE)
        eqfun <- function(x) sum(x[(m + 1):(m + m2)]) + fixed.Ci.sum - 1
        # LB <- rep(0, m + m2)
        # if (log_level)
        LB <- c(rep(-Inf, m), rep(0, m2))
        eqB <- 0
    } else {
        eqfun <- function(x) sum(x[(m + 1):(2 * m)]) - 1
        # LB <- rep(0, 2 * m)
        # if (log_level)
        LB <- c(rep(-Inf, m), rep(0, m))
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
.calcObjectiveFunction <- function(fixed.Ci, m, M, rjLog) {
    ### TODO: kann ich m ausrechnen aus M oder fixed.Ci??
    if (!is.null(fixed.Ci)) {
        m2 <- m - sum(!is.na(fixed.Ci)) # number of free weights (not fixed)
        fun <- function(x) {
            RiLog_tmp <- x[1:m]
            Ci_tmp <- fixed.Ci
            Ci_tmp[is.na(fixed.Ci)] <- x[(m + 1):(m + m2)]
            res <- .errorEquation(RiLog = RiLog_tmp, Ci = Ci_tmp, M = M,
                                  rjLog = rjLog)$res_squ_err
            return(res)
        }
    } else {
        fun <- function(x) {
            RiLog_tmp <- x[1:m]
            Ci_tmp <- x[(m + 1):(2 * m)]
            res <- .errorEquation(RiLog = RiLog_tmp, Ci = Ci_tmp, M = M,
                                  rjLog = rjLog)$res_squ_err
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
## TODO: some parameter checks for this private function useful?
.minimizeSquaredError <- function(S,
    fixed.Ci = NULL,
    verbose = FALSE,
    #log_level = TRUE,
    control = list(),
    ...) {
    M <- S$X  ## biadjacency matrix
    m <- ncol(S$X) ## number of proteins
    n <- nrow(S$X) ## number of peptides
    rjLog <- S$fc  ## given peptide ratios (already log2-transformed)
    checkmate::assertMatrix(M, mode = "integerish")
    stopifnot(all(S$X %in% c(0,1)))
    checkmate::assertNumeric(rjLog)
    checkmate::assertNumeric(fixed.Ci, len = length(rjLog), lower = 0, upper = 1,
                             null.ok = TRUE)
    stopifnot(sum(fixed.Ci, na.rm = TRUE) <= 1)
    if (!verbose) control <- c(control, trace = 0)
    is.Ci.fixed <- !is.null(fixed.Ci)
    if (is.Ci.fixed) which.Ci.fixed <- which(!is.na(fixed.Ci))

    Ci_start <- .initializeCi(fixed.Ci, m)
    RiLog_start <- .initializeRi(M, rjLog, m)
    if (is.Ci.fixed) {
        pars <- c(RiLog_start, Ci_start[-which.Ci.fixed])
    } else {
        pars <- c(RiLog_start, Ci_start)
    }
    ## initial error term
    RES <- .errorEquation(RiLog = RiLog_start, Ci = Ci_start, M = M, rjLog = rjLog)

    ### TODO: do we need tracking? -> yes, but maybe separate function?
    track_colnames <- c("iter", "squ_err", paste0("RLog", 1:m), paste0("C", 1:m))
    Tracking <- matrix(c(0, RES$res_squ_err, RiLog_start, Ci_start), nrow = 1)
    Tracking <- as.data.frame(Tracking)
    colnames(Tracking) <- track_colnames

    fun <- .calcObjectiveFunction(fixed.Ci, m, M, rjLog)
    constr <- .calcConstraints(fixed.Ci, m)

    res <- Rsolnp::solnp(pars = pars, fun = fun, LB = constr$LB,
                         eqfun = constr$eqfun, eqB = constr$eqB,
                         control = control)
    # extract optimal Ri and Ci values from optimization result
    RiLog <- res$pars[1:m]
    if (is.Ci.fixed) {
        m2 <- m - sum(!is.na(fixed.Ci)) # number of free weights (not fixed)
        Ci_tmp <- res$pars[(m + 1):(m + m2)]
        Ci <- fixed.Ci
        Ci[is.na(Ci)] <- Ci_tmp
    } else {
        Ci <- res$pars[(m + 1):(2 * m)]
    }
    ## update RES
    RES <- .errorEquation(RiLog = RiLog, Ci = Ci, M = M, rjLog = rjLog)
    Tracking <- rbind(Tracking, c(1, RES$res_squ_err, RiLog, Ci))
    #if (log_level) Ri <- 2^Ri
    result <- list(RiLog = RiLog, Ci = Ci, RES = RES, Tracking = Tracking,
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
    control = list(),
    extend_grid_at_borders = FALSE ) {# ,
    checkmate::assertList(S, types = c("matrix", "numeric"))
    checkmate::assertIntegerish(grid.size, lower = 1)
    checkmate::assertFlag(omit_grid_borders)
    checkmate::assertNumeric(grid.start, lower = 0, upper = 1)
    checkmate::assertNumeric(grid.stop, lower = 0, upper = 1)
    checkmate::assertFlag(verbose)
    checkmate::assertList(control)

    n <- ncol(S$X) ## number of protein groups

    ### TODO: fall mit 1 Protein
    #### Grid überspringen und direkt einzelne Lösung ausspucken
    #### bzw. Grid mit nur einer Zeile


    grid <- seq(grid.start, grid.stop, length.out = grid.size + 1)
    if (extend_grid_at_borders) {
        grid_min <- grid[2] ## 2. element, as first is 0
        grid_max <- grid[length(grid) - 1] ## 2nd to last, as last element is 1
        grid_extend_min <- seq(grid.start, grid_min, length.out = 11)
        grid_extend_max <- seq(grid_max, grid.stop, length.out = 11)
        grid <- sort(unique(c(grid, grid_extend_min, grid_extend_max)))
    }
    if (omit_grid_borders) grid <- grid[-c(1, length(grid))]

    cnames <- c(paste0("RLog", 1:n), paste0("C", 1:n))

    f <- function(j, gridpoint, cnames, S) {
        Ci_tmp <- rep(NA, n)
        Ci_tmp[j] <- gridpoint

        RES <- try({
            .minimizeSquaredError(S,
                                  fixed.Ci = Ci_tmp,
                                  verbose = verbose,
                                  control = control)
        })
        if ("try-error" %in% class(RES)) {
            res_Ri_Ci <- rep(NA, length(cnames))
            error <- NA
        } else {
            res_Ri_Ci <- c(RES$RiLog, RES$Ci)
        }
        names(res_Ri_Ci) <- cnames
        error <- RES$RES$res_squ_err
        result <- c(protein = j, grid = gridpoint, res_Ri_Ci, error = error)
        return(result)
    }
    result <- pbapply::pbmapply(FUN = f,
                      j = rep(1:n, each = length(grid)), gridpoint = grid,
                      MoreArgs = list(cnames = cnames, S = S))
    ## TODO: disable progress bar if "verbose = TRUE"
    #     # if (!verbose)### pboptions(type = "none")
    return(as.data.frame(t(result)))
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
    ratioLog_tol = 1e-6) {

    if (S_is_graph) {
        X <- igraph::as_biadjacency_matrix(S)
        fc <- stats::na.omit(igraph::V(S)$pep_ratio)
        S <- list(X = X, fc = fc)  ### TODO: logFC nennen??
    }



    n <- ncol(S$X) ## number of protein groups

    vapply(1:n, FUN = function(x, error_tol, ratioLog_tol) {
        resProt <- res[res$protein == x, c("error", "RLog1", "C1")]
        colnames(resProt) <- c("error", "RLog", "C")
        return(.analyseResultSingleProt(resProt, error_tol = error_tol, ratioLog_tol = ratioLog_tol))
    }, FUN.VALUES = ,
    error_tol = error_tol, ratioLog_tol = ratioLog_tol)



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


resProt <- res2[res2$protein == 1, c("error", "RLog1", "C1")]
colnames(resProt) <- c("error", "RLog", "C")

resProt2 <- res2[res2$protein == 2, c("error", "RLog2", "C2")]
colnames(resProt) <- c("error", "RLog", "C")

# resProt: result from one specific protein


.analyseResultSingleProt(resProt)
.analyseResultSingleProt(resProt2)


.analyseResultSingleProt <- function(resProt, error_tol = 1e-10,
                                     ratioLog_tol = 1e-6) {
    minError <- min(resProt$error)
    indMinError <- which(abs(minError - resProt$error) <= error_tol)

    RLog_tmp <- resProt$RLog[indMinError]
    C_tmp <- resProt$C[indMinError]
    e_tmp <- resProt$error[indMinError]

    # case 1: single point with minimum error
    if (length(indMinError) == 1) {
        res_tmp <- c(RLog_tmp, NA, NA, C_tmp, NA, NA, 3)
    } else { # case 2: error at least partially constant
        if (abs(diff(range((RLog_tmp)))) <= ratioLog_tol) {  # log2????
            res_tmp <- c(mean(RLog_tmp), NA, NA, NA, min(C_tmp), max(C_tmp), NA) # case 2?
            res_tmp[7] <- ifelse(length(indMinError) == nrow(resProt), 1, 4) # all(R[ind_min_tol] == 0) |???
        } else {
            res_tmp <- c(NA, min(RLog_tmp), max(RLog_tmp), NA, min(C_tmp), max(C_tmp), NA) # case 1?
            res_tmp[7] <- ifelse(length(indMinError) == nrow(resProt), 2, 5)
        }
    }

    res_names <- c("RiLog", "RiLog_min", "RiLog_max", "Ci", "Ci_min", "Ci_max", "case")
    names(res_tmp) <- res_names
    return(res_tmp)

}




