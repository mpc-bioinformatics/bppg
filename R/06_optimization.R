# Functions in this file:
# .errorEquation
# .trackingDataFrame
# .initializeCi
# .initializeRi
# .calcConstraints
# .calcObjectiveFunction
# .minimizeSquaredError
# iterateOverCi
# .calcResultGridpoint
# automatedAnalysisIteratedCi
# .analyseResultSingleProt





#' Function to set up the error equations for the optimization problem
#'
#' @param RiLog       \strong{numeric vector} \cr
#'                    Contains the (estimated) protein ratios (log2-scale).
#' @param Ci          \strong{numeric vector} \cr
#'                    Contains the protein weights (estimated, sum up to 1)
#' @param M           \strong{matrix} \cr
#'                    The biadjaceny matrix of the corresponding graphs.
#' @param rjLog       \strong{numeric vector} \cr
#'                    Contains the measured peptide ratios (log2-scale).
#'
#' @return list containing the following elements:
#' \item{res_Mat}{matrix containing the estimated peptide
#'      ratios using the given Ri and Ci}
#' \item{res_equ}{vector of error terms for each peptide}
#' \item{res_squ_err}{sum of squared error terms}
#' \item{W}{internal weight matrix}
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
    rjLog) {
    m <- length(RiLog) ## number of proteins
    n <- length(rjLog) ## number of peptides
    ## backtransformation
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


#' Funktion to format Tracking information into a dataframe with columnnames.
#'
#' @param i     \strong{numeric} \cr
#'              iterrator, for when this tracking was performed
#' @param RES   \strong{list} \cr
#'              The result from [.errorEquation]
#' @param RiLog \strong{numeric vector} \cr
#'              Contains the (estimated) protein ratios (log2-scale).
#' @param Ci    \strong{numeric vector} \cr
#'              Contains the protein weights (estimated, sum up to 1)
#' @return      \strong{data.frame} \cr
#'              combined information in one data.frame.
.trackingDataFrame <- function(i, RES, RiLog, Ci){
    track_colnames <- c("iter", "squ_err",
        paste0("RLog", seq_along(RiLog)), paste0("C", seq_along(Ci)))
    Tracking <- matrix(c(i, RES$res_squ_err, RiLog, Ci), nrow = 1)
    Tracking <- as.data.frame(Tracking)
    colnames(Tracking) <- track_colnames
    return(Tracking)
}



#' Calulate initial values for Ci (protein weights) for optimization
#'
#' @param fixedCi \strong{numeric vector} \cr
#'                    The fixed protein weights, variable weights set as NA.
#'                    Sum of fixed weights must not exceed 1.
#'                    If NULL, all Cis will be considered as variable.
#'                    This argument is needed to fix Ci on a grid point in the
#'                    iterated_Ci function.
#' @param m \strong{integer(1)} \cr
#'                    number of proteins in the respective graph
#' @details
#' If no Ci are fixed, all Ci are initialized with equal weights (1/m).
#' If at least one Ci is fixed, the remaining weight (1 - sum of fixed Ci)
#' is distributed equally among the non-fixed Ci.
#'
#' @returns Ci_start: vector with initial start values for Ci for the
#'                    optimization step
#'
#' @examples
#' # without fixed Ci -> equal starting weights
#' bppg:::.initializeCi(fixedCi = NULL, m = 3)
#'
#' # with fixed Ci -> fixes these and equal weight for the others
#' bppg:::.initializeCi(fixedCi = c(NA, 0.8, NA), m = 3)
.initializeCi <- function(fixedCi, m) {
    isCiFixed <- !is.null(fixedCi)
    whichCiFixed <- which(!is.na(fixedCi))
    if (sum(fixedCi, na.rm = TRUE) > 1) {
        stop("Sum of fixed Ci exceeds 1.")
    }
    ## Initialization of Ci:
    if (!isCiFixed) {
        ## the algorithm starts with equal weights for each protein
        Ci_start <- rep(1 / m, m)
    } else {
        ## if at least one Ci is fixed, the algorithm distributes the remaining
        ## weight equally among the non-fixed proteins
        m2 <- m - length(whichCiFixed)
        ## sum of fixed Ci (as all Ci have to sum up tp 1)
        fixedCiSum <- sum(fixedCi, na.rm = TRUE)
        Ci_start <- fixedCi
        ## starting values for the remaining Ci values
        Ci_start[is.na(Ci_start)] <- (1 - fixedCiSum) / m2
    }
    return(Ci_start)
}




#' Calulate initial values for Ri (protein ratios) for optimization
#'
#' @param M \strong{matrix} \cr
#'                     Biadjacency matrix of the graph.
#' @param rjLog \strong{numeric} \cr
#'                   Vector of measured peptide ratios (log2-scale).
#'
#' @details
#' For each protein, the mean of the peptide ratios (log2-scale) of the
#' peptides belonging to this protein is calculated.
#' If there are unique peptides for this protein, only those are used.
#'
#' @returns RiLog_start: vector with inital start values for Ri for the
#'                    optimization step (log2-scale)
#'
#' @examples
#' M <- matrix(c(1, 1, 1, 0, 1, 1), ncol = 2, byrow = FALSE)
#' rjLog <- log2(c(0.6, 1.2, 1.5))
#' bppg:::.initializeRi(M, rjLog)
#'
.initializeRi <- function(M, rjLog) {
    m <- ncol(M)
    RiLog_start <- rep(NA, m)
    for (j in seq_len(m)) { # for each protein
        belongsToProt <- (M[, j] == 1) # peptides belonging to protein j
        uniquePep <- (rowSums(M) == 1) & (M[, j] == 1)
        if (any(uniquePep)) {
            RiLog_start[j] <- mean(rjLog[uniquePep & belongsToProt],
                na.rm = TRUE)
        } else {
            RiLog_start[j] <- mean(rjLog[belongsToProt], na.rm = TRUE)
        }
    }
    return(RiLog_start)
}




#' Calculate equality and inequality constraints for the optimization step
#'
#' @param fixedCi \strong{numeric vector} \cr
#'                    The fixed protein weights, variable weights set as NA.
#'                    Sum of fixed weights must not exceed 1.
#'                    If NULL, all Cis will be considered as variable.
#'                    This argument is needed to fix Ci on a grid point in the
#'                    iterated_Ci function.
#' @param m \strong{integer(1)} \cr
#'          number of proteins in the respective graph
#'
#' @returns list containing the following elements (see also
#'                  \code{\link[Rsolnp]{solnp}}):
#' \item{eqfun}{function for calculating equality constraints}
#' \item{eqB}{vector of equality bounds for the variables}
#' \item{LB}{lower bound for the variables}
#'
#' @details
#' x is a vector containing first the RiLog values and then the Ci values
#' (length 2*m if no Ci are fixed).
#'
#' For the equality constraint (eqfun with corresponding bound eqB), the
#' difference of the sum of the Ci values and 1 is calculated. The bound is set
#' to 0, i.e. forcing the sum of the Ci to be 1. In case of fixed Ci values, the
#' sum of the fixed and the free Ci values is considered.
#'
#' The lower bound (LB) for the RiLog values is set to -Inf (no constraint)
#' and to 0 for the Ci.
#'
#' @examples
#' # without fixed Ci
#' bppg:::.calcConstraints(fixedCi = NULL, m = 3)
#'
#' # with fixed Ci
#' bppg:::.calcConstraints(fixedCi = c(NA, 0.8, NA), m = 3)
#'
.calcConstraints <- function(fixedCi, m) {
    if (!is.null(fixedCi)) {
        m2 <- m - sum(!is.na(fixedCi)) # number of free weights (not fixed)
        fixedCiSum <- sum(fixedCi, na.rm = TRUE)
        eqfun <- function(x) sum(x[(m + 1):(m + m2)]) + fixedCiSum - 1
        LB <- c(rep(-Inf, m), rep(0, m2))
        eqB <- 0
    } else {
        eqfun <- function(x) sum(x[(m + 1):(2 * m)]) - 1
        LB <- c(rep(-Inf, m), rep(0, m))
        eqB <- 0
    }
    return(list(eqfun = eqfun, eqB = eqB, LB = LB))
}



#' Calculate objective function for the optimization step
#'
#' @param fixedCi    \strong{numeric vector} \cr
#'                    The fixed protein weights, variable weights set as NA.
#'                    Sum of fixed weights must not exceed 1.
#'                    If NULL, all Cis will be considered as variable.
#'                    This argument is needed to fix Ci on a grid point in the
#'                    iterated_Ci function.
#' @param M           \strong{matrix} \cr
#'                    Biadjacency matrix of the graph.
#' @param rjLog          \strong{numeric vector} \cr
#'                    Contains the measured peptide ratios (log2-scale).
#'
#' @returns objective function that will be minimized
#'
#' @details
#' The objective function is the sum of squared error terms calculated in
#' \code{\link[bppg]{.errorEquation}}.
#'
#' x is a vector containing first the RiLog values and then the Ci values
#' (length 2*m if no Ci are fixed). If there are fixed Ci and m2 non-fixed ones,
#' x has the length m + m2.
#'
#' @examples
#'
#' M <- matrix(c(1, 1, 1, 0, 1, 1), ncol = 2, byrow = FALSE)
#' rjLog <- log2(c(0.6, 1.2, 1.5))
#'
#' fun <- bppg:::.calcObjectiveFunction(fixedCi = NULL, M = M, rjLog = rjLog)
#'
#'
.calcObjectiveFunction <- function(fixedCi, M, rjLog) {
    m <- ncol(M)
    if (!is.null(fixedCi)) {
        m2 <- m - sum(!is.na(fixedCi)) # number of free Ci (not fixed)
        fun <- function(x) {
            RiLog_tmp <- x[seq_len(m)]
            Ci_tmp <- fixedCi
            Ci_tmp[is.na(fixedCi)] <- x[(m + 1):(m + m2)]
            res <- .errorEquation(RiLog = RiLog_tmp, Ci = Ci_tmp, M = M,
                rjLog = rjLog)$res_squ_err
            return(res)
        }
    } else {
        fun <- function(x) {
            RiLog_tmp <- x[seq_len(m)]
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
#' @param G           \strong{igraph object} \cr
#'                    bipartite peptide-protein graph
#' @param fixedCi    \strong{numeric vector} \cr
#'                    The fixed protein weights, variable weights set as NA.
#'                    Sum of fixed weights must not exceed 1.
#'                    If NULL, all Cis will be considered as variable.
#'                    This argument is needed to fix Ci on a grid point in the
#'                    iterated_Ci function.
#' @param verbose     \strong{logical} \cr
#'                    If \code{TRUE}, additional information on each iteration
#'                    of the optimization is printed
#'                    (see \code{\link[Rsolnp]{solnp}}).
#' @param control     \strong{list} \cr
#'                    The control parameters for solnp
#'                    (see \code{\link[Rsolnp]{solnp}}).
#'
#' @return list containing the following elements:
#' \item{RiLog}{estimated protein ratios (Log2-scale)}
#' \item{Ci}{estimated protein weights}
#' \item{RES}{final result of \code{\link[bppg]{.errorEquation}}, which also
#' contains the final, minimal error term}
#' \item{Tracking}{Tracking of Ri, Ci and error term for the
#'  different iterations (-1: no optimization)}
#' \item{outer.iter}{Number of outer iterations needed for the optimization
#'  algorithm to converge or stop (see also \code{\link[Rsolnp]{solnp}})}
#' \item{convergence}{Indicates whether the solver has converged (0) or
#' not (1 or 2) (see also \code{\link[Rsolnp]{solnp}})}
#'
#' @examples
#' file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
#' graphs <- readRDS(file)
#' G <- graphs$"1_2"[[2]]
#' bppg:::.minimizeSquaredError(G)
#'
#' @importFrom igraph as_biadjacency_matrix vertex_attr
#' @importFrom stats na.omit
#' @importFrom checkmate assertNumeric
#' @importFrom Rsolnp solnp
.minimizeSquaredError <- function(G,
    fixedCi = NULL,
    verbose = FALSE,
    control = list()
    ) {
    M <- igraph::as_biadjacency_matrix(G)
    m <- ncol(M) ## number of proteins
    n <- nrow(M) ## number of peptides
    rjLog <- stats::na.omit(igraph::vertex_attr(G, "pep_logRatio"))
    if (is.null(rjLog)) stop("G does not contain peptide ratios.")
    checkmate::assertNumeric(rjLog)
    if (m == 1){
        RiLog <- mean(rjLog, na.rm = TRUE)
        RES <- .errorEquation(RiLog = c(RiLog),
            Ci = c(1.0), M = M, rjLog = rjLog)
        Tracking <- .trackingDataFrame(-1, RES, c(RiLog), c(1.0))
        result <- list(RiLog = c(RiLog), Ci = 1.0, RES = RES,
            Tracking = Tracking, outer.iter = 0, convergence = 0)
        return(result)
    }

    checkmate::assertNumeric(fixedCi, len = m, lower = 0, upper = 1,
        null.ok = TRUE)
    stopifnot(sum(fixedCi, na.rm = TRUE) <= 1)
    if (!verbose) control <- c(control, trace = 0)
    isCiFixed <- !is.null(fixedCi)
    Ci_start <- .initializeCi(fixedCi, m)
    RiLog_start <- .initializeRi(M, rjLog)
    if (isCiFixed) {
        whichCiFixed <- which(!is.na(fixedCi))
        pars <- c(RiLog_start, Ci_start[-whichCiFixed])
    } else {
        pars <- c(RiLog_start, Ci_start)
    }
    ## initial error term
    RES <- .errorEquation(RiLog = RiLog_start, Ci = Ci_start, M = M,
        rjLog = rjLog)
    Tracking <- .trackingDataFrame(0, RES, RiLog_start, Ci_start)

    fun <- .calcObjectiveFunction(fixedCi, M, rjLog)
    constr <- .calcConstraints(fixedCi, m)
    res <- Rsolnp::solnp(pars = pars, fun = fun, LB = constr$LB,
        eqfun = constr$eqfun, eqB = constr$eqB, control = control)
    # extract optimal Ri and Ci values from optimization result
    RiLog <- res$pars[seq_len(m)]
    if (isCiFixed) {
        m2 <- m - sum(!is.na(fixedCi)) # number of free weights (not fixed)
        Ci_tmp <- res$pars[(m + 1):(m + m2)]
        Ci <- fixedCi
        Ci[is.na(Ci)] <- Ci_tmp
    } else {
        Ci <- res$pars[(m + 1):(2 * m)]
    }
    ## update RES
    RES <- .errorEquation(RiLog = RiLog, Ci = Ci, M = M, rjLog = rjLog)
    Tracking <- rbind(Tracking, c(1, RES$res_squ_err, RiLog, Ci))
    result <- list(RiLog = RiLog, Ci = Ci, RES = RES, Tracking = Tracking,
        outer.iter = res$outer.iter, convergence = res$convergence)
    return(result)
}



#' Calculate optimal solution for a single gridpoint
#'
#' @param j             \strong{integer(1)} \cr
#'                      index of protein
#' @param gridpoint     \strong{numeric(1)} \cr
#'                      Ci value to fix for this protein
#' @param cnames        \strong{character} \cr
#'                      column names for result dataframe
#' @param G             \strong{igraph object} \cr
#'                      bipartite peptide-protein graph
#' @param n             \strong{integer(1)} \cr
#'                      number of proteins in the graph
#' @param ...           additional arguments for [.minimizeSquaredError()],
#'                     e.g. verbose, control
#'
#' @returns Vector with the following elements:
#' \item{protein}{index of protein (j)}
#' \item{grid}{gridpoint}
#' \item{res_Ri_Ci}{multiple values: estimated RiLog and Ci values}
#' \item{error}{error term}
.calcResultGridpoint <- function(j, gridpoint, cnames, G, n, ...) {
    Ci_tmp <- rep(NA, n)
    Ci_tmp[j] <- gridpoint

    RES <- try({
        .minimizeSquaredError(G, fixedCi = Ci_tmp, ...)
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



#' Iterate over possible Ci values
#'
#' @param G                        \strong{igraph object} \cr
#'                                 Bipartite peptide-protein graph
#' @param gridStart               \strong{integer} \cr
#'                                 The start of the grid (default is 0).
#' @param gridStop                \strong{integer} \cr
#'                                 The end of the grid (default is 1).
#' @param gridSize                \strong{integer} \cr
#'                                 The number of grid points for the Cis.
#' @param extend_grid_at_borders   \strong{logical} \cr
#'                                 If \code{TRUE}, the grid will be extend close
#'                                 to the borders (0 and 1).
#'                                 While exact values of 0 and 1 may cause
#'                                 numerical problems, values close to those may
#'                                 be valuable to get a better estimate of the
#'                                 protein ratios. Default is \code{FALSE}.
#' @param omit_grid_borders        \strong{logical} \cr
#'                                 If \code{TRUE} (default), omit exact value of
#'                                 1 and 0 from the grid (recommended, as they
#'                                 may cause numerical issues).
#' @param verbose                  \strong{logical} \cr
#'                                 If \code{TRUE}, print additional information
#'                                 (see \code{\link[Rsolnp]{solnp}} function).
#' @param control                  \strong{list} \cr
#'                                 The \code{control} object to be passed to the
#'                                 [.minimizeSquaredError()] function.
#'
#' @return
#' A data.frame containing the optimal Ci and Ri values together with the reached
#' minimal error term for each grid point.
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
#' The table can be further processed with
#' [bppg::automatedAnalysisIteratedCi()].
#'
#' @examples
#' file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
#' graphs <- readRDS(file)
#' G <- graphs$"1_2"[[2]]
#' # small example with a small grid size
#' iterateOverCi(G, gridSize = 100)
#'
#' @importFrom checkmate assertClass assertFlag assertIntegerish assertList assertNumeric checkTRUE
#' @importFrom igraph is_bipartite V
#' @importFrom  pbapply pbmapply pboptions
iterateOverCi <- function(G,
    gridStart = 0,
    gridStop = 1,
    gridSize = 1000,
    omit_grid_borders = TRUE,
    extend_grid_at_borders = FALSE,
    verbose = FALSE,
    control = list()) {
    checkmate::assertClass(G, classes = c("igraph"))
    checkmate::checkTRUE(igraph::is_bipartite(G))
    checkmate::assertIntegerish(gridSize, lower = 1)
    checkmate::assertFlag(omit_grid_borders)
    checkmate::assertFlag(extend_grid_at_borders)
    checkmate::assertNumeric(gridStart, lower = 0, upper = 1)
    checkmate::assertNumeric(gridStop, lower = 0, upper = 1)
    checkmate::assertFlag(verbose)
    checkmate::assertList(control)

    n <- sum(igraph::V(G)$type) ## number of protein groups
    if (n == 1) { # special case for only a single protein group in the graph
        RES <-  .minimizeSquaredError(G, fixedCi = NULL, verbose = verbose,
            control = control)
        result <- data.frame(protein = 1, grid = 1, RLog1 = RES$RiLog,
            C1 = 1, error = RES$RES$res_squ_err)
        return(result)
    } else { # n > 1
        grid <- seq(gridStart, gridStop, length.out = gridSize + 1)
        if (extend_grid_at_borders) {
            grid_min <- grid[2] ## 2nd element, as first is 0
            grid_max <- grid[length(grid) - 1] ## 2nd to last, last element is 1
            grid_extend_min <- seq(gridStart, grid_min, length.out = 11)
            grid_extend_max <- seq(grid_max, gridStop, length.out = 11)
            grid <- sort(unique(c(grid, grid_extend_min, grid_extend_max)))
        }
        if (omit_grid_borders) grid <- grid[-c(1, length(grid))]
        cnames <- c(paste0("RLog", seq_len(n)), paste0("C", seq_len(n)))
        if (!verbose) {
            pbo <- pbapply::pboptions(type = "none")
            on.exit(pbapply::pboptions(pbo), add = TRUE)
        }
        result <- pbapply::pbmapply(FUN = .calcResultGridpoint,
            j = rep(seq_len(n), each = length(grid)), gridpoint = grid,
                MoreArgs = list(cnames = cnames, G = G, n = n,
                    verbose = verbose, control = control))
        return(as.data.frame(t(result)))
    }
}




#' Analyse results for a single protein
#'
#' @param protNr \strong{integer(1)} \cr
#'                                      Index of the protein to be analysed.
#' @param resProt \strong{data.frame} \cr
#'                                      The data.frame resulting from the
#'                                      [bppg::iterateOverCi()] function,
#'                                      filtered for a specific protein.
#' @param error_tol \strong{numeric(1)} \cr
#'                  tolerance for the error term.
#' @param ratioLog_tol \strong{numeric(1)} \cr
#'                      tolerance for the log protein ratios.
#'
#' @returns Vector with minimal error, estimate for Ri (single value
#' or min/max), estimate for Ci (single value or min/max).
#'
.analyseResultSingleProt <- function(protNr, resProt,
    error_tol = 1e-10, ratioLog_tol = 1e-6) {
    proteins <- unique(resProt$protein)
    minError <- min(resProt$error)
    indMinError <- which(abs(minError - resProt$error) <= error_tol)

    RLog_tmp <- resProt$RLog[indMinError]
    C_tmp <- resProt$C[indMinError]
    e_tmp <- resProt$error[indMinError]

    if (length(indMinError) == 1) {
        # case 3: single solution with one optimal data point
        res_tmp <- c(minError, RLog_tmp, NA, NA, C_tmp, NA, NA, 3)
    } else {
        if (abs(diff(range((RLog_tmp)))) <= ratioLog_tol) {
            # case 1 and 4: RLog is almost constant
            RLog_tmp_mean_rounded <- round(mean(RLog_tmp),
                digits = ceiling(-log10(ratioLog_tol)))
            res_tmp <- c(mean(e_tmp), RLog_tmp_mean_rounded, NA, NA, NA,
                min(C_tmp), max(C_tmp), NA)
            # all(R[ind_min_tol] == 0) |???
            res_tmp[8] <- ifelse(length(indMinError) == nrow(resProt), 1, 4)
        } else {
            # case 2 and 5: range solution for Rlog
            res_tmp <- c(mean(e_tmp), NA, min(RLog_tmp), max(RLog_tmp), NA,
                min(C_tmp), max(C_tmp), NA)
            res_tmp[8] <- ifelse(length(indMinError) == nrow(resProt), 2, 5)
        }
    }

    res_names <- c("error_min", "RiLog", "RiLog_min", "RiLog_max",
        "Ci", "Ci_min", "Ci_max", "case")
    names(res_tmp) <- res_names
    return(res_tmp)

}








#' Extract protein ratio solutions from the result of iterateOverCi()
#'
#' @param G                             \strong{igraph object} \cr
#'                                      An igraph graph of the bipartite
#'                                      peptide-protein graph with peptide
#'                                      ratios
#' @param res                           \strong{data.frame} \cr
#'                                      The data.frame resulting from the
#'                                      [bppg::iterateOverCi()] function.
#' @param use_results_from_other_proteins   \strong{logical} \cr
#'                                          If \code{TRUE}, the results from
#'                                          other proteins within the same graph
#'                                          will be used to calculate the
#'                                          optimal solution for each protein
#'                                          node. Default is TRUE
#' @param verbose                       \strong{logical} \cr
#'                                      If \code{TRUE}, additional information
#'                                      will be printed.
#' @param job                           \strong{BatchExperiment job object} \cr
#'                                      Is used to print the job id and
#'                                      parameters in the output
#' @param error_tol                    \strong{numeric(1)} \cr
#'                                      Tolerance for a constant error term.
#'                                      Default is 1e-10.
#' @param ratioLog_tol                 \strong{numeric(1)} \cr
#'                                      Tolerance for a constant log2 protein
#'                                      ratios. The default is 1e-6.
#'
#' @return A SummarizedExperiment object with one row for each protein in the
#'  assay.
#' @export
#'
#' @seealso [bppg::iterateOverCi()], [.minimizeSquaredError()]
#'
#' @examples
#' file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
#' graphs <- readRDS(file)
#' G <- graphs$"1_2"[[2]]
#' # small example with a small grid size
#' res <- iterateOverCi(G, gridSize = 100)
#' automatedAnalysisIteratedCi(G, res)
automatedAnalysisIteratedCi <- function(G,
                                        res,
                                        use_results_from_other_proteins = TRUE,
                                        verbose = FALSE,
                                        job = NULL,
                                        error_tol = 1e-10,
                                        ratioLog_tol = 1e-6) {

    n <- sum(igraph::V(G)$type) ## number of protein groups
    accessions <- igraph::V(G)$name[igraph::V(G)$type]

    if (!is.null(job)) {
        graphID <- job$pars$prob.pars$k
        comparison <- job$prob.name
        job.id <- job$job.id
    } else {
        graphID <- NA
        comparison <- NA
        job.id <- NA
    }

    # filter res for potential NaNs in the error column (could occur if a Ci
    # is estimated as 0)
    ind_error_NA <- which(is.na(res$error))
    if (length(ind_error_NA) > 0) {
        if (verbose) {
            message(paste0(length(ind_error_NA),
                           " grid points with NA or NaN error term were removed."))
        }
        res <- res[-ind_error_NA, ]
    }

    f <- function(x, res, error_tol, ratioLog_tol,
        use_results_from_other_proteins) {

        cols <- c("protein", "error", paste0("RLog", x), paste0("C", x))

        resProt <- res[res$protein == x, cols]
        colnames(resProt) <- c("protein", "error", "RLog", "C")

        if (use_results_from_other_proteins) {
            resProt2 <- res[res$protein != x, cols]
            colnames(resProt2) <- c("protein", "error", "RLog", "C")
            # remove too extreme Ci
            resProt2 <- resProt2[resProt2$C > 0.01 & resProt2$C < 0.99,]
            resProt <- rbind(resProt, resProt2)
        }
        return(.analyseResultSingleProt(x, resProt, error_tol = error_tol,
            ratioLog_tol = ratioLog_tol))
    }

    RES <- vapply(seq_len(n), FUN = f,
        FUN.VALUE = c("error_min" = 0, "RiLog" = 0, "RiLog_min" = 0,
        "RiLog_max" = 0, "Ci" = 0, "Ci_min" = 0, "Ci_max" = 0, "case" = 0),
        res = res, error_tol = error_tol, ratioLog_tol = ratioLog_tol,
        use_results_from_other_proteins = use_results_from_other_proteins)

    # comparison=rep(comparison, n),
    RES_info <- data.frame(accession = accessions,
        graphID = rep(graphID, n), proteinNr = seq_len(n))

    RES <- cbind(RES_info, as.data.frame(t(RES)))
    rownames(RES) <- RES$accession
    RES <- RES[, -1] # remove accession, as it is now in rownames

    RES_SE <- SummarizedExperiment::SummarizedExperiment(
        assays = list(results = RES),
        rowData = data.frame(accession = accessions),
        colData = data.frame(colnames = colnames(RES)))
    return(RES_SE)
}








