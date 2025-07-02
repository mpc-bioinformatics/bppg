
#' Function to set up the equations for the optimization problem
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
#'                    If \code{TRUE}, the Ri are given on log2-level and need to be back-transformed here (this may allow a symmetric behaviour during optimization)
#'
#' @return list containing the following elements:
#' \item{res_Mat}{matrix containing the estimated peptide ratios using Ri and Ci}
#' \item{res_equ}{vector of error terms for each peptide}
#' \item{res_squ_err}{sum of squared error terms}
#' \item{W}{internal weight matrix}
#'
#' @export
#'
#' @examples
#' Ri <- c(0.5, 1.3)
#' Ci <- c(0.3, 0.7)
#' M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
#' rj <- c(0.6, 1.2)
#' equation(Ri, Ci, M, rj)
equation <- function(Ri,
                     Ci,
                     M,
                     rj,
                     #error.type = "multiplicative",
                     #error.trans = "square",
                     log_level = FALSE) {

  m <- length(Ri) # number of proteins
  n <- length(rj) # number of peptides

  if (log_level) Ri <- 2^Ri   ## backtransformation if necessary
  W <- sweep(M, MARGIN = 2, Ci, '*')  ## multiply the delta values (from the biadjacency matrix) with the corresponding weights Ci
  W_sum <- rowSums(W) # sum of the weights per peptide
  W <- as.matrix(sweep(W, 1, W_sum, "/")) # divide the weights by the sum of the weights per peptide

  res_Mat <- sweep(W, MARGIN = 2, Ri, '*') # multiply Ri with the corresponding weight

  # error term per peptide (on log-scale)
  res_equ <- log(rj) - log(rowSums(res_Mat))

  # sum of squared error terms
  res_squ_err <- sum(res_equ^2)

  return(list(res_Mat = res_Mat, res_equ = res_equ, res_squ_err = res_squ_err, W = W))
}




#' Function to set up the optimization problem and minimize the sum of squared error terms
#'
#' @param S           \strong{list} \cr
#'                    A list of biadjacency matrix of the bipartite peptide-protein graph (X) and measured peptide ratios (fc).
#' @param fixed.Ci    \strong{numeric vector} \cr
#'                    The fixed protein weights, variable weights set as NA. Sum of fixed weights must not exceed 1.
#'                    If NULL, all Cis will be considered as variable.
#'                    This argument is needed to fix Ci on a grid point in the iterated_Ci function.
#' @param verbose     \strong{logical} \cr
#'                    If \code{TRUE}, additional information on each iteration of the optimization is printed (see also rsolnp function in package Rsolnp).
#' @param reciprocal  \strong{logical} \cr
#'                    If \code{TRUE}, the reciprocal of the peptide ratios is used for the optimization.
#' @param log_level   \strong{logical} \cr
#'                    If \code{TRUE}, the Ri are log2-transformed before optimization, allowing a symmetric consideration of Ri < 0 and > 0.
#' @param control     \strong{list} \cr
#'                    The control parameters for solnp.
#' @param min_ci      \strong{float} \cr
#'                    Lower bound for concentration
#' @param ...         Additional parameters to solnp.
#'
#' @return list containing the following elements:
#' \item{Ri}{estimated protein ratios}
#' \item{Ci}{estimated protein weights}
#' \item{RES}{final result of equation(), which also contains the final, minimal error term}
#' \item{Tracking}{Tracking of Ri, Ci and error term for the different iterations}
#' \item{outer.iter}{Number of outer iterations needed for the optimization algorithm to converge or stop}
#' \item{convergence}{Indicates whether the solver has converged (0) or not (1 or 2).}
#'
#' @export
#'
#' @examples
#' M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
#' rj <- c(0.6, 1.2)
#' S <- list(X = M, fc = rj)
#' minimize_squared_error(S)
#'
minimize_squared_error <- function(S,
                                      #error.type = "multiplicative",
                                      fixed.Ci = NULL,
                                      verbose = FALSE,
                                      #error.trans = "square",
                                      reciprocal = FALSE,
                                      log_level = TRUE,
                                      control = list(),
                                      min_ci = 0.0001,
                                      ...) {

  is.Ci.fixed <- !is.null(fixed.Ci)
  if (is.Ci.fixed) which.Ci.fixed <- which(!is.na(fixed.Ci)) # assesses which Cis are fixed by the user
  if (sum(fixed.Ci, na.rm = TRUE) > 1) stop("Sum of chosen Ci values exceeds 1!")

  m <- ncol(S$X) # number of proteins
  n <- nrow(S$X) # number of peptides
  rj <- S$fc     # given peptide ratios
  if (reciprocal) {
    rj <- 1/rj
  }
  M <- S$X  # biadjacency matrix

  ### Initialization of Ci:
  if (is.Ci.fixed) {
    # if at least one Ci is fixed, the algorithm distributes the remaining weight equally among the non-fixed proteins
    m2 <- m - length(which.Ci.fixed)             # number of non-fixed Ci values
    fixed.Ci.sum <- sum(fixed.Ci, na.rm = TRUE)  # sum of fixed Ci (as all Ci have to sum up tp 1)
    Ci_start <- fixed.Ci
    Ci_start[is.na(Ci_start)] <- (1-fixed.Ci.sum)/m2 # starting values for the remaining Ci values
  } else {
    Ci_start <- rep(1/m, m) # if no Ci is fixed, the algorithm starts with equal weights for each protein
  }

  ### initialization of Ri as the geometric mean of (if possible only unique) peptide ratios
  Ri <- rep(NA, m)
  for (j in 1:m) {
    tmp <- S$X * S$fc
    tmp[tmp == 0] <- NA    # 0 -> peptide is not present in the protein
    unique <- (rowSums(M)==1) & (M[,j] == 1)
    if (any(unique)) {
      Ri[j] <- 2^mean(log2(tmp[unique,j]), na.rm = TRUE)
    } else {
      Ri[j] <- 2^mean(log2(tmp[,j]), na.rm = TRUE)
    }
  }
  if (log_level) Ri <- log2(Ri)

  ### initial error term
  RES <- equation(Ri = Ri, Ci = Ci_start, M = M, rj = rj,
                  log_level = log_level)

  track_colnames <- c("iter", "squ_err", paste0("R", 1:m), paste0("C", 1:m))
  Tracking <- matrix(c(0, RES$res_squ_err,Ri, Ci_start), nrow = 1)
  Tracking <- as.data.frame(Tracking)
  colnames(Tracking) <- track_colnames

  iter <- 1

  ## starting parameters
  if (is.Ci.fixed) {
    pars <- c(Ri, Ci_start[-which.Ci.fixed])
  } else {
    pars <- c(Ri, Ci_start)
  }

  ### function to optimize
  if (is.Ci.fixed) {
    fun <- function(x) {
      Ri_tmp <- x[1:m]
      Ci_tmp <- fixed.Ci
      Ci_tmp[is.na(fixed.Ci)] <- x[(m+1):(m+m2)]
      equation(Ri = Ri_tmp, Ci = Ci_tmp, M = M, rj = rj,
                  #error.type = error.type, error.trans = error.trans,
                  log_level = log_level)$res_squ_err
    }
  } else {
    fun <- function(x) {
      Ri_tmp <- x[1:m]
      Ci_tmp <- x[(m+1):(2*m)]
      equation(Ri = Ri_tmp, Ci = Ci_tmp, M = M, rj = rj,
                  #error.type = error.type, error.trans = error.trans,
                  log_level = log_level)$res_squ_err
    }
  }

  ### constraints
  if (is.Ci.fixed) {
    eqfun <- function(x) sum(x[(m + 1):(m + m2)]) + fixed.Ci.sum - 1  ## sum Ci = 1
    eqB <- 0 
    LB <- c(rep(0, m), rep(min_ci / 10, m2)) # min value in grid for Ci
    if (log_level) LB <- c(rep(-Inf, m), rep(min_ci / 10, m2))
  } else {
    eqfun <- function(x) sum(x[(m + 1):(2 * m)]) - 1  ## sum Ci = 1
    eqB <- 0
    LB <- c(rep(0, m), rep(min_ci / 10, m)) # min value in grid for Ci
    if (log_level) LB <- c(rep(-Inf, m), rep(min_ci / 10, m))
  }


  if (verbose) {
    control <- c(control, trace = 0)
  }


  ### Optimization
  #res <- Rsolnp::solnp(pars = pars, fun = fun, LB = LB, eqfun = eqfun, eqB = eqB, control = control, ...)
  # now in c implementation
  res <- Rsolnp::csolnp(pars = pars, fn = fun, lower = LB, eq_fn = eqfun, eq_b = eqB, control = control, ...)


  outer.iter <- res$outer.iter
  convergence <- res$convergence

  Ri <- res$pars[1:m]

  if (is.Ci.fixed) {
    Ci_tmp <- res$pars[(m+1):(m+m2)]
    Ci <- fixed.Ci
    Ci[is.na(Ci)] <- Ci_tmp
  } else {
    Ci <- res$pars[(m+1):(2*m)]
  }

  ## update RES
  RES <- equation(Ri = Ri, Ci = Ci, M = M, rj = rj,
                     #error.type = error.type, error.trans = error.trans,
                     log_level = log_level)

  Tracking <- rbind(Tracking, c(iter, RES$res_squ_err, Ri, Ci))


  if (log_level) {
    Ri <- 2^Ri
  }
  if (reciprocal) {
    Ri <- 1/Ri
  }


  result <- list(Ri = Ri, Ci = Ci, RES = RES, Tracking = Tracking, outer.iter = outer.iter, convergence = convergence)
  class(result) <- "res_min_squ_error"

  return(result)

}
