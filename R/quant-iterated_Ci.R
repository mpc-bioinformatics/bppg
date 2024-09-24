

### TODO: progress bar!
### TODO: use graphs instead of submatrix

#' Iterate over possible Ci values
#'
#' @param S list of biadjacency matrix (M) and
#' @param grid.size number of grid points for the Cis
#' @param omit_grid_borders if TRUE, omit exact value of 1 and 0 from the grid (recommended, as they may cause numerical issues)
#' @param grid.start start of the grid (default is 0)
#' @param grid.stop end of the grid (default is 1)
#' @param verbose if TRUE, print additional information (see solnp function)
#' @param verbose_opt verbose argument of the minimize_squared_error function
#' @param control control object to be passed to minimize_squared_error
#' @param extend_grid_at_borders extend the grid close to the borders (0 and 1). While exact values of 0 and 1 may cause numerical problems, values close to those may be valuable to get a better estimate of the protein ratios
#' @param log_level passed to minimize_squared_error
#'
#' @return
#' A dataframe containing the optimal Ci and Ri values together with the reached minimal error term.
#'
#' @export
#'
#' @details
#' With minimize_squared_error() each protein node receives one estimate for the protein ratio. However, in some cases, there are multiple possible values for the protein ratios that lead to the same, minimal error term.
#' To get a better coverage of the optimal solutions, this function iterates over a grid of possible weights Ci for each protein node.
#' For each protein node, the Ci value if fixed on a point on the grid, while the other Ci values and all Ri values are optimized using minimize_squared_error().
#' The resulting table can be used to assess a range of possible solutions for the protein ratios.
#'
#'
#' @examples
#' M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
#' rj <- c(0.6, 1.2)
#' S <- list(X = M, fc = rj)
#' minimize_squared_error(S)
iterate_over_Ci <- function(S,
                            grid.size = 1000,
                            omit_grid_borders = TRUE,
                            grid.start = 0,
                            grid.stop = 1,
                            verbose = FALSE,
                            verbose_opt = FALSE,
                            control = list(),
                            extend_grid_at_borders = FALSE,
                            log_level = TRUE) {

  n <- ncol(S$X)# number of protein groups

  grid <- seq(grid.start, grid.stop, length.out = grid.size+1)

  if (extend_grid_at_borders) {
    grid_min <- grid[2] # 2. Element, da erstes 0
    grid_max <- grid[length(grid)-1] # 1. Element, da letztes 1

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

  pb <- pbapply::startpb(0, n*length(grid))
  for (j in 1:n) {
    for (i in seq_along(grid)) {

      pbapply::setpb(pb, (j-1) * length(grid) + i)

      if (verbose) print(paste0("j = ", j, " i = ", i))
      Ci_tmp <- rep(NA, n)
      Ci_tmp[j] <- grid[i]

      RES <- try({
        minimize_squared_error(S,
                               #error.type = "multiplicative",
                               fixed.Ci = Ci_tmp,
                               verbose = verbose_opt,
                               #error.trans = "square",
                               reciprocal = FALSE,
                               control = control,
                               log_level = log_level)
      })
      if ("try-error" %in% class(RES)) {
        if(grepl("reached elapsed time limit", RES)) {
          stop(paste0("protein ", j, " grid point ", i))
        } else {
          next
        }
      }
      err_tmp[i] <- RES$RES$res_squ_err
      Ris_tmp[i, ] <- RES$Ri
      Cis_tmp[i,] <- RES$Ci

    }

    result_tmp <- data.frame(protein = rep(j, length(grid)), grid = grid, Ris_tmp, Cis_tmp, error = err_tmp)
    result <- rbind(result, result_tmp)
  }

  invisible(NULL)
  pbapply::closepb(pb)

  return(result)
}




#' Extract protein ratio solutions from the result of iterate_over_Ci()
#'
#' @param S igraph object or list of biadjacency matrix (M) and fold changes (fc)
#' @param res result of iterate_over_Ci()
#' @param use_results_from_other_proteins if TRUE, use the results from other proteins within the same graph to calculate the optimal solution for each protein node
#' @param verbose if TRUE, print additional information
#' @param job job object from the BatchExperiment. Used to print the job id and parameters in the output
#' @param S_is_graph set to TRUE,if S is an igraph object and FALSE if S is a list with the biadjacency matrix and fold changes
#'
#' @return data frame
#' @export
#'
#' @examples # TODO
automated_analysis_iterated_Ci <- function(S,
                                           res,
                                           use_results_from_other_proteins = FALSE,
                                           verbose = FALSE,
                                           job = NULL,
                                           S_is_graph = FALSE) {

  if (S_is_graph & !is.null(S)) {
    X <- igraph::as_incidence_matrix(S)
    fc <- stats::na.omit(igraph::V(S)$pep_ratio)
    S <- list(X = X, fc = fc)
  }


  if (!is.null(S)) {
    nr_proteins <- ncol(S$X)
    if (nr_proteins == 1 ) res <- NULL
  } else {
    nr_proteins <- 1
    res <- NULL
  }


  error_optimal <- NA
  Ri_optimal <- rep(NA, nr_proteins)
  Ci_optimal <- rep(NA, nr_proteins)




  if (nr_proteins == 1 & is.null(res)) {
    # graph contains only one protein node and was skipped, so optimal solution has to be calculated here
    solution_optimal <- bppg::minimize_squared_error(S,
                                                     #error.type = "multiplicative",
                                                     fixed.Ci = NULL,
                                                     verbose = FALSE, #error.trans = "square",
                                                     reciprocal = FALSE, log_level = TRUE,
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
      D_tmp$jop.id <- job$job.id
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
      D_tmp$jop.id <- job$job.id
    }

    if (verbose) print(paste0("Protein ", i))

    X_tmp <- res[res$protein == i,]
    X_tmp <- X_tmp[!is.na(X_tmp$error),]
    C_tmp <- X_tmp[, paste0("C", i)]

    ### use results from other proteins but only when Ci is not too extreme
    X_tmp3 <- res[res$protein != i,]
    X_tmp3 <- X_tmp3[!is.na(X_tmp3$error),]
    C_tmp3 <- X_tmp3[, paste0("C", i)]
    X_tmp3 <- X_tmp3[X_tmp3$protein == i | (C_tmp3 >= 0.01 & C_tmp3 < 0.99),]
    R_3 <- X_tmp3[, paste0("R", i)]
    C_3 <- X_tmp3[, paste0("C", i)]


    error <- X_tmp$error
    R <- X_tmp[, paste0("R", i)]
    C <- X_tmp[, paste0("C", i)]

    D_tmp$min_error <- min(error_optimal, min(error), na.rm = TRUE)
    D_tmp$error_optimal <- error_optimal


    #### 1st check: is error constant?
    if (abs(diff(range(X_tmp$error))) <= 1e-10) {
      ### constant error, i.e. single solution or interval with lower and upper border

      D_tmp$error_constant <- "yes"

      if (verbose) print(paste0("Error is nearly constant (abs. diff. = ", abs(diff(range(error))), ")."))
      if (verbose) print(paste0("Mean error is ", mean(error, na.rm = TRUE), "."))
      if (verbose) print(paste0("Minimal error is ", min(error, na.rm = TRUE) , "."))
      ind_min <- which.min(error)

      #### 2nd check: Is Ri constant too?
      if (abs(diff(range(R))) > 1e-4) {
        if (verbose) print(paste0("Range for R", i, ": ", BBmisc::collapse(range(R), sep = " - "), "."))
        D_tmp$Ri <- NA
        D_tmp$Ri_min <- min(R)
        D_tmp$Ri_max <- max(R)
        D_tmp$case <- 1
      } else {
        if (verbose) print(paste0("Constant Solution for R", i, ": ", bppg::geom_mean(R)))
        D_tmp$Ri <- bppg::geom_mean(R)
        D_tmp$Ri_min <- NA
        D_tmp$Ri_max <- NA
        D_tmp$case <- 2
      }
      if (verbose) print(paste0("Range for C", i, ": ", BBmisc::collapse(range(C), sep = " - "), "."))
      D_tmp$Ci_min <- min(C)
      D_tmp$Ci_max <- max(C)

    } else {
      error <- X_tmp$error
      R <- X_tmp[, paste0("R", i)]
      C <- X_tmp[, paste0("C", i)]

      if (verbose & !is.na(error_optimal) & any(error < error_optimal)) {
        warning(paste0("Lower than optimal error detected. Difference = ", abs(min(error) - error_optimal)))
      }

      ind_min <- which.min(error)
      ind_min_tol <- which(abs(min(error) - error) <= 1e-10)  ### area in which the error is almost constant

      if (length(ind_min_tol) <= 1) {

        D_tmp$error_constant <- "no"
        D_tmp$Ri <- R[ind_min_tol]
        D_tmp$Ri_min <- NA
        D_tmp$Ri_max <- NA
        D_tmp$Ci <- C[ind_min_tol]
        D_tmp$Ci_min <- NA
        D_tmp$Ci_max <- NA
        D_tmp$case <- 3

        if (verbose) print(paste0("Single optimal solution with error ", min(error, na.rm = TRUE), "."))

        if (verbose) print(paste0("Single optimal solution is R", i, "= ", Ri_optimal[i], " and C", i, " = ", Ci_optimal[i], "."))
      } else {  ### multiple data points with nearly constant error

        D_tmp$error_constant <- "partially"

        if (verbose) print(paste0("Multiple points with optimal error."))

        if (all(R[ind_min_tol] == 0) | abs(diff(range(log2(R[ind_min_tol])))) < 1e-4) {  ## Ri is constant but Ci is not
          D_tmp$Ri <- bppg::geom_mean(R[ind_min_tol])
          D_tmp$Ri_min <- NA
          D_tmp$Ri_max <- NA
          D_tmp$case <- 4
          if (verbose) print(paste0("Constant Solution for R", i, ": ", bppg::geom_mean(R[ind_min_tol])))
        } else {  # Ri is not constant
          D_tmp$Ri <- NA
          D_tmp$Ri_min <- min(R[ind_min_tol])
          D_tmp$Ri_max <- max(R[ind_min_tol])
          D_tmp$case <- 5
          if (verbose) print(paste0("Range for R", i, ": ", BBmisc::collapse(range(R[ind_min_tol]), sep = " - "), "."))
        }
        if (verbose) print(paste0("Range for C", i, ": ", BBmisc::collapse(range(C[ind_min_tol]), sep = " - "), "."))
        D_tmp$Ci <- NA
        D_tmp$Ci_min <- min(C[ind_min_tol])
        D_tmp$Ci_max <- max(C[ind_min_tol])
      }
    }


    if (use_results_from_other_proteins) {
      ### see if solution can be enhanced by data form the other proteins
      ind <- which(abs(min(error) - X_tmp3$error) <= 1e-10)
      if (length(ind) > 0) {
        if (!is.na(D_tmp$Ri)) {
          if (abs(diff(range(log2(c(D_tmp$Ri, R_3[ind]))))) >= 1e-4) {  # if there was 1 solution before but there are multiple now

            D_tmp$Ri_min <- min(c(D_tmp$Ri, R_3[ind]))
            D_tmp$Ri_max <- max(c(D_tmp$Ri, R_3[ind]))
            D_tmp$Ri <- NA
          }
        } else {                                          # if the range of solutions can be enhanced
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














###################################################################################
###################################################################################
####### OLD version of the function
#' Extract protein ratio solutions from the result of iterate_over_Ci()
#'
#' @param S igraph object or list of biadjacency matrix (M) and fold changes (fc)
#' @param res result of iterate_over_Ci()
#' @param use_results_from_other_proteins if TRUE, use the results from other proteins within the same graph to calculate the optimal solution for each protein node
#' @param verbose if TRUE, print additional information
#' @param job job object from the BatchExperiment. Used to print the job id and parameters in the output
#' @param calc_op_sol if TRUE, calculate the optimal solution for each protein node without iterating over the Ci values
#' @param S_is_graph set to TRUE,if S is an igraph object and FALSE if S is a list with the biadjacency matrix and fold changes
#'
#' @return data frame
#'
#' @examples # TODO
automated_analysis_iterated_Ci_OLD <- function(S,
                                           res,
                                           use_results_from_other_proteins = FALSE,
                                           verbose = FALSE,
                                           job = NULL,
                                           calc_op_sol = TRUE,
                                           S_is_graph = FALSE) {

  # transform graph into incidence matrix is necessary
  if (S_is_graph & !is.null(S)) {
    X <- igraph::as_incidence_matrix(S)
    fc <- stats::na.omit(igraph::V(S)$pep_ratio)
    S <- list(X = X, fc = fc)
  }

  # calculate number of protein nodes
  if (!is.null(S)) {
    nr_proteins <- ncol(S$X)
  } else {
    nr_proteins <- 1
    res <- NULL
  }

  # calculate optimal solution without iterating the Ci
  if (calc_op_sol) {
    solution_optimal <- minimize_squared_error(S, #error.type = "multiplicative",
                                               fixed.Ci = NULL,
                                               verbose = FALSE, #error.trans = "square",
                                               reciprocal = FALSE, log_level = TRUE,
                                               control = list(trace = 0))

    error_optimal <- solution_optimal$RES$res_squ_err
    Ri_optimal <- solution_optimal$Ri
    Ci_optimal <- solution_optimal$Ci
  } else {
    error_optimal <- NA
    Ri_optimal <- rep(NA, nr_proteins)
    Ci_optimal <- rep(NA, nr_proteins)
  }



  if (nr_proteins == 1 & is.null(res)) {
    ### case 1: only one protein node. These cases were skipped during batch calculation,
    ### so know the optimal solution has to be calculated

    solution_optimal <- minimize_squared_error(S, error.type = "multiplicative", fixed.Ci = NULL,
                                                  verbose = FALSE, error.trans = "square",
                                                  reciprocal = FALSE, log_level = TRUE,
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

    if(!is.null(job)) {
      D_tmp$graphID <- job$pars$prob.pars$k
      D_tmp$comparison <- job$prob.name
      D_tmp$job.id <- job$job.id
    }

    if (verbose) print(paste0("Protein ", i))

    X_tmp <- res[res$protein == i,]  ### results only for the specific protein
    X_tmp <- X_tmp[!is.na(X_tmp$error),]
    C_tmp <- X_tmp[, paste0("C", i)]

    ### use also results from the other proteins in the same graph, but only if the Ci are not too extreme
    X_tmp3 <- res[res$protein != i,]
    X_tmp3 <- X_tmp3[!is.na(X_tmp3$error),]
    C_tmp3 <- X_tmp3[, paste0("C", i)]
    X_tmp3 <- X_tmp3[X_tmp3$protein == i |(C_tmp3 >= 0.01 & C_tmp3 < 0.99),]
    R_3 <- X_tmp3[, paste0("R", i)]
    C_3 <- X_tmp3[, paste0("C", i)]


    error <- X_tmp$error
    R <- X_tmp[, paste0("R", i)]
    C <- X_tmp[, paste0("C", i)]

    D_tmp$min_error <- min(error_optimal, min(error), na.rm = TRUE)
    D_tmp$error_optimal <- error_optimal


    #### 1. Check: Is the error constant over whole Ci in (0,1) (up to a tolerance)?
    if (abs(diff(range(X_tmp$error))) <= 1e-10) {
      ### error is constant, so either we have a single solution or a range

      D_tmp$error_constant <- "yes"

      if (verbose) print(paste0("Error is nearly constant (abs. diff. = ", abs(diff(range(error))), ")."))
      if (verbose) print(paste0("Mean error is ", mean(error, na.rm = TRUE), "."))
      if (verbose) print(paste0("Minimal error is ", min(error, na.rm = TRUE) , "."))
      ind_min <- which.min(error)

      #### 2. check: Is Ri constant (tolerance in log-scale)
      if (abs(diff(range(log2(R)))) > 1e-4) { ### Vergleich mit 1e-5??
        if (verbose) print(paste0("Range for R", i, ": ", BBmisc::collapse(range(R), sep = " - "), "."))
        D_tmp$Ri <- NA
        D_tmp$Ri_min <- min(R)
        D_tmp$Ri_max <- max(R)
        D_tmp$case <- 1
      } else {
        if (verbose) print(paste0("Constant Solution for R", i, ": ", bppg::geom_mean(R)))
        D_tmp$Ri <- bppg::geom_mean(R)
        D_tmp$Ri_min <- NA
        D_tmp$Ri_max <- NA
        D_tmp$case <- 2
      }

      if (verbose) print(paste0("Range for C", i, ": ", BBmisc::collapse(range(C), sep = " - "), "."))
      D_tmp$Ci_min <- min(C)
      D_tmp$Ci_max <- max(C)

    } else {
      ### single solution with minimal error or a range of Cis that is possible (error partly constant)

      error <- X_tmp$error
      R <- X_tmp[, paste0("R", i)]
      C <- X_tmp[, paste0("C", i)]

      ### check, if the optimization without any constraints really leads to the optimal result
      if (verbose & !is.na(error_optimal) & any(error < error_optimal)) {
        warning(paste0("Lower than optimal error detected. Difference = ", abs(min(error) - error_optimal)))
      }
      ### TODO: dont warn if error is extremely small (near 0 is!) or distance of the two errors is very very low

      ind_min <- which.min(error)
      ind_min_tol <- which(abs(min(error) - error) <= 1e-10)  ### the error is constant in which area?

      if (length(ind_min_tol) <= 1) {   ## at least one data point close to the minimal error term

        D_tmp$error_constant <- "no"
        D_tmp$Ri <- R[ind_min_tol]
        D_tmp$Ri_min <- NA
        D_tmp$Ri_max <- NA
        D_tmp$Ci <- C[ind_min_tol]
        D_tmp$Ci_min <- NA
        D_tmp$Ci_max <- NA
        D_tmp$case <- 3

        if (verbose) print(paste0("Single optimal solution with error ", min(error, na.rm = TRUE), "."))

        #### TODO: adjust message to not use Ri_optimal and Ci_optimal
        if (verbose) print(paste0("Single optimal solution is R", i, "= ", Ri_optimal[i], " and C", i, " = ", Ci_optimal[i], "."))
      } else {  ### multiple data points with constant error term

        D_tmp$error_constant <- "partially"

        if (verbose) print(paste0("Multiple points with optimal error."))

        if (abs(diff(range(log2(R[ind_min_tol])))) < 1e-4) {  ##  Ri is constant but Ci is not
          D_tmp$Ri <- bppg::geom_mean(R[ind_min_tol])
          D_tmp$Ri_min <- NA
          D_tmp$Ri_max <- NA
          D_tmp$case <- 4
          if (verbose) print(paste0("Constant Solution for R", i, ": ", bppg::geom_mean(R[ind_min_tol])))
        } else {  # Ri is not constant
          D_tmp$Ri <- NA
          D_tmp$Ri_min <- min(R[ind_min_tol])
          D_tmp$Ri_max <- max(R[ind_min_tol])
          D_tmp$case <- 5
          if (verbose) print(paste0("Range for R", i, ": ", BBmisc::collapse(range(R[ind_min_tol]), sep = " - "), "."))
        }
        if (verbose) print(paste0("Range for C", i, ": ", BBmisc::collapse(range(C[ind_min_tol]), sep = " - "), "."))
        D_tmp$Ci <- NA
        D_tmp$Ci_min <- min(C[ind_min_tol])
        D_tmp$Ci_max <- max(C[ind_min_tol])
      }
    }


    if (use_results_from_other_proteins) {
      ### see if solution can be extended by information on other proteins
      ind <- which(abs(min(error) - X_tmp3$error) <= 1e-10)
      if (length(ind) > 0) {
        if (!is.na(D_tmp$Ri)) {
          if (abs(diff(range(log2(c(D_tmp$Ri, R_3[ind]))))) >= 1e-4) {  # if there was 1 solution before but now there are multiple

            D_tmp$Ri_min <- min(c(D_tmp$Ri, R_3[ind]))
            D_tmp$Ri_max <- max(c(D_tmp$Ri, R_3[ind]))
            D_tmp$Ri <- NA
          }
        } else {                                          # if solution range can be extended
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







