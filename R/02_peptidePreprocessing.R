#' Functions in this file:
#' .foldChange
#' aggregateReplicates
#' calculatePeptideRatios


#' Calculate peptide ratios for pairwise comparisons of groups (Y/X).
#'
#' @param D       \strong{data.frame} \cr
#'                The data set.
#' @param X       \strong{character} \cr
#'                The column name of group1.
#' @param Y       \strong{character} \cr
#'                The column name of group2.
#' @param useNA   \strong{logical} \cr
#'                If \code{TRUE},results 0 and Inf are possible, otherwise
#'                ratio is NA if value for X or Y is NA
#'
#' @return The fold changes (Y/X).
#'
#'
#' @examples ## TODO
#'

.foldChange <- function(D, X, Y, useNA = FALSE) {
    FC <- D[, Y] / D[, X]

    if (useNA) {
        FC[is.na(D[, Y]) & !is.na(D[, X])] <- 0
        FC[is.na(D[, X]) & !is.na(D[, Y])] <- Inf
    }

    return(FC)
}


#' Aggregate replicates of the same experimental group.
#'
#' @param D               \strong{data.frame} \cr
#'                        The data set containing the peptide intensities.
#' @param missing.limit   \strong{numeric} \cr
#'                        The proportion of missing values that is allowed 
#'                        (e.g. 0 means no missings allowed).
#' @param method          \strong{character} \cr
#'                        The method of aggregation. Options are 
#'                        "mean", "sum" or "median"
#' @param group           \strong{character factor} \cr
#'                        The groups for aggregation.
#' @param id_cols         \strong{integer vector} \cr
#'                        The column numbers that contain peptide sequences etc
#'                        (everything except intensities).
#' @param imp_method      \strong{character} \cr
#'                        Chosen imputation optional approach, current method: 
#'                        "min_2_impute"
#'
#' @return A data set with aggregated intensities.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "peptides.txt", package = "bppg")
#' D <- readMqPeptideTable(path = file, LFQ = TRUE, remove_contaminants = FALSE)
#' group <- factor(rep(1:9, each = 3))
#' aggregateReplicates(D, group = group)

aggregateReplicates <- function(D, group, missing.limit = 0, method = "mean",
                                 id_cols = 1, imp_method = NULL) {

  id <- D[, id_cols, drop = FALSE]
  intensities <- D[, -(id_cols)]

  min_row <- apply(intensities, 1, min, na.rm = TRUE)

  res <- NULL
  mask_imputed <- NULL  # track imputed vales
  for (i in 1:length(levels(group))) {
    X_tmp <- intensities[, group == levels(group)[i]]

    FUN <- switch(method,
                  mean  = rowMeans,
                  sum = rowSums,
                  median = robustbase::rowMedians)

    X_tmp <- as.matrix(X_tmp)
    res_tmp <- FUN(X_tmp, na.rm = TRUE)

    missingx <- apply(X_tmp, 1, function(x) mean(is.na(x)))
    mask_tmp <- c(missingx > missing.limit | missingx == 1)
    res_tmp[mask_tmp] <- NA

    # apply imputation on missing values
    if (!is.null(imp_method)){
      FUN <- switch(imp_method,
                    min_2_imp = min_2_impute)

      vals_imp <- FUN(X_tmp, min_row)
      res_tmp[mask_tmp] <- vals_imp[mask_tmp]  # only replace missing values
    }

    res <- cbind(res, res_tmp)
    mask_imputed <- cbind(mask_imputed, mask_tmp)
  }

  # mask to remove rows with only imputed values, all under missing.limit
  all_imputed <- apply(mask_imputed, 1, all)

  res <- as.data.frame(res)
  colnames(res) <- levels(group)
  res <- data.frame(id, res)[!all_imputed, ]

  mask_imputed <- as.data.frame(mask_imputed)
  colnames(mask_imputed) <- levels(group)
  mask_imputed <- data.frame(id, mask_imputed)[!all_imputed, ]
  
  return(list(agg = res, imp = mask_imputed))
}

#' Calculation of peptide ratios from aggregated intensities.
#'
#' @param aggr_intensities   \strong{data.frame} \cr
#'                           The result from function [aggregateReplicates()].
#' @param id_cols            \strong{integer vector} \cr
#'                           The column numbers that contain peptide sequences
#'                           etc (everything except intensities).
#' @param group_levels       \strong{character factor} \cr
#'                           The levels of groups in the right order.
#' @param type               \strong{character} \cr
#'                           Options "ratio" or "difference". Use "difference"
#'                           if values are already on log-scale.
#' @param log_base           \strong{numeric} \cr
#'                           The log base.
#'
#' @return A data set with peptide ratios.
#' @export
#'
#' @examples ## TODO
#'

calculatePeptideRatios <-  function(data, id_cols = 1,
                                    group_levels = NULL, type = "ratio",
                                    log_base = 10) {

  aggr_intensities <- data$agg
  mask_impute <- data$imp

  id <- aggr_intensities[, id_cols, drop = FALSE]
  aggr_intensities <- aggr_intensities[, -(id_cols)]
  mask_impute <- mask_impute[, -(id_cols)]

  if (is.null(group_levels)) {
    group_levels <- factor(colnames(aggr_intensities), levels = colnames(aggr_intensities))
  }


  peptide_ratios <- list()
  k <- 1
  # peptide_ratios <- NULL
  # imputed_ratios <- NULL
  for (i in 1:(length(group_levels) - 1)) {
    for (j in (i + 1):length(group_levels)) {

      col1 <- which(colnames(aggr_intensities) == group_levels[i])
      col2 <- which(colnames(aggr_intensities) == group_levels[j])

      name <- paste0("ratio_", group_levels[i], "_", group_levels[j])

      if (type == "ratio") {
        FC <- foldChange(D = aggr_intensities, X = col1, Y = col2)
      }
      if (type == "difference") {
        FC <- aggr_intensities[, col2] - aggr_intensities[, col1]
        FC <- log_base^FC
      }

      dub_fc_mask <- apply(mask_impute[, c(col1, col2)], 1,
                          function(x) x[1] & x[2])
      FC[dub_fc_mask] <- NA   #remove ratio of two imputed values
      imp_fc_mask <- apply(mask_impute[, c(col1, col2)], 1,
                          function(x) x[1] | x[2])

      peptide_ratios[[k]] <- data.frame(id, FC, imputed = imp_fc_mask)
      colnames(peptide_ratios[[k]])[c(1, ncol(peptide_ratios[[k]]) - 1)] <- c("Sequence", name)
      names(peptide_ratios)[[k]] <- name
      k <- k + 1
    }
  }

  return(peptide_ratios)
}