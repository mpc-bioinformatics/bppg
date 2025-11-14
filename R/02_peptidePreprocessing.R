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
#' @examples 
#' D <- data.frame(s1 = c(1,4), s2 = c(2,5), s3 = c(3,6))
#' X <- "s1"
#' Y <- "s2"
#' bppg:::.foldChange(D,X,Y)

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
    id_cols = 1) {

    id <- D[, id_cols, drop = FALSE]
    intensities <- D[, -(id_cols)]

    FUN <- switch(method,
        mean  = rowMeans,
        sum = rowSums,
        median = robustbase::rowMedians)

    res <- NULL
    # TODO kann das in ein apply verwandelt werden? split(intensities, group)?
    for (i in 1:length(levels(group))) {

        X_tmp <- intensities[, group == levels(group)[i]]
        X_tmp <- as.matrix(X_tmp)

        res_tmp <- FUN(X_tmp, na.rm = TRUE)
        missingx <- apply(X_tmp, 1, function(x) mean(is.na(x)))
        res_tmp[missingx > missing.limit | missingx == 1] <- NA
        res <- cbind(res, res_tmp)
    }

    res <- as.data.frame(res)
    colnames(res) <- levels(group)
    res <- data.frame(id, res)
    return(res)
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

calculatePeptideRatios <- function(aggr_intensities, id_cols = 1,
    group_levels = NULL, type = "ratio",
    log_base = 10) {

    id <- aggr_intensities[, id_cols, drop = FALSE]
    aggr_intensities <- aggr_intensities[, -(id_cols)]

    if (is.null(group_levels)) {
        group_levels <- factor(colnames(aggr_intensities),
            levels = colnames(aggr_intensities))
    }

    peptide_ratios <- NULL ## could we allocate that

    ## TODO VAPPLY?
    for (i in 1:(length(group_levels) - 1)) {
        for (j in (i + 1):length(group_levels)) {
            col1 <- which(colnames(aggr_intensities) == group_levels[i])
            col2 <- which(colnames(aggr_intensities) == group_levels[j])

            name <- paste0("ratio_", group_levels[i], "_", group_levels[j])

            if (type == "ratio") {
                FC <- .foldChange(D = aggr_intensities, X = col1, Y = col2)
            }
            if (type == "difference") {
                FC <- aggr_intensities[, col2] - aggr_intensities[, col1]
                FC <- log_base^FC
            }

            peptide_ratios <- cbind(peptide_ratios, FC)
            colnames(peptide_ratios)[ncol(peptide_ratios)] <- name
        }
    }

    return(data.frame(id, peptide_ratios))
}
