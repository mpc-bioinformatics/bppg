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
#'                If \code{TRUE}, results 0 and Inf are possible, otherwise
#'                ratio is NA, if value for X or Y is NA
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
#' @param D              \strong{SummarizedExperiment} \cr
#'                       The data experiment containing the peptide intensities.
#' @param group          \strong{character factor} \cr
#'                       The groups for aggregation, if not already in 
#'                       SummarizedExperiment::colData(D)$group. 
#' @param missing.limit  \strong{numeric} \cr
#'                       The proportion of missing values that is allowed 
#'                       (e.g. 0 means no missings allowed).
#' @param method         \strong{character} \cr
#'                       The method of aggregation. Options are 
#'                       "mean", "sum" or "median"
#' @param seq_col         \strong{character} \cr
#'                       The column name containaining the peptide sequences 
#'                       in the rowData of the SummarizedExperiment. 
#'                       Default is "Sequence"
#' @param imp_method      \strong{character} \cr
#'                        Chosen imputation optional approach, current method: 
#'                        "min_2_impute"
#'
#' @return A SummarizedExperiment with aggregated intensities ($intensities).
#' @export
#'
#' @examples
#' file <- system.file("extdata", "peptides.txt", package = "bppg")
#' D <- readMqPeptideTable(path = file, LFQ = TRUE, remove_contaminants = FALSE)
#' group <- factor(rep(1:9, each = 3))
#' aggregateReplicates(D, group = group)

aggregateReplicates <- function(D, group = NULL, missing.limit = 0, 
    method = "mean", seq_col = "Sequence", imp_method = NULL) {
    checkmate::assertClass(D, "SummarizedExperiment") 
    checkmate::assertDataFrame(SummarizedExperiment::assays(D)$intensities, 
        all.missing=FALSE)
    checkmate::assertFactor(group, null.ok = TRUE)
    checkmate::assertNumber(missing.limit, lower = 0, upper = 1)
    checkmate::assertCharacter(method, pattern = "mean|sum|median")
    checkmate::assertCharacter(seq_col)
    checkmate::assertCharacter(imp_method, pattern = "min_2_impute", null.ok = TRUE)

    id <- SummarizedExperiment::rowData(D)[, seq_col]
    intensities <- SummarizedExperiment::assays(D)$intensities

    if (is.null(group)) {
        group <- factor(SummarizedExperiment::colData(D)$group)
    }
    min_row <- apply(intensities, 1, min, na.rm = TRUE) # only needed for min2impute, there is rowMeans also TODO
    mask_impute <- NULL  # track imputed vales, missleading name

    FUN <- switch(method,
        mean  = rowMeans,
        sum = rowSums,
        median = robustbase::rowMedians)

    # Track missingness and later imputed values
    mask_impute <- vapply(1:length(levels(group)), function(i) {
        X_tmp <- intensities[, group == levels(group)[i]]
        X_tmp <- as.matrix(X_tmp)

        missingx <- rowMeans(is.na(X_tmp))
        mask_tmp <- c(missingx > missing.limit | missingx == 1)
        return(mask_tmp)
    }, FUN.VALUE = logical(length(id)))

    # doofe kombi mit impute hier drin,, vielelicht lieber raus bewegen?
    # also das ganze vapply?
    res <- vapply(1:length(levels(group)), function(i, mask_impute) {
        X_tmp <- intensities[, group == levels(group)[i]]
        X_tmp <- as.matrix(X_tmp)

        res_tmp <- FUN(X_tmp, na.rm = TRUE)
        res_tmp[mask_impute[, i]] <- NA

        # apply imputation on missing values
        if (!is.null(imp_method)) {
            FUN <- switch(imp_method,
                          min_2_imp = min_2_impute)
            vals_imp <- FUN(X_tmp, min_row)
            # only replace missing values
            res_tmp[mask_impute[, i]] <- vals_imp[mask_impute[, i]]
        }
        return(res_tmp)
    }, FUN.VALUE = numeric(length(id)), mask_impute)

    res <- as.data.frame(res)
    colnames(res) <- levels(group)
    rownames(res) <- id

    if (!is.null(imp_method)) {
        mask_impute <- as.data.frame(mask_impute)
        colnames(mask_impute) <- levels(group)
        rownames(mask_impute) <- id
        all_imputed <- apply(mask_impute, 1, all)
        res <- SummarizedExperiment::SummarizedExperiment(
            assays = list(intensities = res, maskImputation = mask_impute),
            colData = data.frame(group = colnames(res)),
            rowData = SummarizedExperiment::rowData(D),
            metadata = list(imputed = TRUE))
        res <- res[!all_imputed, ]

    } else {
        res <- SummarizedExperiment::SummarizedExperiment(
            assays = list(intensities = res),
            colData = data.frame(group = colnames(res)),
            rowData = SummarizedExperiment::rowData(D),
            metadata = list(imputed = FALSE))
    }
    return(res)
}

#' Calculation of peptide ratios from aggregated intensities.
#'
#' @param D                  \strong{SummarizedExperiment} \cr
#'                           The result from function [aggregateReplicates()].
#' @param group_levels       \strong{character factor} \cr
#'                           The levels of groups in the right order.
#'
#' @return A SummarizedExperiment with log2 peptide ratios (logRatios).
#' @export
#'
#' @examples
#' file <- system.file("extdata", "peptides.txt", package = "bppg")
#' D <- readMqPeptideTable(path = file, LFQ = TRUE, remove_contaminants = FALSE)
#' group <- factor(rep(1:9, each = 3))
#' dAgg <- aggregateReplicates(D, group = group)
#' calculatePeptideRatios(dAgg)

calculatePeptideRatios <- function(D, group_levels = NULL) {
    checkmate::assertClass(D, "SummarizedExperiment")
    checkmate::assertDataFrame(SummarizedExperiment::assays(
        D)$intensities, all.missing=FALSE)
    checkmate::assert_logical(S4Vectors::metadata(D)$imputed)
    checkmate::assertVector(group_levels, unique = TRUE, null.ok = TRUE)

    aggr_intensities <- SummarizedExperiment::assays(D)$intensities
    # aber nur falls das nicht null ist 

    if (is.null(group_levels)) {
        group_levels <- SummarizedExperiment::colData(D)$group
    }

    # create pairwise groups for ratio calculation
    groupCombinations <- combn(group_levels, 2)
    peptide_log_ratios <- vapply(seq_len(ncol(groupCombinations)), function(i) {
        log2(.foldChange(D = aggr_intensities, X = groupCombinations[1, i],
             Y = groupCombinations[2, i]))
    }, numeric(nrow(aggr_intensities)))

    peptide_log_ratios <- data.frame(peptide_log_ratios)
    colnames(peptide_log_ratios) <- paste0("logRatio_", 
        groupCombinations[1, ], "_", groupCombinations[2, ])
    rownames(peptide_log_ratios) <- rownames(aggr_intensities)

    if (S4Vectors::metadata(D)$imputed) {
        mask_impute <- SummarizedExperiment::assays(D)$maskImputation

        fakeFCMask <- vapply(seq_len(ncol(groupCombinations)), function(i) {
            mask_impute[, groupCombinations[1, i]] & 
                mask_impute[, groupCombinations[2, i]] 
        }, logical(nrow(aggr_intensities)))
        peptide_log_ratios[fakeFCMask] <- NA # remove ratio of imputed values
        imputedFCs <- vapply(seq_len(ncol(groupCombinations)), function(i) {
            mask_impute[, groupCombinations[1, i]] | 
                mask_impute[, groupCombinations[2, i]] 
        }, logical(nrow(aggr_intensities)))

        res <- SummarizedExperiment::SummarizedExperiment(
            assays = list(logRatios = peptide_log_ratios,
                maskImputation = imputedFCs), 
            colData = data.frame(comparison = colnames(peptide_log_ratios)),
            rowData = SummarizedExperiment::rowData(D),
            metadata = list(imputed = TRUE))
    } else {
        res <- SummarizedExperiment::SummarizedExperiment(
            assays = list(logRatios = peptide_log_ratios), 
            colData = data.frame(comparison = colnames(peptide_log_ratios)),
            rowData = SummarizedExperiment::rowData(D),
            metadata = FALSE)
    }
    return(res)
}
