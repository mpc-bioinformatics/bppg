# Functions in this file:
# .foldChange
# aggregateReplicates
# calculatePeptideRatios
# normalizePeptideIntensities


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

aggregateReplicates <- function(D, 
    group = NULL, 
    missing.limit = 0, 
    method = "mean",
    seq_col = "Sequence", 
    imp_method = NULL) {
    checkmate::assertClass(D, "SummarizedExperiment")
    checkmate::assertDataFrame(SummarizedExperiment::assays(D)$intensities,
        all.missing=FALSE)
    checkmate::assertFactor(group, null.ok = TRUE)
    checkmate::assertNumber(missing.limit, lower = 0, upper = 1)
    checkmate::assertCharacter(method, pattern = "mean|sum|median")
    checkmate::assertCharacter(seq_col)

    id <- SummarizedExperiment::rowData(D)[, seq_col]
    intensities <- SummarizedExperiment::assays(D)$intensities

    if(is.null(group)){
        group <- factor(SummarizedExperiment::colData(D)$group)
    }
    min_row <- apply(intensities, 1, min, na.rm = TRUE)
    mask_imputed <- NULL  # track imputed vales

    FUN <- switch(method,
        mean  = rowMeans,
        sum = rowSums,
        median = robustbase::rowMedians)

    res <- vapply(1:length(levels(group)), function(i) {
        X_tmp <- intensities[, group == levels(group)[i]]
        X_tmp <- as.matrix(X_tmp)

        res_tmp <- FUN(X_tmp, na.rm = TRUE)
        missingx <- rowMeans(is.na(X_tmp))
        mask_tmp <- c(missingx > missing.limit | missingx == 1)
        res_tmp[mask_tmp] <- NA

        # apply imputation on missing values
        if (!is.null(imp_method)){
            FUN <- switch(imp_method,
                            min_2_imp = min_2_impute)

            vals_imp <- FUN(X_tmp, min_row)
            res_tmp[mask_tmp] <- vals_imp[mask_tmp]  # only replace missing values
        }
        mask_imputed <- cbind(mask_imputed, mask_tmp) ### how can I solve this differently?
        return(res_tmp)
    }, numeric(length(id)))

    res <- as.data.frame(res)
    colnames(res) <- levels(group)
    rownames(res) <- id

    mask_imputed <- as.data.frame(mask_imputed)
    colnames(mask_imputed) <- levels(group)
    all_imputed <- apply(mask_imputed, 1, all) # remove rows with only imputed, would that not be the on OFF case?
    mask_imputed <- data.frame(id, mask_imputed)[!all_imputed, ]

    res <- SummarizedExperiment::SummarizedExperiment(
        assays = list(intensities = res, maskImputation = mask_imputed),
        colData = data.frame(group = colnames(res)),
        rowData = SummarizedExperiment::rowData(D))
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
    checkmate::assertVector(group_levels, unique = TRUE, null.ok = TRUE)

    aggr_intensities <- SummarizedExperiment::assays(D)$intensities
    mask_impute <- SummarizedExperiment::assay(D)$maskImputation  # muss hier was mit der id verändert werden?

    if (is.null(group_levels)) {
        group_levels <- SummarizedExperiment::colData(D)$group
    }

    # create pairwise groups for ratio calculation
    groupCombinations <- utils::combn(group_levels, 2)
    peptide_log_ratios <- vapply(seq_len(ncol(groupCombinations)), function(i) {
        log2(.foldChange(D = aggr_intensities, X = groupCombinations[1, i],
            Y = groupCombinations[2, i]))
    }, numeric(nrow(aggr_intensities)))



    peptide_log_ratios <- data.frame(peptide_log_ratios)
    colnames(peptide_log_ratios) <- paste0("logRatio_", groupCombinations[1, ], "_",
        groupCombinations[2, ])
    rownames(peptide_log_ratios) <- rownames(aggr_intensities)

    dub_fc_mask <- apply(mask_impute[, c(col1, col2)], 1,
                          function(x) x[1] & x[2])
    peptide_log_ratios[dub_fc_mask] <- NA   #remove ratio of two imputed values
    imp_fc_mask <- apply(mask_impute[, c(col1, col2)], 1,
                          function(x) x[1] | x[2])


    res <- SummarizedExperiment::SummarizedExperiment(
        assays = list(logRatios = peptide_log_ratios, maskImputation = imp_fc_mask), 
        colData = data.frame(comparison = colnames(peptide_log_ratios)),
        rowData = SummarizedExperiment::rowData(D))
    return(res)
}







#' Normalization of peptide intensities
#'
#' @param D \strong{SummarizedExperiment} \cr
#'          Dataset containing peptide intensities, e.g. the result of
#'          [readMqPeptideTable].
#' @param method \strong{character} \cr
#'          The method of normalization. Options are "nonorm"
#'          (no normalization), "median", "loess",  "quantile" or "lts"
#'          normalization. Default is "loess"
#' @param lts.quantile \strong{numeric} \cr
#'          The quantile for the lts normalization. Default is 0.8.
#' @param log_base \strong{numeric} \cr
#'          The base for log-transformation. Default is 2.
#' @returns \strong{SummarizedExperiment} \cr
#'          Dataset containing normalized peptide intensities.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "peptides.txt", package = "bppg")
#' D <- readMqPeptideTable(path = file, LFQ = TRUE, remove_contaminants = FALSE)
#' D_norm <- normalizePeptideIntensities(D, method = "loess")
normalizePeptideIntensities <- function(D, method = "loess", lts.quantile = 0.8,
                                        log_base = 2) {

    DATA <- SummarizedExperiment::assays(D)$intensities


    if (method %in% c("loess", "quantile", "median")) {

        log_DATA <- log(DATA, base = log_base)
        #### choose normalization function
        fun <- limma::normalizeBetweenArrays
        args <- switch(method,
                       "loess" = list(object = log_DATA, method = "cyclicloess"),
                       "quantile" = list(object = log_DATA, method = "quantile"),
                       "median" = list(object = log_DATA, method = "scale"))

        DATA_norm <- do.call(fun, args)
        DATA_norm <- as.data.frame(DATA_norm)
        DATA_norm <- log_base^DATA_norm # re-transform
    }

    if (method == "lts") {
        # Does not need log-transformation, as it does a glog trans
        # (similar to log2)
        DATA_norm <- vsn::vsn2(as.matrix(DATA), lts.quantile = lts.quantile)
        DATA_norm <- DATA_norm@hx
        DATA_norm <- as.data.frame(DATA_norm)
        DATA_norm <- log_base^DATA_norm # re-transform
    }

    if (method == "nonorm") {
        DATA_norm <- DATA
    }
    res <- SummarizedExperiment::SummarizedExperiment(
        assays = DATA_norm,
        colData = SummarizedExperiment::colData(D),
        rowData = SummarizedExperiment::rowData(D))
    return(res)
}

