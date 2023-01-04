

### TODO: read in data directly from MaxQuant and filter unnecessary columns and decoys
### TODO: Normalization





#' Import of MaxQuant's peptide.txt-table
#'
#' @param path Path to the peptides.txt table
#' @param LFQ If TRUE, LFQ intensities are used, if FALSE, raw (unnormalized) intensities
#' @param remove_contaminants If TRUE, peptide sequences from potential contaminants are removed
#' @param rename_columns Rename columns? If TRUE, "Intensity." or "LFQ.intensity." are removed
#' @param zeroToNA If TRUE, zeros are converted to NAs.
#'
#' @return Dataframe with sequences and intensities
#' @export
#'
#' @examples
#' file <- system.file("extdata", "peptides.txt", package = "bppg")
#' D <- read_MQ_peptidetable(path = file, LFQ = TRUE, remove_contaminants = FALSE)
read_MQ_peptidetable <- function(path, LFQ = FALSE, remove_contaminants = FALSE,
                                 rename_columns = TRUE, zeroToNA = TRUE) {

  D <- utils::read.table(path, sep = "\t", header = TRUE)

  ### remove decoy entries:
  ind_decoy <- D$Reverse == "+"
  D <- D[!ind_decoy,]
  print(paste0("Removed ", sum(ind_decoy), " decoy sequences."))

  ind_cont <- D$Potential.contaminant == "+"
  if(remove_contaminants) {
    D <- D[!ind_cont,]
    print(paste0("Removed ", sum(ind_cont), " contaminant sequences."))
  }


  ## search for itensity columns or LFQ values
  if(LFQ) {
    intensities <- D[, grep("LFQ", colnames(D))]
    if (rename_columns) colnames(intensities) <- stringr::str_replace(colnames(intensities), "LFQ.intensity.", "")
  } else {
    intensities <- D[, grep("Intensity.", colnames(D))]
    if (rename_columns) colnames(intensities) <- stringr::str_replace(colnames(intensities), "Intensity.", "")
  }

  if(zeroToNA) {
    intensities[intensities == 0] <- NA
  }


  RES <- data.frame(Sequence = D$Sequence, intensities)

  return(RES)
}







#' aggregates replicates of the same experimental group
#'
#' @param D data set with peptide intensities
#' @param missing.limit proportion of missing values that is allowed (e.g. 0 means no missings allowed)
#' @param method "mean", "sum" oder "median"
#' @param group groups for aggregation as a factor
#' @param id_cols column numbers that contain peptide sequences etc (everything except intensities)
#'
#' @return data set with aggregated intensities
#' @export
#'
#' @examples
#' file <- system.file("extdata", "peptides.txt", package = "bppg")
#' D <- read_MQ_peptidetable(path = file, LFQ = TRUE, remove_contaminants = FALSE)
#' group <- factor(rep(1:9, each = 3))
#' aggregate_replicates(D, group = group)
aggregate_replicates <- function(D, group, missing.limit = 0, method = "mean",
                                 id_cols = 1) {


  peptides <- D$Sequence
  intensities <- D[,-(id_cols)]

  ### TODO: Aus Spaltennamen Gruppe selber erschließen?

  res <- NULL
  for (i in 1:length(levels(group))) {

    X_tmp <- intensities[, group == levels(group)[i]]

    FUN <- switch(method,
                  mean  = rowMeans,
                  sum = rowSums,
                  median = robustbase::rowMedians)

    X_tmp <- as.matrix(X_tmp)

    res_tmp <- FUN(X_tmp, na.rm = TRUE)

   # if (!use0) {
      missingx <- apply(X_tmp, 1, function(x) mean(is.na(x)))
      res_tmp[missingx > missing.limit | missingx == 1] <- NA
  #  }

    res <- cbind(res, res_tmp)
  }

  res <- as.data.frame(res)
  colnames(res) <- levels(group)
  res <- cbind(D[, id_cols], res)
  return(res)
}



#' calculates peptide ratios for pairwise comparisons of groups (Y/X)
#'
#' @param D data set
#' @param X group1 (column name)
#' @param Y group2 (column name)
#' @param useNA if TRUE, results 0 and Inf are possible, otherwise ratio is NA if value for X or Y is NA
#'
#' @return fold changes (Y/X)
#' @export
#'
#' @examples
#' ### TODO
foldChange <- function(D, X, Y, useNA = FALSE) {
  FC <- D[, Y] / D[, X]

  if(useNA) {
    FC[is.na(D[,Y]) & !is.na(D[,X])] <- 0
    FC[is.na(D[,X]) & !is.na(D[,Y])] <- Inf
  }

  return(FC)
}



#' Calculation of peptide ratios from aggregated intensities
#'
#' @param aggr_intensities result from function aggregate_replicates
#' @param id_cols column numbers that contain peptide sequences etc (everything except intensities)
#' @param group_levels levels of groups in the right order
#'
#' @return data set with peptide ratios
#' @export
#'
#' @examples
#' ## TODO
calculate_peptide_ratios <- function(aggr_intensities, id_cols = 1, group_levels = NULL) {

  id <- aggr_intensities[,id_cols]
  aggr_intensities <- aggr_intensities[,-(id_cols)]

  if(is.null(group_levels)) {
    group_levels <- factor(colnames(aggr_intensities), levels = colnames(aggr_intensities))
  }


  peptide_ratios <- NULL

  for (i in 1:(length(group_levels)-1)) {
    for (j in (i+1):length(group_levels)) {

      col1 <- which(colnames(aggr_intensities) == group_levels[i])
      col2 <- which(colnames(aggr_intensities) == group_levels[j])

      name <- paste0("ratio_", group_levels[j], "_", group_levels[i])
      FC <- foldChange(D = aggr_intensities, X = col1, Y = col2)
      peptide_ratios <- cbind(peptide_ratios, FC)
      colnames(peptide_ratios)[ncol(peptide_ratios)] <- name
    }
  }

  return(cbind(id, peptide_ratios))
}


