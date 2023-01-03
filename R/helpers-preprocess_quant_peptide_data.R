

### TODO: read in data directly from MaxQuant and filter unnecessary columns and decoys
### TODO: Normalization





#' Import of MaxQuant's peptide.txt-table
#'
#' @param path Path to the peptides.txt table
#' @param LFQ If TRUE, LFQ intensities are used, if FALSE, raw (unnormalized) intensities
#' @param remove_contaminants If TRUE, peptide sequences from potential contaminants are removed
#'
#' @return Dataframe with sequences and intensities
#' @export
#'
#' @examples
#' file <- system.file("extdata", "peptides.txt", package = "bppg")
#' D <- read_MQ_peptidetable(path = file, LFQ = TRUE, remove_contaminants = FALSE)
read_MQ_peptidetable <- function(path, LFQ = FALSE, remove_contaminants = FALSE,
                                 rename_columns = TRUE, zeroToNA = TRUE) {

  D <- read.table(path, sep = "\t", header = TRUE)

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
#' @param intensities data set with peptide intensities
#' @param missing.limit proportion of missing values that is allowed (e.g. 0 means no missings allowed)
#' @param method "mean", "sum" oder "median"
#' @param group groups for aggregation as a factor
#'
#' @return data set with aggregated intensities
#' @export
#'
#' @examples
#' file <- system.file("extdata", "peptides.txt", package = "bppg")
#' D <- read_MQ_peptidetable(path = file, LFQ = TRUE, remove_contaminants = FALSE)
#' group <- factor(rep(1:9, each = 3))
#' aggregate_replicates(D)
aggregate_replicates <- function(D, group,  missing.limit = 0, method = "mean") {


  peptides <- D$Sequence
  intensities <- D[,-1]  ## TODO: es koennten weitere Spalten vorhanden sein!

  ### TODO: Aus Spaltennamen Gruppe selber erschließen!

  res <- NULL
  for (i in 1:length(levels(group))) {

    X_tmp <- D[, group == levels(group)[i]]

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
  res <- cbind(pep_sequence = peptides, res)
  return(res)
}



#### calculates peptide ratios for pairwise comparisons of groups
## D: dataset
## X: group1 (column name)
## Y: group2 (column name)
## useNA: if TRUE, results 0 and Inf are possible, otherwise ratio is NA if value for X or Y is NA
### result: fold changes (Y/X)
foldChange <- function(D, X, Y, useNA = FALSE) {
  FC <- D[, Y] / D[, X]

  if(useNA) {
    FC[is.na(D[,Y]) & !is.na(D[,X])] <- 0
    FC[is.na(D[,X]) & !is.na(D[,Y])] <- Inf
  }

  return(FC)
}



calculate_peptide_ratios <- function(aggr_intensities) {

  peptides <- aggr_intensities[,1]
  aggr_intensities <- aggr_intensities[,-1]

  peptide_ratios <- NULL
  for (i in 1:(ncol(aggr_intensities)-1)) {
    for (j in 2:ncol(aggr_intensities)) {
      name <- paste0("ratio_", colnames(aggr_intensities)[j], "_", colnames(aggr_intensities)[i])
      FC <- foldChange(D = aggr_intensities, X = colnames(aggr_intensities)[i], Y = colnames(aggr_intensities)[j])
      peptide_ratios <- cbind(peptide_ratios, FC)
      colnames(peptide_ratios)[ncol(peptide_ratios)] <- name
    }
  }

  return(cbind(pep_sequence = peptides, peptide_ratios))
}


