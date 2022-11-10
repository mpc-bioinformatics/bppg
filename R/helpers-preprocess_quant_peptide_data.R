

### TODO: read in data directly from MaxQuant and filter unnecessary columns and decoys
### TODO: Normalization


#### data set with peptide sequence and intensities plus group information




# file <- "inst/extdata/peptides_D2_prep.txt"
#
#   #system.file("extdata", "peptides_D2_prep.txt", package = "bppg")
# D <- read.table(file, header = TRUE, sep = "\t")
# peptides <- D[,1]
# intensities <- D[,-1]
#
# intensities[intensities == 0] <- NA
#
# group <- factor(limma::strsplit2(colnames(intensities), "_")[,1],
#                 levels = c("UPS50amol", "UPS125amol", "UPS250amol", "UPS500amol", "UPS12500amol", "UPS2500amol", "UPS5000amol", "UPS25000amol", "UPS50000amol"))
#
# intensities <- cbind(pep_sequence = peptides, intensities)
#
# aggr_intensities <- aggregate_replicates(intensities, group)
#
# calculate_peptide_ratios(aggr_intensities)



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
aggregate_replicates <- function(intensities, group,  missing.limit = 0, method = "mean") {


  peptides <- intensities[,1]
  intensities <- intensities[,-1]


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


