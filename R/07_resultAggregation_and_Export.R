

#' Combine relative quantification results over multiple pairwise comparisons
#' into a single SummarizedExperiments object.
#'
#' @param compResultList \strong{list} \cr
#'      List of results to be concatenated. Each list element is the result of
#'      the optimization step for one pairwise comparison, i.e. a
#'      SummarizedExperiment object. Ideally, the element names should
#'      correspond to the comparison names.
#'
#' @returns A SummarizedExperiment object with one assay per comparison. The row
#' names of the assay matrices are the union of all row names across the input
#' SummarizedExperiment objects. For rows that are not present in a given
#' comparison, NA values are introduced. Row names are sorted alphabetically.
#' @export
#'
#' @examples
#' # prepared result list. Check "inst/scripts/filter_fasta_and_quantData.R" for
#' # details on how it was generated.
#' file <- system.file("extdata", "resultsList.rds", package = "bppg")
#' resultList <- readRDS(file)
#' SE <- combineComparisons(resultList)
#'
combineComparisons <- function(compResultList) {

    # get union of all row names and sort them
    all_rows <- Reduce(union, lapply(compResultList, function(x) {
        return(rownames(SummarizedExperiment::assay(x)))
    }), init = NULL)

    # helper function to harmonize row names
    # (introduce NA rows for missing proteins)
    harmonizeRowNames <- function(se) {
        df <- SummarizedExperiment::assays(se)$results

        # introduce NA rows for all protein groups not present in this
        # comparison
        missing <- setdiff(all_rows, rownames(df))
        na_rows <- matrix(NA, nrow = length(missing), ncol = ncol(df),
            dimnames = list(missing, colnames(df)))
        df_tmp <- rbind(df, na_rows)

        df_tmp <- df_tmp[order(rownames(df_tmp)), ] # sort by rowname
        return(df_tmp)
    }

    df_list_harmonized <- lapply(compResultList, harmonizeRowNames)
    names(df_list_harmonized) <- names(compResultList)

    SE <- SummarizedExperiment::SummarizedExperiment(
        assays = df_list_harmonized,
        rowData = data.frame(accession = rownames(df_list_harmonized[[1]])),
        colData = data.frame(colnames = colnames(df_list_harmonized[[1]])))

    return(SE)
}






#' Export SummarizedExperiment object to an Excel file,
#' with one assay per sheet.
#'
#' @param SE \strong{SummarizedExperiment object} \cr Object to be exported.
#'      E.g. the result of \code{\link{combineComparisons}}.
#' @param file \strong{character} \cr file path to the output Excel file.
#'      If the file already exists, it will be overwritten.
#'
#' @returns nothing, but an Excel file is written containing the assay matrices.
#' Each assay is written to a separate sheet, and the sheet name corresponds to
#' the assay name. Row names of the assay matrices are included
#' in the Excel file.
#' @export
#'
#' @examples
#' # prepared result list. Check "inst/scripts/filter_fasta_and_quantData.R" for
#' # details on how it was generated.
#' file <- system.file("extdata", "resultsList.rds", package = "bppg")
#' SE <- combineComparisons(file)
#' exportSE(SE, "results.xlsx")
#'
#' @importFrom openxlsx addWorksheet createWorkbook saveWorkbook writeData
#' @importFrom SummarizedExperiment assay assayNames
exportSE <- function(SE, file) {

    wb <- openxlsx::createWorkbook()

    for (assay_name in SummarizedExperiment::assayNames(SE)) {
        openxlsx::addWorksheet(wb, assay_name)
        mat <- as.data.frame(SummarizedExperiment::assay(SE, assay_name))
        openxlsx::writeData(wb, sheet = assay_name, x = mat, rowNames = TRUE,
            keepNA = TRUE)
    }
    ## TODO: cbind rowData if available?

    openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
    return(invisible(NULL))
}








