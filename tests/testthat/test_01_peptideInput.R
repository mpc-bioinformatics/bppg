test_that("read FragPipe Output table", {
    file <- system.file("extdata", "peptides_FP_filtered.tsv", package = "bppg")

    D1 <- bppg::readFpPeptideTable(file)
    D2 <- bppg::readFpPeptideTable(file,
        further_columns_to_keep = c("Protein", "Mapped.Proteins"))

    expect_snapshot(D1)
    expect_snapshot(ceiling(SummarizedExperiment::assays(D1)$intensities))
    expect_snapshot(as.data.frame(tail(
        SummarizedExperiment::rowData(D1), n = 1000)))
    expect_snapshot(as.data.frame(tail(
        SummarizedExperiment::colData(D1), n = 1000)))

    expect_snapshot(D2)
    expect_snapshot(ceiling(SummarizedExperiment::assays(D2)$intensities))
    expect_snapshot(as.data.frame(tail(
        SummarizedExperiment::rowData(D2), n = 1000)))
    expect_snapshot(as.data.frame(tail(
        SummarizedExperiment::colData(D2), n = 1000)))
})


test_that("read MaxQuant Output table", {
    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")

    D1 <- bppg::readMqPeptideTable(file)
    D2 <- bppg::readMqPeptideTable(file,
                                   LFQ = TRUE,
                                   further_columns_to_keep = c("Proteins", "Score"))
    expect_snapshot(D1)
    expect_snapshot(SummarizedExperiment::assays(D1)$intensities)
    expect_snapshot(as.data.frame(tail(
        SummarizedExperiment::rowData(D1), n = 1000)))
    expect_snapshot(as.data.frame(tail(
        SummarizedExperiment::colData(D1), n = 1000)))

    expect_snapshot(D2)
    expect_snapshot(SummarizedExperiment::assays(D2)$intensities)
    expect_snapshot(as.data.frame(tail(
        SummarizedExperiment::rowData(D2), n = 1000)))
    expect_snapshot(as.data.frame(tail(
        SummarizedExperiment::colData(D2), n = 1000)))
})
