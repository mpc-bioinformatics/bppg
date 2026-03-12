
test_that("read MaxQuant Output table", {
    D1 <- bppg::readMqPeptideTable(test_path("testfiles/peptides.txt"))
    D2 <- bppg::readMqPeptideTable(test_path("testfiles/peptides.txt"),
                                   LFQ = TRUE,
                                   further_columns_to_keep = c("Proteins", "Score"))
    expect_snapshot(SummarizedExperiment::assays(D1)$intensities)
    expect_snapshot(D1)
    expect_snapshot(SummarizedExperiment::assays(D2)$intensities)
    expect_snapshot(D2)
})


test_that("normalize peptide data", {
    D <- bppg::readMqPeptideTable(test_path("testfiles/peptides.txt"), LFQ = FALSE)
    D_norm_loess <- bppg::normalizePeptideIntensities(D, method = "loess")
    D_norm_lts <- bppg::normalizePeptideIntensities(D, method = "lts")
    expect_snapshot(SummarizedExperiment::assays(D_norm_loess)$intensities)
    expect_snapshot(SummarizedExperiment::assays(D_norm_lts)$intensities)
})

