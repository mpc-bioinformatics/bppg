
test_that("read MaxQuant Output table", {
    D1 <- bppg::readMqPeptideTable(test_path("testfiles/peptides.txt"))
    D2 <- bppg::readMqPeptideTable(test_path("testfiles/peptides.txt"),
                                   LFQ = TRUE,
                                   further_columns_to_keep = c("Proteins", "Score"))
    expect_snapshot(D1)
    expect_snapshot(D2)
})