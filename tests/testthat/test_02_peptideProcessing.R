test_that("test aggregateReplicates", {

    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    D <- bppg::readMqPeptideTable(file)
    D_norm <- bppg::normalizePeptideIntensities(D, method = "nonorm")

    D1 <- bppg::aggregateReplicates(D = D_norm, group = factor(rep(1:9, each = 3)))

    D2 <- bppg::aggregateReplicates(D = D_norm,
                                    group = factor(rep(1:9, each = 3)),
                                    missing.limit = 0.35,
                                    method = "median")

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



test_that("test calculatePeptideRatios", {

    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    D <- bppg::readMqPeptideTable(file)
    D_norm <- bppg::normalizePeptideIntensities(D, method = "nonorm")
    D1 <- bppg::aggregateReplicates(D = D_norm, group = factor(rep(1:9, each = 3)))
    D1 <- bppg::calculatePeptideRatios(D = D1)

    expect_snapshot(D1)
    expect_snapshot(SummarizedExperiment::assays(D1)$logRatios)
    expect_snapshot(SummarizedExperiment::rowData(D1))
    expect_snapshot(SummarizedExperiment::colData(D1))
})


test_that("normalize peptide data", {
    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    D <- bppg::readMqPeptideTable(file)
    D_norm_loess <- bppg::normalizePeptideIntensities(D, method = "loess")
    D_norm_lts <- bppg::normalizePeptideIntensities(D, method = "lts")

    expect_snapshot(D_norm_loess)
    expect_snapshot(SummarizedExperiment::assays(D_norm_loess)$intensities_norm)
    #    variant = Sys.info()[["sysname"]])
    expect_snapshot(SummarizedExperiment::rowData(D_norm_loess))
    expect_snapshot(SummarizedExperiment::colData(D_norm_loess))

    expect_snapshot(D_norm_lts)
    expect_snapshot(SummarizedExperiment::assays(D_norm_lts)$intensities_norm)
    #    variant = Sys.info()[["sysname"]])
    expect_snapshot(SummarizedExperiment::rowData(D_norm_lts))
    expect_snapshot(SummarizedExperiment::colData(D_norm_lts))
})

