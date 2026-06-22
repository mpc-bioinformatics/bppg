

test_that("result aggregation", {
    file <- system.file("extdata", "resultsList.rds", package = "bppg")
    resultsList <- readRDS(file)

    X <- combineComparisons(resultsList)

    expect_snapshot(X)
    expect_snapshot(SummarizedExperiment::assays(X)$"1_2")
    expect_snapshot(SummarizedExperiment::rowData(X))
    expect_snapshot(SummarizedExperiment::colData(X))

    ## Add test for protOrigin im combined result
})
