

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



test_that("result export", {
    file <- system.file("extdata", "resultsList.rds", package = "bppg")
    resultsList <- readRDS(file)

    temp_dir <- tempfile(pattern = "test_dir")
    dir.create(temp_dir)
    on.exit(unlink(temp_dir, recursive = TRUE))

    SE <- combineComparisons(resultsList)
    exportSE(SE, file.path(temp_dir,"results.xlsx"))

    expect_true(file.exists(file.path(temp_dir, "results.xlsx")))
})
