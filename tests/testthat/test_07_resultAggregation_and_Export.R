

test_that("result aggregation", {
    file <- system.file("extdata", "quantGraphsForTesting.rds", package = "bppg")
    graphs <- readRDS(file)

    RES <- list()
    for (comp in seq_along(graphs)) {
        comparison <- names(graphs)[comp]
        graphs_tmp <- graphs[[comp]]
        RES_tmp <- NULL
        for (i in 1:length(graphs_tmp)) {
            G <- graphs_tmp[[i]]
            res <- iterateOverCi(G, gridSize = 100)
            if (is.null(RES_tmp)) {
                RES_tmp <- automatedAnalysisIteratedCi(G, res)
            } else {
                RES_tmp <- rbind(RES_tmp, automatedAnalysisIteratedCi(G, res))
            }
        }
        RES <- c(RES, RES_tmp)
    }
    names(RES) <- names(graphs)

    resultsList <- RES

    X <- combineComparisons(resultsList)
     
    expect_snapshot(X)
    expect_snapshot(SummarizedExperiment::assays(X)$sample1_sample2)
    expect_snapshot(SummarizedExperiment::assays(X)$sample1_sample3)
    expect_snapshot(SummarizedExperiment::assays(X)$sample2_sample3)
    expect_snapshot(SummarizedExperiment::rowData(X))
    expect_snapshot(SummarizedExperiment::colData(X))
})
