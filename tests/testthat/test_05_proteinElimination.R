test_that("test proteinElimination", {

    file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
    graphs <- readRDS(file)
    G <- graphs[[1]][[4]]

    res <- proteinElimination(G = G)

    ## prepare result object for snapshot (cant deal with random igraph ids)
    res$protsOriginIDs <- names(res$protsOriginIDs)
    res$res_best$G <- lapply(res$res_best$G, igraph::as_edgelist)

    ## flooring error term here because of small Windows/Linux differences
    res$resDF$error <- floor(res$resDF$error * 10^5) / 10^5

    expect_snapshot(res)

})
