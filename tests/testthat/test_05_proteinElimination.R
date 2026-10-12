test_that("test proteinElimination", {

    file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
    graphs <- readRDS(file)
    G <- graphs[[1]][[4]]

    res <- proteinElimination(G = G)

    print(names(res))
    print(class(res$res_best$G[[1]]))
    print(is.null(res$G_plot))

    p <- plotBipartiteGraph(
        res$G_plot,
        three_shapes = TRUE
    )

    ggplot2::ggsave(
        filename = "protein_elimination.png",
        plot = p,
        width = 10,
        height = 8,
        dpi = 300
    )


    ## prepare result object for snapshot (cant deal with random igraph ids)
    res$protsOriginIDs <- names(res$protsOriginIDs)
    res$res_best$G <- lapply(res$res_best$G, igraph::as_edgelist)

    ## flooring error term here because of small Windows/Linux differences
    res$resDF$error <- floor(res$resDF$error * 10^4) / 10^4

    expect_snapshot(res)

})
