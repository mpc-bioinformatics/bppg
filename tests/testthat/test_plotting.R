test_that("plot a bipartite graph", {
    # Skip this test on continuous integration systems like GitHub Actions
    # The function expect_snapshot_file is otherwise too strict
    # and there is no way to get a few pixel of tolerance
    testthat::skip_on_ci()

    # Create a temporary directory so no permanent files are put on a package users directory
    temp_dir <- tempfile(pattern = "test_dir")
    dir.create(temp_dir)
    on.exit(unlink(temp_dir, recursive = TRUE))

    # Create bipartite graph
    file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
    G <- readRDS(file)[[1]][[3]]
    file <- system.file("extdata", "quantGraphs_collpept.rds", package = "bppg")
    G_coll <- readRDS(file)[[1]][[3]]

    # Plot graph (default settings)
    p1 <- plotBipartiteGraph(G_coll)
    ggplot2::ggsave(plot = p1,
        filename = file.path(temp_dir, "bipartitGraph.png"),
        width = 10, height = 7, dpi = 300, device = "png"
    )
    
    p2 <- plotBipartiteGraph(G,
        node_labels_proteins = "accessions",
        node_labels_peptides = "pep_logRatios",
        round_digits = 2)
    ggplot2::ggsave(plot = p2,
        filename = file.path(temp_dir, "bipartitGraph2.png"),
        width = 10, height = 7, dpi = 300, device = "png"
    )

    p3 <- plotBipartiteGraph(G_coll,
        node_labels_proteins = "numbers_noord",
        node_labels_peptides = "pep_ratios_mean",
        round_digits = 2)
    ggplot2::ggsave(plot = p3,
        filename = file.path(temp_dir, "bipartitGraph3.png"),
        width = 10, height = 7, dpi = 300, device = "png"
    )

    p4 <- plotBipartiteGraph(G_coll,
        node_labels_proteins = "numbers_noord",
        node_labels_peptides = "")
    ggplot2::ggsave(plot = p4,
        filename = file.path(temp_dir, "bipartitGraph4.png"),
        width = 10, height = 7, dpi = 300, device = "png"
    )


    expect_snapshot_file(path = file.path(temp_dir, "bipartitGraph.png"),
                         name = "plotBipartitGraph.png", variant = Sys.info()[["sysname"]])
    expect_snapshot_file(path = file.path(temp_dir, "bipartitGraph2.png"),
                         name = "plotBipartitGraph2.png", variant = Sys.info()[["sysname"]])
    expect_snapshot_file(path = file.path(temp_dir, "bipartitGraph3.png"),
                         name = "plotBipartitGraph3.png", variant = Sys.info()[["sysname"]])
    expect_snapshot_file(path = file.path(temp_dir, "bipartitGraph4.png"),
                         name = "plotBipartitGraph4.png", variant = Sys.info()[["sysname"]])

})

## only test file existence, not snapshots which can be tested also on github actions
test_that("plot a bipartite graph file existence", {

    library(igraph)

    # Create a temporary directory so no permanent files are put on a package users directory
    temp_dir <- tempfile(pattern = "test_dir")
    dir.create(temp_dir)
    on.exit(unlink(temp_dir, recursive = TRUE))

    # Create bipartite graph
    file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
    G <- readRDS(file)[[1]][[3]]
    file <- system.file("extdata", "quantGraphs_collpept.rds", package = "bppg")
    G_coll <- readRDS(file)[[1]][[3]]

    # Plot graph (default settings)
    plotBipartiteGraph(G_coll,
        output_path = file.path(temp_dir, "bipartitGraph.png"),
        width = 10, height = 7, dpi = 300, device = "png")

    plotBipartiteGraph(G,
        node_labels_proteins = "accessions",
        node_labels_peptides = "pep_logRatios",
        output_path = file.path(temp_dir, "bipartitGraph2.png"),
        width = 10, height = 7, dpi = 300, device = "png")

    plotBipartiteGraph(G_coll,
        node_labels_proteins = "numbers_noord",
        node_labels_peptides = "pep_ratios_mean",
        round_digits = 2,
        output_path = file.path(temp_dir, "bipartitGraph3.png"),
        width = 10, height = 7, dpi = 300, device = "png")

    plotBipartiteGraph(G_coll,
        node_labels_proteins = "numbers_noord",
        node_labels_peptides = "",
        output_path = file.path(temp_dir, "bipartitGraph4.png"),
        width = 10, height = 7, dpi = 300, device = "png")



    expect_true(file.exists(file.path(temp_dir, "bipartitGraph.png")))
    expect_true(file.exists(file.path(temp_dir, "bipartitGraph2.png")))
    expect_true(file.exists(file.path(temp_dir, "bipartitGraph3.png")))
    expect_true(file.exists(file.path(temp_dir, "bipartitGraph4.png")))

})
