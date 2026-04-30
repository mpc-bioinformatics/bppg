test_that("plot a bipartite graph", {
    # Skip this test on continuous integration systems like GitHub Actions
    # The function expect_snapshot_file is otherwise too strict
    # and there is no way to get a few pixel of tolerance
    testthat::skip_on_ci()

    library(igraph)

    # Create a temporary directory so no permanent files are put on a package users directory
    temp_dir <- tempfile(pattern = "test_dir")
    dir.create(temp_dir)
    on.exit(unlink(temp_dir, recursive = TRUE))

    # Create bipartite graph
    edges <- c("prot_1", "pep_1",
                "prot_1", "pep_2",
                "prot_1", "pep_3",
                "prot_2", "pep_2",
                "prot_2", "pep_3",
                "prot_3", "pep_2",
                "prot_3", "pep_3",
                "prot_3", "pep_4")
    types <- c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE)
    G <- igraph::make_graph(edges = edges, directed = FALSE)
    V(G)$type <- types

    V(G)$pep_ratio <- rep(NA, 7)
    V(G)$pep_ratio[!V(G)$type] <- c(0.177569, 0.1, 0.5, 0.6)

    V(G)$pep_ratio_aggr <- rep(NA, 7)
    V(G)$pep_ratio_aggr[!V(G)$type] <- rev(c(0.177569, 0.1, 0.5, 0.6))


    # Plot graph (default settings)
    png(filename = file.path(temp_dir, "bipartitGraph.png"))
    plotBipartiteGraph(G)
    dev.off()

    png(filename = file.path(temp_dir, "bipartitGraph2.png"))
    plotBipartiteGraph(G,
                        node_labels_proteins = "accessions",
                        node_labels_peptides = "pep_ratios",
                        useCanonicalPermutation = TRUE,
                        three_shapes = TRUE,
                        use_edge_attributes = TRUE,
                        round_digits = 2)
    dev.off()

    png(filename = file.path(temp_dir, "bipartitGraph3.png"))
    plotBipartiteGraph(G,
                        node_labels_proteins = "numbers_noord",
                        node_labels_peptides = "pep_ratios_aggr",
                        useCanonicalPermutation = TRUE,
                        three_shapes = TRUE,
                        use_edge_attributes = TRUE,
                        round_digits = 2)
    dev.off()

    png(filename = file.path(temp_dir, "bipartitGraph4.png"))
    plotBipartiteGraph(G,
                        node_labels_proteins = "numbers_noord",
                        node_labels_peptides = "")
    dev.off()


    expect_snapshot_file(path = file.path(temp_dir, "bipartitGraph.png"), name = "plotBipartitGraph.png", variant = Sys.info()[["sysname"]])
    expect_snapshot_file(path = file.path(temp_dir, "bipartitGraph2.png"), name = "plotBipartitGraph2.png", variant = Sys.info()[["sysname"]])
    expect_snapshot_file(path = file.path(temp_dir, "bipartitGraph3.png"), name = "plotBipartitGraph3.png", variant = Sys.info()[["sysname"]])
    expect_snapshot_file(path = file.path(temp_dir, "bipartitGraph4.png"), name = "plotBipartitGraph4.png", variant = Sys.info()[["sysname"]])

})

## only test file existence, not snapshots which can be tested also on github actions
test_that("plot a bipartite graph file existence", {

    library(igraph)

    # Create a temporary directory so no permanent files are put on a package users directory
    temp_dir <- tempfile(pattern = "test_dir")
    dir.create(temp_dir)
    on.exit(unlink(temp_dir, recursive = TRUE))

    # Create bipartite graph
    edges <- c("prot_1", "pep_1",
                "prot_1", "pep_2",
                "prot_1", "pep_3",
                "prot_2", "pep_2",
                "prot_2", "pep_3",
                "prot_3", "pep_2",
                "prot_3", "pep_3",
                "prot_3", "pep_4")
    types <- c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE)
    G <- igraph::make_graph(edges = edges, directed = FALSE)
    V(G)$type <- types

    V(G)$pep_ratio <- rep(NA, 7)
    V(G)$pep_ratio[!V(G)$type] <- c(0.177569, 0.1, 0.5, 0.6)

    V(G)$pep_ratio_aggr <- rep(NA, 7)
    V(G)$pep_ratio_aggr[!V(G)$type] <- rev(c(0.177569, 0.1, 0.5, 0.6))


    # Plot graph (default settings)
    png(filename = file.path(temp_dir, "bipartitGraph.png"))
    plotBipartiteGraph(G)
    dev.off()

    png(filename = file.path(temp_dir, "bipartitGraph2.png"))
    plotBipartiteGraph(G,
                        node_labels_proteins = "accessions",
                        node_labels_peptides = "pep_ratios",
                        useCanonicalPermutation = TRUE,
                        three_shapes = TRUE,
                        use_edge_attributes = TRUE,
                        round_digits = 2)
    dev.off()

    png(filename = file.path(temp_dir, "bipartitGraph3.png"))
    plotBipartiteGraph(G,
                        node_labels_proteins = "numbers_noord",
                        node_labels_peptides = "pep_ratios_aggr",
                        useCanonicalPermutation = TRUE,
                        three_shapes = TRUE,
                        use_edge_attributes = TRUE,
                        round_digits = 2)
    dev.off()

    png(filename = file.path(temp_dir, "bipartitGraph4.png"))
    plotBipartiteGraph(G,
                        node_labels_proteins = "numbers_noord",
                        node_labels_peptides = "")
    dev.off()


    expect_true(file.exists(file.path(temp_dir, "bipartitGraph.png")))
    expect_true(file.exists(file.path(temp_dir, "bipartitGraph2.png")))
    expect_true(file.exists(file.path(temp_dir, "bipartitGraph3.png")))
    expect_true(file.exists(file.path(temp_dir, "bipartitGraph4.png")))

})
