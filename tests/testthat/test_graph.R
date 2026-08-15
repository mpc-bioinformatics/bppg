test_that("plotBipartiteGraph ", {

    graphs <- readRDS(
        testthat::test_path(
            "testfiles",
            "graphs_coll_pept_prot_test.rds"
        )
    )

    expect_true(length(graphs) > 0)

    for (i in seq_along(graphs)) {

        p <- plotBipartiteGraph(
            graphs[[i]],
            three_shapes = TRUE,
            node_labels_proteins = "letters",
            node_labels_peptides = "numbers"
        )

        expect_s3_class(p, "ggplot")
    }
})

test_that("save", {

    graphs <- readRDS(
        testthat::test_path(
            "testfiles",
            "graphs_coll_pept_prot_test.rds"
        )
    )

    output_dir <- testthat::test_path("testplots")
    dir.create(output_dir, showWarnings = FALSE)

    for (i in seq_along(graphs)) {

        out_file <- file.path(
            output_dir,
            paste0("graph_", i, ".jpg")
        )

        plotBipartiteGraph(
            graphs[[i]],
            three_shapes = TRUE,
            node_labels_proteins = "letters",
            node_labels_peptides = "numbers",
            save_path = out_file
        )

        expect_true(file.exists(out_file))
    }
})