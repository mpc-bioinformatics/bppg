## TODO add test for new contracting functions

test_that("generation of graphs from edgelist", {
    library(igraph)

    graphs_coll_pept_prot <- readRDS(testthat::test_path("testfiles/graphs_coll_pept_prot_test.rds"))
    graphs_coll_prot <- readRDS(testthat::test_path("testfiles/graphs_coll_prot_test.rds"))

    # with collapsing of peptide and protein nodes
    edgelist_coll_pept_prot <- readRDS(testthat::test_path("testfiles/edgelist_coll_pept_prot_test.rds"))

    res <- bppg::generateGraphsFromEdgelist(edgelist_coll_pept_prot)


    # with collapsing of only protein nodes
    edgelist_coll_prot <- readRDS(test_path("testfiles/edgelist_coll_prot_test.rds"))
    res2 <- bppg::generateGraphsFromEdgelist(edgelist_coll_prot)

    expect_true(bppg:::.isomorphicBipartite(res[[1]], graphs_coll_pept_prot[[1]]))
    expect_true(bppg:::.isomorphicBipartite(res[[2]], graphs_coll_pept_prot[[2]]))
    expect_true(bppg:::.isomorphicBipartite(res[[3]], graphs_coll_pept_prot[[3]]))

    expect_true(bppg:::.isomorphicBipartite(res2[[1]], graphs_coll_prot[[1]]))
    expect_true(bppg:::.isomorphicBipartite(res2[[2]], graphs_coll_prot[[2]]))
    expect_true(bppg:::.isomorphicBipartite(res2[[3]], graphs_coll_prot[[3]]))

    # NOTE: snapshots don't work here because of random graph ids
})

test_that("test generateQuantGraphs", {
    # Create a temporary directory so no permanent files are put on a package users directory
    temp_dir <- tempfile(pattern = "test_dir")
    dir.create(temp_dir)
    on.exit(unlink(temp_dir, recursive = TRUE))

    # Create a ratio table and edgelist
    set.seed(8)
    ratio_table <- data.frame(peptides = paste0("pep_", 1:10),
                                ratio_sample1_sample2 = round(runif(10, min = 0.9, max = 1.1), digits = 3),
                                ratio_sample1_sample3 = round(runif(10, min = 0.9, max = 1.1), digits = 3),
                                ratio_sample2_sample3 = round(runif(10, min = 0.9, max = 1.1), digits = 3))
    for (i in 2:4) {
        ratio_table[sample(1:10, size = 2), i] <- NA # Insert some NAs
    }

    proteins <- rep(paste0("prot_", 1:5), times = c(4,2,3,4,4))
    peptides <- c(paste0("pep_", 1:4), paste0("pep_", 3:4), paste0("pep_", 5:7), paste0("pep_", 7:10), paste0("pep_", 7:10))
    edgelist <- data.frame(protein = proteins, peptide = peptides)

    rownames(ratio_table) <- ratio_table$peptides

    expData <- SummarizedExperiment::SummarizedExperiment(
        assays = list(logRatios = ratio_table[, -1]),
        rowData = data.frame(peptides = rownames(ratio_table)),
        colData = data.frame(comparison = colnames(ratio_table[, -1]))
    )    
    # Compute function
    graphs <- bppg::generateQuantGraphs(exp_peptide_ratios = expData,
                                    fasta_edgelist = edgelist,
                                    outpath = temp_dir,
                                    seq_column = "peptides",
                                    collProtNodes = TRUE,
                                    collPeptNodes = TRUE,
                                    suffix = "")

    # Check result attributes
    expect_true(file.exists(file.path(temp_dir, "edgelist_filtered_.xlsx")))
    expect_equal(unname(lapply(graphs, length)), list(3,2,2))
    expect_equal(names(graphs), c("sample1_sample2", "sample1_sample3", "sample2_sample3"))

    for (i in 1:3) {
        for (j in seq_along(graphs[[i]])) {
        expect_snapshot(igraph::as_edgelist(graphs[[i]][[j]]))
        expect_snapshot(igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio"))
        }
    }

})
