

test_that("generation of graphs from edgelist", {
    library(igraph)

    graphs_coll_pept_prot <- readRDS(testthat::test_path("testfiles/graphs_coll_pept_prot_test.rds"))
    graphs_coll_prot <- readRDS(testthat::test_path("testfiles/graphs_coll_prot_test.rds"))

    # with collapsing of peptide and protein nodes
    edgelist_coll_pept_prot <- readRDS(testthat::test_path("testfiles/edgelist_coll_pept_prot_test.rds"))

    res <- bppg:::.generateGraphsFromEdgelist(edgelist_coll_pept_prot)
    expect_true(bppg:::.isomorphicBipartite(res[[1]], graphs_coll_pept_prot[[1]]))
    expect_true(bppg:::.isomorphicBipartite(res[[2]], graphs_coll_pept_prot[[2]]))
    expect_true(bppg:::.isomorphicBipartite(res[[3]], graphs_coll_pept_prot[[3]]))

    # with collapsing of only protein nodes
    edgelist_coll_prot <- readRDS(test_path("testfiles/edgelist_coll_prot_test.rds"))
    res2 <- bppg:::.generateGraphsFromEdgelist(edgelist_coll_prot)
    expect_true(bppg:::.isomorphicBipartite(res2[[1]], graphs_coll_prot[[1]]))
    expect_true(bppg:::.isomorphicBipartite(res2[[2]], graphs_coll_prot[[2]]))
    expect_true(bppg:::.isomorphicBipartite(res2[[3]], graphs_coll_prot[[3]]))

})
