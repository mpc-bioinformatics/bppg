



test_that("collapsing of edgelists", {
    edgelist_coll_pept_prot <- readRDS(testthat::test_path("testfiles/edgelist_coll_pept_prot_test.rds"))

    edgelist <- readRDS(testthat::test_path("testfiles/edgelist_test.rds"))
    res <- bppg:::.collapseEdgelist(edgelist,
                                    collapse_protein_nodes = TRUE,
                                    collapse_peptide_nodes = TRUE)

    expect_equal(res, edgelist_coll_pept_prot)


    edgelist_coll_prot <- readRDS(testthat::test_path("testfiles/edgelist_coll_prot_test.rds"))

    edgelist <- readRDS(testthat::test_path("testfiles/edgelist_test.rds"))
    res2 <- bppg:::.collapseEdgelist(edgelist,
                                     collapse_protein_nodes = TRUE,
                                     collapse_peptide_nodes = FALSE)

    expect_equal(res2, edgelist_coll_prot)
})

