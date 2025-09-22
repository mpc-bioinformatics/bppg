



test_that("collapsing of edgelists", {
    #edgelist_coll_pept_prot <- readRDS(testthat::test_path("testfiles/edgelist_coll_pept_prot_test.rds"))

    edgelist <- readRDS(testthat::test_path("testfiles/edgelist_test.rds"))
    edgelist_collapsed <- bppg:::.collapseEdgelist(edgelist,
                                    collProtNodes  = TRUE,
                                    collPeptNodes  = TRUE)
    expect_snapshot(edgelist_collapsed)

    edgelist_collapsed2 <- bppg:::.collapseEdgelist(edgelist,
                                     collProtNodes  = TRUE,
                                     collPeptNodes  = FALSE)
    expect_snapshot(edgelist_collapsed2)

    # this shouldnt change the edgelist
    edgelist_collapsed3 <- bppg:::.collapseEdgelist(edgelist,
                                                    collProtNodes  = FALSE,
                                                    collPeptNodes  = FALSE)
    expect_equal(edgelist_collapsed3, edgelist)
})

