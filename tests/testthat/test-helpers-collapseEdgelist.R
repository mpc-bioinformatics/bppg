
test_that("collapsing of edgelists", {
    edgelist <- readRDS(testthat::test_path("testfiles/edgelist_test.rds"))
    edgelist_collapsed <- bppg:::.collapseEdgelist(edgelist,
                                    collProtNodes  = TRUE,
                                    collPeptNodes  = TRUE)

    edgelist_collapsed2 <- bppg:::.collapseEdgelist(edgelist,
                                     collProtNodes  = TRUE,
                                     collPeptNodes  = FALSE)

    # this shouldnt change the edgelist
    edgelist_collapsed3 <- bppg:::.collapseEdgelist(edgelist,
                                                    collProtNodes  = FALSE,
                                                    collPeptNodes  = FALSE)

    expect_snapshot(edgelist_collapsed)
    expect_snapshot(edgelist_collapsed2)
    expect_equal(edgelist_collapsed3, edgelist)
})

