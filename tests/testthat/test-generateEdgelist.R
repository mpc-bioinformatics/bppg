
test_that("generation of an edgelist", {
    edgelist <- readRDS(testthat::test_path("testfiles/edgelist_test.rds"))

    digested_proteins <- readRDS(testthat::test_path("testfiles/digested_proteins_test.rds"))
    res <- generateEdgelist(digested_proteins)

    expect_equal(res, edgelist)
})

