test_that("test .addUniquenessAttributes", {

  library(igraph)

  # W shaped graph
  M <- matrix(c(1,1,0,0,1,1), nrow = 2, byrow = TRUE)
  G <- igraph::graph_from_biadjacency_matrix(M)
  G_new <- bppg:::.addUniquenessAttributes(G)

  # M shaped graph
  V(G)$type <- !V(G)$type
  G_new2 <- bppg:::.addUniquenessAttributes(G)

  expect_equal(igraph::V(G_new)$nr_unique_peptides, c(NA, NA, 0, 0, 0))
  expect_equal(igraph::V(G_new)$nr_shared_peptides, c(NA, NA, 1, 2, 1))
  expect_equal(igraph::V(G_new)$uniqueness, c(FALSE, FALSE, NA, NA, NA))

  expect_equal(igraph::V(G_new2)$nr_unique_peptides, c(1, 1, NA, NA, NA))
  expect_equal(igraph::V(G_new2)$nr_shared_peptides, c(1, 1, NA, NA, NA))
  expect_equal(igraph::V(G_new2)$uniqueness, c(NA, NA, TRUE, FALSE, TRUE))

})

test_that("geometric mean", {

  x <- c(1,5,4.5,2)
  n <- length(x)

  res <- exp(mean(log(x)))
  res_prod <- prod(x)^(1/n)

  expect_equal(res, bppg:::.geomMean(x, useprod = FALSE))
  expect_equal(res_prod, bppg:::.geomMean(x, useprod = TRUE))
})