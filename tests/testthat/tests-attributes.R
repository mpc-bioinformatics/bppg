

test_that("test add_uniqueness_attributes", {

  library(igraph)

  # W shaped graph
  M <- matrix(c(1,1,0,0,1,1), nrow = 2, byrow = TRUE)
  G <- igraph::graph_from_biadjacency_matrix(M)
  G_new <- bppg::add_uniqueness_attributes(G)

  expect_equal(igraph::V(G_new)$nr_unique_peptides, c(NA, NA, 0, 0, 0))
  expect_equal(igraph::V(G_new)$nr_shared_peptides, c(NA, NA, 1, 2, 1))
  expect_equal(igraph::V(G_new)$uniqueness, c(FALSE, FALSE, NA, NA, NA))



  # M shaped graph
  V(G)$type <- !V(G)$type
  G_new2 <- bppg::add_uniqueness_attributes(G)

  expect_equal(igraph::V(G_new2)$nr_unique_peptides, c(1, 1, NA, NA, NA))
  expect_equal(igraph::V(G_new2)$nr_shared_peptides, c(1, 1, NA, NA, NA))
  expect_equal(igraph::V(G_new2)$uniqueness, c(NA, NA, TRUE, FALSE, TRUE))

})
