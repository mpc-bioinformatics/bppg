## tests for isomorphisms

test_that("isomorphic_bipartite works as intended", {
  library(igraph)

  # N shape graph
  M <- matrix(c(1,1,0,1), nrow = 2, byrow = TRUE)
  G <- igraph::graph_from_biadjacency_matrix(M)

  # W shape graph
  M2 <- matrix(c(1,1,0,0,1,1), nrow = 2, byrow = TRUE)
  G2 <- igraph::graph_from_biadjacency_matrix(M2)

  # M shape graph
  G3 <- G2
  V(G3)$type <- !V(G3)$type

  expect_true(bppg::isomorphic_bipartite(G, G))
  expect_false(bppg::isomorphic_bipartite(G, G2))
  expect_false(bppg::isomorphic_bipartite(G, G3))


  # M + 1 graph
  M4 <- matrix(c(1,1,0,0, 1, 0,0,1,1), nrow = 3, byrow = TRUE)
  G4 <- igraph::graph_from_biadjacency_matrix(M4)

  G5 <- G4
  igraph::V(G5)$type <- !igraph::V(G4)$type

  expect_false(bppg::isomorphic_bipartite(G4, G5))

})



test_that("generation of a prototype list", {
  library(igraph)

  # N shape graph
  M <- matrix(c(1,1,0,1), nrow = 2, byrow = TRUE)
  G <- igraph::graph_from_biadjacency_matrix(M)

  # W shape graph
  M2 <- matrix(c(1,1,0,0,1,1), nrow = 2, byrow = TRUE)
  G2 <- igraph::graph_from_biadjacency_matrix(M2)

  # M shape graph
  G3 <- G2
  V(G3)$type <- !V(G3)$type


  G_test <- list(G3, G3, G, G, G, G2)

  proto_list <- bppg::generatePrototypeList(G_test, sort_by_nr_edges = TRUE)
  proto_list2 <- bppg::generatePrototypeList(G_test, sort_by_nr_edges = FALSE)


  expect_type(proto_list, "list")
  expect_type(proto_list2, "list")

  expect_length(proto_list, 2)
  expect_length(proto_list2, 2)

  expect_length(proto_list$graphs, 3)
  expect_length(proto_list2$graphs, 3)

  expect_equal(proto_list$counter, c(3,2,1))
  expect_equal(proto_list2$counter, c(2,3,1))


})


