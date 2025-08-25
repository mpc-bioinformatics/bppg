## tests for isomorphisms

test_that(".isomorphicBipartite works as intended", {
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

  expect_true(bppg::.isomorphicBipartite(G, G))
  expect_false(bppg::.isomorphicBipartite(G, G2))
  expect_false(bppg::.isomorphicBipartite(G, G3))


  # M + 1 graph
  M4 <- matrix(c(1,1,0,0, 1, 0,0,1,1), nrow = 3, byrow = TRUE)
  G4 <- igraph::graph_from_biadjacency_matrix(M4)

  G5 <- G4
  igraph::V(G5)$type <- !igraph::V(G4)$type

  expect_false(bppg::.isomorphicBipartite(G4, G5))

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

  proto_list <- bppg::.generatePrototypeList(G_test, sort_by_nr_edges = TRUE)
  proto_list2 <- bppg::.generatePrototypeList(G_test, sort_by_nr_edges = FALSE)


  expect_type(proto_list, "list")
  expect_type(proto_list2, "list")

  expect_length(proto_list, 2)
  expect_length(proto_list2, 2)

  expect_length(proto_list$graphs, 3)
  expect_length(proto_list2$graphs, 3)

  expect_equal(proto_list$counter, c(3,2,1))
  expect_equal(proto_list2$counter, c(2,3,1))


})



test_that("test plotting of the top10 isomorphism classes", {

  # Create a temporary directory so no permanent files are put on a package users directory
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))


  # Make graphs
  G_list <- list()
  set.seed(2)
  for (i in 4:15) {
    nr_repeats <- sample(1:3, size = 1)

    # Make easy edgelist with edges 1-2, 2-3, ..., (i-1)-i
    edge_list <- matrix(rbind(1:(i-1), 2:i), ncol = 2, byrow = TRUE)
    G <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
    V(G)$name <- 1:i
    V(G)$type <- (as.numeric(V(G)$name) %% 2 != 0) # True for odd, False for even

    G_list <- c(G_list, replicate(nr_repeats, G, simplify = FALSE))
  }

  proto_list <- bppg::.generatePrototypeList(G_list, sort_by_nr_edges = TRUE)

  plot_top10_iso_classes(prototypelist = proto_list,
                         out_file = file.path(temp_dir, "proto_list.pdf"))

  expect_true(file.exists(file.path(temp_dir, "proto_list.pdf")))

})

