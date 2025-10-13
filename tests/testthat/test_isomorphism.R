test_that("test .directBipartiteGraph", {
  library(igraph)

  # Create bipartite graph
  edges <- c("prot_1", "pep_1",
             "pep_3" , "prot_1",
             "pep_2" , "prot_1",
             "prot_2", "pep_2",
             "prot_2", "pep_3",
             "prot_3", "pep_2",
             "pep_3" , "prot_3",
             "prot_3", "pep_4")
  types <- c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE)
  bipartite_graph <- igraph::make_graph(edges = edges, directed = FALSE)
  V(bipartite_graph)$type <- types

  # from peptide to protein
  directed_graph <- bppg:::.directBipartiteGraph(bipartite_graph, from_type = FALSE)
  pred_res <- cbind(c("pep_1", "pep_3", "pep_2", "pep_2", "pep_3", "pep_2", "pep_3", "pep_4"),
                    c("prot_1", "prot_1", "prot_1", "prot_2", "prot_2", "prot_3", "prot_3", "prot_3"))

  # from protein to peptide
  directed_graph2 <- bppg:::.directBipartiteGraph(bipartite_graph, from_type = TRUE)
  pred_res2 <- cbind(c("prot_1", "prot_1", "prot_1", "prot_2", "prot_2", "prot_3", "prot_3", "prot_3"),
                    c("pep_1", "pep_3", "pep_2", "pep_2", "pep_3", "pep_2", "pep_3", "pep_4"))


  expect_true(igraph::is_directed(directed_graph))
  expect_equal(igraph::as_edgelist(directed_graph), pred_res)
  expect_true(igraph::is_directed(directed_graph2))
  expect_equal(igraph::as_edgelist(directed_graph2), pred_res2)

})

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

  # M + 1 graph
  M4 <- matrix(c(1,1,0,0, 1, 0,0,1,1), nrow = 3, byrow = TRUE)
  G4 <- igraph::graph_from_biadjacency_matrix(M4)

  # G5 is same as G4 but with changes role of peptides/proteins
  G5 <- G4
  igraph::V(G5)$type <- !igraph::V(G4)$type

  # G6 is same as G4 but rows and columns are in different order
  M6 <- M4
  M6 <- M4[c(3,2,1), c(2,1,3)]
  G6 <- igraph::graph_from_biadjacency_matrix(M6)

  expect_true(bppg:::.isomorphicBipartite(G, G)) # N = N
  expect_false(bppg:::.isomorphicBipartite(G, G2)) # N != W
  expect_false(bppg:::.isomorphicBipartite(G, G3)) # N != M
  expect_false(bppg:::.isomorphicBipartite(G2, G3)) # W != M
  expect_false(bppg:::.isomorphicBipartite(G5, G4))
  expect_true(bppg:::.isomorphicBipartite(G6, G4))

})

test_that("generation of a prototype list", {
    library(igraph)

    # N shape graph
    M <- matrix(c(1,1,0,1), nrow = 2, byrow = TRUE)
    G <- igraph::graph_from_biadjacency_matrix(M)

    # N shape graph type 2 (different order of nodes)
    M2 <- M[c(2,1), c(2,1)]
    G2 <- igraph::graph_from_biadjacency_matrix(M2)

    # W shape graph
    M3 <- matrix(c(1,1,0,0,1,1), nrow = 2, byrow = TRUE)
    G3 <- igraph::graph_from_biadjacency_matrix(M3)

    # W shape graph type 2 (different order of nodes)
    M4 <- M3[, c(2,1,3)]
    G4 <- igraph::graph_from_biadjacency_matrix(M4)

    # M shape graph
    G5 <- G3
    igraph::V(G5)$type <- !igraph::V(G5)$type

    G_list <- list(G3, G3, G5, G, G, G2)
    # just to test the case when the last element will form its own class:
    G_list2 <- list(G3, G3, G, G, G2, G5)

    proto_list <- bppg:::.generatePrototypeList(G_list, sort_by_nr_edges = TRUE)
    proto_list2 <- bppg:::.generatePrototypeList(G_list, sort_by_nr_edges = FALSE)
    proto_list3 <- bppg:::.generatePrototypeList(G_list2, sort_by_nr_edges = FALSE)

    expect_type(proto_list, "list")
    expect_type(proto_list2, "list")

    expect_length(proto_list, 2)
    expect_length(proto_list2, 2)

    expect_equal(names(proto_list), c("graphs", "counter"))

    expect_length(proto_list$graphs, 3)
    expect_length(proto_list2$graphs, 3)
    expect_length(proto_list3$graphs, 3)

    expect_equal(proto_list$counter, c(3,2,1))
    expect_equal(proto_list2$counter, c(2,1,3))
    expect_equal(proto_list3$counter, c(2,3,1))

})
