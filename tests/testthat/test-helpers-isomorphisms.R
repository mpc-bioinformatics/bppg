## tests for isomorphisms

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




