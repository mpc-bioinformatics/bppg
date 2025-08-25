test_that("test .directBipartiteGraph", {

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


  directed_graph <- .directBipartiteGraph(bipartite_graph, from_type = FALSE)
  pred_res <- cbind(c("pep_1", "pep_3", "pep_2", "pep_2", "pep_3", "pep_2", "pep_3", "pep_4"),
                    c("prot_1", "prot_1", "prot_1", "prot_2", "prot_2", "prot_3", "prot_3", "prot_3"))
  expect_equal(as_edgelist(directed_graph), pred_res)


  directed_graph <- .directBipartiteGraph(bipartite_graph, from_type = TRUE)
  pred_res <- cbind(c("prot_1", "prot_1", "prot_1", "prot_2", "prot_2", "prot_3", "prot_3", "prot_3"),
                    c("pep_1", "pep_3", "pep_2", "pep_2", "pep_3", "pep_2", "pep_3", "pep_4"))
  expect_equal(as_edgelist(directed_graph), pred_res)

})


test_that("test .convertToBipartiteGraph", {

  M <- matrix(c(1,0,0,1,1,1), nrow = 2, byrow = TRUE)

  g <- bppg::.convertToBipartiteGraph(M)

  #pred_res <- cbind(c(1,2,2,2), c(3,3,4,5))

  expect_equal(as_edgelist(g), cbind(c(1,2,2,2), c(3,3,4,5)))
  expect_equal(V(g)$type, c(FALSE, FALSE, TRUE, TRUE, TRUE))

})
