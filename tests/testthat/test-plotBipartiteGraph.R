test_that("plot a bipartite graph", {

  # Skip this test on continuous integration systems like GitHub Actions
  # The function expect_snapshot_file is otherwise too strict
  # and there is no way to get a few pixel of tolerance
  testthat::skip_on_ci()

  # Create a temporary directory so no permanent files are put on a package users directory
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  # Create bipartite graph
  edges <- c("prot_1", "pep_1",
             "prot_1", "pep_2",
             "prot_1", "pep_3",
             "prot_2", "pep_2",
             "prot_2", "pep_3",
             "prot_3", "pep_2",
             "prot_3", "pep_3",
             "prot_3", "pep_4")
  types <- c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE)
  bipartite_graph <- igraph::make_graph(edges = edges, directed = FALSE)
  V(bipartite_graph)$type <- types


  # Plot graph
  png(filename = file.path(temp_dir, "bipartitGraph.png"))
  plotBipartiteGraph(bipartite_graph)
  dev.off()


  expect_snapshot_file(path = file.path(temp_dir, "bipartitGraph.png"), name = "plotBipartitGraph.png" )

})
