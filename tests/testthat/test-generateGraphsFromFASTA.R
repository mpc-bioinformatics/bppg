
test_that("generate graphs from fasta",{
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))


  file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

  res <- bppg::generateGraphsFromFASTA(fasta = fasta,
                                       save_intermediate = TRUE,
                                       result_path = paste0(temp_dir, "\\"))

  expect_snapshot(igraph::as_edgelist(res[[1]]))
  expect_snapshot(igraph::as_edgelist(res[[2]]))
  expect_snapshot(igraph::as_edgelist(res[[3]]))

  expect_true(file.exists(file.path(temp_dir, "edgelist_.txt")))
  expect_true(file.exists(file.path(temp_dir, "edgelist_collprotpept_.txt")))
  expect_true(file.exists(file.path(temp_dir, "subgraphs_collprotpept_.rds")))

})


