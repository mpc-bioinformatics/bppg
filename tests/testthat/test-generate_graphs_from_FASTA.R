
test_that("generate graphs from fasta",{

  file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

  res <- bppg::generateGraphsFromFASTA(fasta = fasta)

  expect_snapshot(igraph::as_edgelist(res[[1]]))
  expect_snapshot(igraph::as_edgelist(res[[2]]))
  expect_snapshot(igraph::as_edgelist(res[[3]]))

})



