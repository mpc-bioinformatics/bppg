
test_that("digestion of a FASTA file", {
  digested_proteins <- readRDS(test_path("digested_proteins_test.rds"))

  file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
  names(fasta) <- limma::strsplit2(names(fasta), "\\|")[,2]
  res <- digest_fasta(fasta)

  expect_equal(res, digested_proteins)
})


test_that("generation of an edgelist", {
  edgelist <- readRDS(test_path("edgelist_test.rds"))

  digested_proteins <- readRDS(test_path("digested_proteins_test.rds"))
  res <- generate_edgelist(digested_proteins)

  expect_equal(res, edgelist)
})

