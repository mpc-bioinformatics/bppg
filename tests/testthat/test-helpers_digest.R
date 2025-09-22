# test of digestFASTA() is in the file test-graph_generation_FASTA.R


test_that("test .digest2", {

  file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

  digested_proteins1 <- bppg:::.digest2(sequence = fasta[[1]],
                                enzyme = "trypsin",
                                missed = 0,
                                remove_initial_M = TRUE)
  expect_snapshot(digested_proteins1)

  digested_proteins2 <- bppg:::.digest2(sequence = fasta[[1]],
                                enzyme = "trypsin.strict",
                                missed = 2,
                                remove_initial_M = FALSE)
  expect_snapshot(digested_proteins2)

  ## test if undefined enzyme is being used
  expect_error(bppg:::.digest2(fasta[[1]], enzyme = "LysC"))

  expect_warning(bppg:::.digest2("ABCD")) # warning for no cleavage site
  expect_warning(bppg:::.digest2("ABCRD", missed = 2)) # warning for not possible 2 missed cleavage

})


test_that("digestion of a FASTA file", {
  digested_proteins <- readRDS(testthat::test_path("testfiles/digested_proteins_test.rds"))

  file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
  names(fasta) <- limma::strsplit2(names(fasta), "\\|")[,2]
  res <- digestFASTA(fasta)

  expect_equal(res, digested_proteins)
})



