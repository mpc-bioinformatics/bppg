# test of digestFASTA() is in the file test-graph_generation_FASTA.R


test_that("test .digest2", {

  file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

  digested_proteins1 <- .digest2(sequence = fasta[[1]],
                                enzyme = "trypsin",
                                missed = 0,
                                remove_initial_M = TRUE)

  digested_proteins2 <- .digest2(sequence = fasta[[1]],
                                enzyme = "trypsin.strict",
                                missed = 2,
                                remove_initial_M = FALSE)


  expect_equal(nrow(digested_proteins1), 122)
  expect_equal(nrow(digested_proteins2), 381)

  expect_equal(digested_proteins1[1:5, 1], c("MHPDATDSGGAGPSPAR",
                                             "AAGAGGRPVSGFR",
                                             "GER",
                                             "RPESPGDAEAAAAAAPGAPGGR",
                                             "SWWKPVAVAALAAVALSFLGPGSGEAAGAAGLSSVLFR"))
  expect_equal(digested_proteins2[1:5, 1], c("MHPDATDSGGAGPSPAR",
                                             "AAGAGGR",
                                             "PVSGFR",
                                             "GER",
                                             "R"))
})
