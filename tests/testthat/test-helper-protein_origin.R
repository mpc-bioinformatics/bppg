test_that("test .getProteinOrigin", {

  file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

  accessions <- c(paste0(names(fasta)[1:3], collapse = ";"),
                  paste0(names(fasta)[4:5], collapse = ";"),
                  paste0(names(fasta)[6:10], collapse = ";"))

  # Make some accessions from the fasta contaminants/spike-ins
  # and split the rest between different organisms
  contaminants <- names(fasta)[1:2]
  spike_ins <- names(fasta)[3:4]
  organisms <- list(human = names(fasta)[5:7], yeast = names(fasta)[8:10])

  origin <- .getProteinOrigin(accessions = accessions,
                               contaminants = contaminants,
                               spike_ins = spike_ins,
                               organisms = organisms)

  res <- c(origin[[1]], origin[[2]], origin[[3]])

  expect_equal(res, c("Contaminant;Spike-in", "Spike-in;human", "human;yeast"))
})
