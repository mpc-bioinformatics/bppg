library(bppg)

library(seqinr)
file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
res <- digest_fasta(fasta)

saveRDS(res, file = "inst/extdata/digested_proteins_test.rds")
saveRDS(res, file = "tests/testthat/digested_proteins_test.rds")
