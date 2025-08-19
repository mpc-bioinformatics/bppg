test_that("test calculate_proteinnode_info", {

  file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

  g1 <- generateGraphsFromFASTA(fasta = fasta,
                                   collapse_protein_nodes = FALSE,
                                   collapse_peptide_nodes = TRUE)
  g2 <- generateGraphsFromFASTA(fasta = fasta,
                                   collapse_protein_nodes = TRUE,
                                   collapse_peptide_nodes = FALSE)
  g3 <- generateGraphsFromFASTA(fasta = fasta,
                                   collapse_protein_nodes = FALSE,
                                   collapse_peptide_nodes = FALSE)

  graphs <- list(s1_s2 = g1, s1_s3 = g2, s2_s3 = g3)

  D <- calculate_proteinnode_info(graphs)

  expect_equal(nrow(D), 30)
  expect_equal(ncol(D), 7)
  expect_equal(unlist(D[11, ]), c(accessions = "tr|A0A494C0S0|A0A494C0S0_HUMAN",
                                  comparison = "s1_s3",
                                  graphID = "1",
                                  ind_within_graph = "1",
                                  nr_peptides = "288",
                                  nr_unique_peptides = "31",
                                  nr_shared_peptides = "257"))
  expect_equal(unlist(D[22, ]), c(accessions = "tr|H0Y9R5|H0Y9R5_HUMAN",
                                  comparison = "s2_s3",
                                  graphID = "1",
                                  ind_within_graph = "2",
                                  nr_peptides = "44",
                                  nr_unique_peptides = "3",
                                  nr_shared_peptides = "41"))
})
