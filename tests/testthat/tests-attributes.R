

test_that("test .addUniquenessAttributes", {

  library(igraph)

  # W shaped graph
  M <- matrix(c(1,1,0,0,1,1), nrow = 2, byrow = TRUE)
  G <- igraph::graph_from_biadjacency_matrix(M)
  G_new <- bppg::.addUniquenessAttributes(G)

  expect_equal(igraph::V(G_new)$nr_unique_peptides, c(NA, NA, 0, 0, 0))
  expect_equal(igraph::V(G_new)$nr_shared_peptides, c(NA, NA, 1, 2, 1))
  expect_equal(igraph::V(G_new)$uniqueness, c(FALSE, FALSE, NA, NA, NA))



  # M shaped graph
  V(G)$type <- !V(G)$type
  G_new2 <- bppg::.addUniquenessAttributes(G)

  expect_equal(igraph::V(G_new2)$nr_unique_peptides, c(1, 1, NA, NA, NA))
  expect_equal(igraph::V(G_new2)$nr_shared_peptides, c(1, 1, NA, NA, NA))
  expect_equal(igraph::V(G_new2)$uniqueness, c(NA, NA, TRUE, FALSE, TRUE))

})


test_that("test .addAveragePepRatio", {

  library(igraph)

  # Create graph
  edges <- c("prot_1", "pep_1",
             "prot_1", "pep_2",
             "prot_1", "pep_3",
             "prot_2", "pep_2",
             "prot_3", "pep_2",
             "prot_3", "pep_3",
             "prot_3", "pep_4")
  types <- c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE)
  pep_ratios <- c(NA, "1.05;0.9", NA, "0.95;1.1", NA, NA, "1.15")
  bipartite_graph <- igraph::make_graph(edges = edges, directed = FALSE)
  V(bipartite_graph)$type <- types
  V(bipartite_graph)$pep_ratio <- pep_ratios

  # Check attribute creation
  bipartite_graph <- .addAveragePepRatio(bipartite_graph)

  expect_equal(V(bipartite_graph)$nr_sequences,
               c(1, 2, 1, 2, 1, 1, 1))
  expect_equal(round(V(bipartite_graph)$pep_ratio_aggr, digits = 5),
               c(NA, 0.97211, NA, 1.02225, NA, NA, 1.15000))

})


test_that("test assign_protein_accessions", {

  # Read in fasta and some peptide sequences
  file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
  fasta_names <- sapply(fasta, function(protein) protein[[1]])

  seq <- c("RLSEFQNLHRK", "LINSNSDVEFLKQLRYQIVVEIIQATTISSFPQLK",
           "KLSECVPSLK", "SLYKEIQQCLVGNKGIEVFYK",
           "VPYTKVQLKELEK", "RISATTNLSERQVTIWFQNR",
           "SRAGSWDMDGLRADGGGAGGAPASSSSSSVAAAAASGQCR", "DQPQGSHFWK",
           "WQLGGSAIPLPPSHKFRK", "MEPGNYATLDGAKDIEGLLGAGGGR",
           "IPYSKGQLRELER", "EPGNYATLDGAKDIEGLLGAGGGR",
           "IAEGFSDTALIMVDNTK", "CRDPHHDYCEDWPEAQRISASLLDSR",
           "SYETLVDFDNHLDDIR")


  # The expected output
  expected_res <- c("sp|Q9H3E2|SNX25_HUMAN/tr|A0A494C0S0|A0A494C0S0_HUMAN",
                    "sp|Q9H3E2|SNX25_HUMAN/tr|A0A494C0S0|A0A494C0S0_HUMAN",
                    "sp|Q9H3E2|SNX25_HUMAN/tr|A0A494C0S0|A0A494C0S0_HUMAN",
                    "sp|Q9H3E2|SNX25_HUMAN/tr|A0A494C0S0|A0A494C0S0_HUMAN/tr|H0Y9R5|H0Y9R5_HUMAN",
                    "sp|P31276|HXC13_HUMAN",
                    "sp|P31271|HXA13_HUMAN/sp|P31276|HXC13_HUMAN",
                    "sp|P35453|HXD13_HUMAN",
                    "sp|P35453|HXD13_HUMAN",
                    "sp|Q1XH10|SKDA1_HUMAN",
                    "sp|Q92826|HXB13_HUMAN",
                    "sp|Q92826|HXB13_HUMAN",
                    "sp|Q92826|HXB13_HUMAN",
                    "sp|O43402|EMC8_HUMAN/tr|M0R1B0|M0R1B0_HUMAN",
                    "sp|O43402|EMC8_HUMAN/tr|M0R1B0|M0R1B0_HUMAN",
                    "sp|O43402|EMC8_HUMAN/tr|M0R1B0|M0R1B0_HUMAN")

  expect_equal(assign_protein_accessions(sequence = seq, fasta_vec = fasta_names),
               expected_res)
})


