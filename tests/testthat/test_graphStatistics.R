test_that("test .calculateProteinNodeInfo", {

  file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

  g1 <- generateGraphsFromFASTA(fasta = fasta,
                                   collProtNodes = FALSE,
                                   collPeptNodes = TRUE)
  g2 <- generateGraphsFromFASTA(fasta = fasta,
                                   collProtNodes = TRUE,
                                   collPeptNodes = FALSE)
  g3 <- generateGraphsFromFASTA(fasta = fasta,
                                   collProtNodes = FALSE,
                                   collPeptNodes = FALSE)

  graphs <- list(s1_s2 = g1, s1_s3 = g2, s2_s3 = g3)

  D <- bppg:::.calculateProteinNodeInfo(graphs)

  expect_equal(nrow(D), 21)
  expect_equal(ncol(D), 7)
  expect_equal(unlist(D[10, ]), c(accessions = "sp|Q03483|YD22A_YEAST",
                                  comparison = "s1_s3",
                                  graphID = "2",
                                  ind_within_graph = "2",
                                  nr_peptides = "88",
                                  nr_unique_peptides = "30",
                                  nr_shared_peptides = "58"))
  expect_equal(unlist(D[19, ]), c(accessions = "sp|Q07793|YD13B_YEAST",
                                  comparison = "s2_s3",
                                  graphID = "1",
                                  ind_within_graph = "5",
                                  nr_peptides = "380",
                                  nr_unique_peptides = "254",
                                  nr_shared_peptides = "126"))
})

test_that("subgraph characteristics table", {
    G <- readRDS(test_path("testfiles/graphs_coll_pept_prot_test.rds"))
    res <- bppg:::.calculateSubgraphCharacteristics(S = G,
                                                    fastalevel = TRUE,
                                                    file = NULL)
    expect_snapshot(res)

    # simulate quant data (list of comparisons, each one is a list of graphs)
    G_quant <- list(comp1 = G, comp2 = G)
    res2 <- bppg:::.calculateSubgraphCharacteristics(S = G_quant,
                                                    fastalevel = FALSE,
                                                    file = NULL)
    expect_snapshot(res2)

})