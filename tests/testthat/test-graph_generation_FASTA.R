
test_that("digestion of a FASTA file", {
  digested_proteins <- readRDS(testthat::test_path("testfiles/digested_proteins_test.rds"))

  file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
  names(fasta) <- limma::strsplit2(names(fasta), "\\|")[,2]
  res <- digest_fasta(fasta)

  expect_equal(res, digested_proteins)
})


test_that("generation of an edgelist", {
  edgelist <- readRDS(testthat::test_path("testfiles/edgelist_test.rds"))

  digested_proteins <- readRDS(testthat::test_path("testfiles/digested_proteins_test.rds"))
  res <- generate_edgelist(digested_proteins)

  expect_equal(res, edgelist)
})


test_that("collapsing of edgelists", {
  edgelist_coll_pept_prot <- readRDS(testthat::test_path("testfiles/edgelist_coll_pept_prot_test.rds"))

  edgelist <- readRDS(testthat::test_path("testfiles/edgelist_test.rds"))
  res <- bppg::collapse_edgelist(edgelist,
                                           collapse_protein_nodes = TRUE,
                                           collapse_peptide_nodes = TRUE)

  expect_equal(res, edgelist_coll_pept_prot)


  edgelist_coll_prot <- readRDS(testthat::test_path("testfiles/edgelist_coll_prot_test.rds"))

  edgelist <- readRDS(testthat::test_path("testfiles/edgelist_test.rds"))
  res2 <- bppg::collapse_edgelist(edgelist,
                                 collapse_protein_nodes = TRUE,
                                 collapse_peptide_nodes = FALSE)

  expect_equal(res2, edgelist_coll_prot)
})



test_that("generation of graphs from edgelist", {
  library(igraph)

  graphs_coll_pept_prot <- readRDS(testthat::test_path("testfiles/graphs_coll_pept_prot_test.rds"))
  graphs_coll_prot <- readRDS(testthat::test_path("testfiles/graphs_coll_prot_test.rds"))

  # with collapsing of peptide and protein nodes
  edgelist_coll_pept_prot <- readRDS(testthat::test_path("testfiles/edgelist_coll_pept_prot_test.rds"))

  res <- bppg::generate_graphs_from_edgelist(edgelist_coll_pept_prot)
  expect_true(bppg::isomorphic_bipartite(res[[1]], graphs_coll_pept_prot[[1]]))
  expect_true(bppg::isomorphic_bipartite(res[[2]], graphs_coll_pept_prot[[2]]))
  expect_true(bppg::isomorphic_bipartite(res[[3]], graphs_coll_pept_prot[[3]]))

  # with collapsing of only protein nodes
  edgelist_coll_prot <- readRDS(test_path("testfiles/edgelist_coll_prot_test.rds"))
  res2 <- bppg::generate_graphs_from_edgelist(edgelist_coll_prot)
  expect_true(bppg::isomorphic_bipartite(res2[[1]], graphs_coll_prot[[1]]))
  expect_true(bppg::isomorphic_bipartite(res2[[2]], graphs_coll_prot[[2]]))
  expect_true(bppg::isomorphic_bipartite(res2[[3]], graphs_coll_prot[[3]]))

})



test_that("subgraph characteristics table", {

  expected <- data.frame(
    graph_ID = c(1L,2L,3L),
    nr_protein_nodes = c(7L,1L,2L),
    nr_peptide_nodes = c(13L,1L,2L),
    nr_unique_peptide_nodes = c(7L,1L,1L),
    nr_shared_peptide_nodes = c(6L,0L,1L),
    nr_edges = c(22L,1L,3L),
    nr_protein_accessions = c(7L,1L,2L),
    nr_peptide_sequences = c(476L, 204L,47L),
    nr_prot_node_only_unique_pep = c(0,1,0L),
    nr_prot_node_unique_and_shared_pep = c(7L,0L,1L),
    nr_prot_node_only_shared_pep = c(0L, 0L, 1L),
    comparison = c(1L,1L, 1L)
  )

  G <- readRDS(test_path("testfiles/graphs_coll_pept_prot_test.rds"))

  res <- bppg::calculate_subgraph_characteristics(S = G, #S2, S3,
                                                  fastalevel = TRUE,
                                                  #comparison = NULL,
                                                  file = NULL)
  expect_equal(res, expected)

})




#
# edgelist <- readRDS(test_path("testfiles/edgelist_test.rds"))
# edgelist$pep_ratio <- NA
# # zum testen random Werte für die Peptid-Ratios, aber ein spezielles Peptid muss immer das gleiche Ratio in jeder Zeile haben
# n <- length(unique(edgelist$peptide))
# random <- data.frame(peptide = unique(edgelist$peptide), pep_ratio = abs(rnorm(n)))
# for (i in 1:nrow(edgelist)) {
#   edgelist$pep_ratio[i] <-random$pep_ratio[which(random$peptide == edgelist$peptide[i])]
#
# }
#
#
#
# edgelist_coll <- bppg::collapse_edgelist_quant(edgelist,
#                                                collapse_protein_nodes = TRUE,
#                                                collapse_peptide_nodes = TRUE)

