



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

    res <- bppg:::.calculateSubgraphCharacteristics(S = G,
                                                    fastalevel = TRUE,
                                                    file = NULL)
    expect_equal(res, expected)

})



