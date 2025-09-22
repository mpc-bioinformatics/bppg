

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


# TODO: prototype argument?
