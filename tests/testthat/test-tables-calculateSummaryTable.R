test_that("test summary table", {

    G <- readRDS(test_path("testfiles/graphs_coll_pept_prot_test.rds"))
    tab <- bppg:::.calculateSubgraphCharacteristics(S = G,
                                                    fastalevel = TRUE,
                                                    file = NULL)


    proto_list <- bppg:::.generatePrototypeList(G, sort_by_nr_edges = FALSE)

    res <- bppg:::.calculateSummaryTable(subgraph_char_tab = tab,
                                  isomorph_list = proto_list)
    expect_snapshot(res)

})



