test_that("generation of a prototype list", {
    library(igraph)

    # N shape graph
    M <- matrix(c(1,1,0,1), nrow = 2, byrow = TRUE)
    G <- igraph::graph_from_biadjacency_matrix(M)

    # N shape graph type 2 (different order of nodes)
    M2 <- M[c(2,1), c(2,1)]
    G2 <- igraph::graph_from_biadjacency_matrix(M2)

    # W shape graph
    M3 <- matrix(c(1,1,0,0,1,1), nrow = 2, byrow = TRUE)
    G3 <- igraph::graph_from_biadjacency_matrix(M3)

    # W shape graph type 2 (different order of nodes)
    M4 <- M3[, c(2,1,3)]
    G4 <- igraph::graph_from_biadjacency_matrix(M4)

    # M shape graph
    G5 <- G3
    igraph::V(G5)$type <- !igraph::V(G5)$type

    G_list <- list(G3, G3, G5, G, G, G2)
    # just to test the case when the last element will form its own class:
    G_list2 <- list(G3, G3, G, G, G2, G5)

    proto_list <- bppg:::.generatePrototypeList(G_list, sort_by_nr_edges = TRUE)
    proto_list2 <- bppg:::.generatePrototypeList(G_list, sort_by_nr_edges = FALSE)
    proto_list3 <- bppg:::.generatePrototypeList(G_list2, sort_by_nr_edges = FALSE)

    expect_type(proto_list, "list")
    expect_type(proto_list2, "list")

    expect_length(proto_list, 2)
    expect_length(proto_list2, 2)

    expect_equal(names(proto_list), c("graphs", "counter"))

    expect_length(proto_list$graphs, 3)
    expect_length(proto_list2$graphs, 3)
    expect_length(proto_list3$graphs, 3)

    expect_equal(proto_list$counter, c(3,2,1))
    expect_equal(proto_list2$counter, c(2,1,3))
    expect_equal(proto_list3$counter, c(2,3,1))

})


