

test_that("test proteinElimination", {

    # Create edgelist
    proteins <- rep(paste0("prot_", 1:5), times = c(3,2,3,2,3))
    peptides <- c(paste0("pep_", 1:3), paste0("pep_", 1:2), paste0("pep_", 2:4),
                  paste0("pep_", 4:5), paste0("pep_", 3:5))
    edgelist <- as.matrix(data.frame(protein = proteins, peptide = peptides))

    G <- igraph::graph_from_edgelist(edgelist, directed = FALSE)
    igraph::V(G)$type <- startsWith(igraph::V(G)$name, "prot_")
    igraph::V(G)$pep_ratio[startsWith(igraph::V(G)$name, "prot_")] <- NA
    igraph::V(G)$pep_ratio[startsWith(igraph::V(G)$name, "pep_")] <- c(1.1, 1.2, 1, 0.8, 0.9)


    res <- proteinElimination(G = G)

    ## prepare result object for snapshot (cant deal with random igraph ids)
    res$protein_nodes_list <- names(res$protein_nodes_list)
    res$G_current <- igraph::as_edgelist(res[["G_current"]][[1]])

    expect_snapshot(res)

})










