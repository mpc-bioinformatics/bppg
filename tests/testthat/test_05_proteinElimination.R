test_that("test proteinElimination", {

    # Create edgelist
    proteins <- rep(paste0("prot_", 1:5), times = c(3,2,3,2,3))
    peptides <- c(paste0("pep_", 1:3), paste0("pep_", 1:2), paste0("pep_", 2:4),
                  paste0("pep_", 4:5), paste0("pep_", 3:5))
    edgelist <- as.matrix(data.frame(protein = proteins, peptide = peptides))

    G <- igraph::graph_from_edgelist(edgelist, directed = FALSE)
    igraph::V(G)$type <- startsWith(igraph::V(G)$name, "prot_")
    pep_logRatio <- rep(NA, igraph::vcount(G))
    pep_logRatio[startsWith(igraph::V(G)$name, "prot_")] <- NA
    pep_logRatio[startsWith(igraph::V(G)$name, "pep_")] <- log2(c(1.1, 1.2, 1, 0.8, 0.9))

    G <- igraph::set_vertex_attr(G, "pep_logRatio", value = pep_logRatio)

    res <- proteinElimination(G = G)

    ## prepare result object for snapshot (cant deal with random igraph ids)
    res$protnodes_list <- names(res$protnodes_list)
    res$res_best$G <- lapply(res$res_best$G, igraph::as_edgelist)

    expect_snapshot(res, variant = Sys.info()[["sysname"]])

})
