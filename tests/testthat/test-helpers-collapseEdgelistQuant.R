test_that("test .collapseEdgelistQuant", {

    # Create edgelist (proteins, peptides and pep_ratios and compute the collapsing
    set.seed(4)
    peptides <- c(rep(paste0("pep_", 1:2), each = 2), rep(paste0("pep_", 3:3), each = 4), rep(paste0("pep_", 4:5), each = 3))
    ratios <- round(runif(5, min = 0.9, max = 1.1), digits = 2)
    pep_ratios <- unlist(mapply(rep, ratios, each = c(2, 2, 4, 3, 3)))
    proteins <- c(paste0("prot_", 1:2), paste0("prot_", 1:2), paste0("prot_", 2:5), paste0("prot_", 3:5), paste0("prot_", 3:5))

    edgelist <- data.frame(protein = proteins, peptide = peptides, pep_ratio = pep_ratios)

    collapsed_edgelist <- bppg:::.collapseEdgelistQuant(edgelist = edgelist, collProtNodes = TRUE, collPeptNodes = TRUE)

    collapsed_edgelist2 <- bppg:::.collapseEdgelistQuant(edgelist = edgelist, collProtNodes = TRUE, collPeptNodes = FALSE)
    expect_snapshot(collapsed_edgelist)
    expect_snapshot(collapsed_edgelist2)

})
