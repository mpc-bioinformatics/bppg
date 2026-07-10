## TODO add test for new contracting functions
test_that("test mapping for graph contraction", {
    file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta", package = "bppg")
    fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
    names(fasta) <- limma::strsplit2(names(fasta), "\\|")[,2]
    res <- digestFASTA(fasta)

    mappingCollProtPept <- bppg:::.getContractMapping(res,
        collProtNodes = TRUE, collPeptNodes = TRUE)
    mappingCollProt <- bppg:::.getContractMapping(res,
        collProtNodes = TRUE, collPeptNodes = FALSE)
    mappingCollPept <- bppg:::.getContractMapping(res,
        collProtNodes = FALSE, collPeptNodes = TRUE)

    expect_snapshot(mappingCollProtPept)
    expect_snapshot(mappingCollProt)
    expect_snapshot(mappingCollPept)
})

test_that("generation of graphs from edgelist", {
    library(igraph)

    file_fasta <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta", package = "bppg")
    file1 <- system.file("extdata", "theoGraphs_collpeptprot.rds", package = "bppg")
    file2 <- system.file("extdata", "theoGraphs_collprot.rds", package = "bppg")

    graphs_coll_pept_prot <- readRDS(file1)
    graphs_coll_prot <- readRDS(file2)

    fasta <- seqinr::read.fasta(file = file_fasta, seqtype = "AA", as.string = TRUE)
    edgelist <- digestFASTA(fasta)

    res <- bppg::generateGraphsFromEdgelist(edgelist, collProtNodes = TRUE, collPeptNodes = TRUE)
    res2 <- bppg::generateGraphsFromEdgelist(edgelist, collProtNodes = TRUE, collPeptNodes = FALSE)

    for (i in 1:4) {
        expect_true(bppg:::.isomorphicBipartite(res[[i]], graphs_coll_pept_prot[[i]]))
        expect_true(bppg:::.isomorphicBipartite(res2[[i]], graphs_coll_prot[[i]]))
    }

    # NOTE: snapshots don't work here because of random graph ids
})

test_that("test generateQuantGraphs", {
    file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta", package = "bppg")
    fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
    edgelist <- digestFASTA(fasta)

    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    group <- factor(rep(1:9, each = 3))
    D <- readMqPeptideTable(path = file, group = group, LFQ = TRUE, remove_contaminants = FALSE)
    D_norm <- bppg::normalizePeptideIntensities(D)
    dAgg <- aggregateReplicates(D_norm)
    exp_peptide_ratios <- calculatePeptideRatios(dAgg)

    graphs <- bppg::generateQuantGraphs(exp_peptide_ratios = exp_peptide_ratios,
                                    fasta_edgelist = edgelist,
                                    collProtNodes = TRUE,
                                    collPeptNodes = TRUE)

    ## imputed case

    imputionMask <- is.na(ratio_table)
    ratio_table[imputionMask] <- min(ratio_table, na.rm = TRUE)

    imputionMask <- as.data.frame(imputionMask)
    rownames(imputionMask) <- rownames(ratio_table)

    impExpData <- SummarizedExperiment::SummarizedExperiment(
        assays = list(logRatios = ratio_table,
            maskImputation = imputionMask),
        rowData = data.frame(peptides = rownames(ratio_table)),
        colData = data.frame(comparison = colnames(ratio_table)),
        metadata = list(imputed = TRUE)
    )
    # Compute function
    graphsImp <- bppg::generateQuantGraphs(exp_peptide_ratios = impExpData,
        fasta_edgelist = edgelist,
        outpath = temp_dir,
        seq_column = "peptides",
        collProtNodes = TRUE,
        collPeptNodes = TRUE,
        suffix = "")

    graphs2 <- bppg::generateQuantGraphs(exp_peptide_ratios = exp_peptide_ratios,
                                        fasta_edgelist = edgelist,
                                        collProtNodes = TRUE,
                                        collPeptNodes = FALSE)

    testfile1 <- system.file("extdata", "quantGraphs.rds", package = "bppg")
    testgraphs <- readRDS(testfile1)
    testfile2 <- system.file("extdata", "quantGraphs_collpept.rds", package = "bppg")
    testgraphs2 <- readRDS(testfile2)

    # Check result attributes
    expect_equal(lengths(graphs), lengths(testgraphs))
    expect_equal(lengths(graphs2), lengths(testgraphs2))

    for (i in 1:3) {
        for (j in seq_along(graphs[[i]])) {
            expect_snapshot(igraph::as_edgelist(graphs[[i]][[j]]))
            expect_snapshot(igraph::as_edgelist(graphs2[[i]][[j]]))
            expect_snapshot(igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio"))
            expect_snapshot(igraph::vertex_attr(graphs2[[i]][[j]], "pep_logRatio"))
        for (j in seq_along(graphsImp[[i]])){
            expect_snapshot(igraph::vertex_attr(graphsImp[[i]][[j]],
                    "pep_ratio_mean"))
            expect_snapshot(igraph::vertex_attr(graphsImp[[i]][[j]],
                    "anyImputed"))
            expect_snapshot(igraph::vertex_attr(graphsImp[[i]][[j]], "imputed"))
        }
    }
})


test_that("test protOrigin", {
    file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta", package = "bppg")
    fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
    names(fasta) <- limma::strsplit2(names(fasta), "\\|")[,2]
    # fake ProtOrigin just for testing
    protOrigin <- as.list(c(rep("yeast", 7), rep("spike_in", 2)))
    edgelist <- digestFASTA(fasta, protOrigin = protOrigin)

    graphs1 <- bppg::generateGraphsFromEdgelist(edgelist)
    graphs2 <- bppg::generateGraphsFromEdgelist(edgelist, collPeptNodes = TRUE)
    graphs3 <- bppg::generateGraphsFromEdgelist(edgelist, collPeptNodes = TRUE,
        collProtNodes = TRUE)

    for (i in seq_along(graphs1)) {
        expect_snapshot(igraph::V(graphs1[[i]])$name[
            igraph::V(graphs1[[i]])$type])
        expect_snapshot(igraph::V(graphs1[[i]])$protOrigin)
        expect_snapshot(igraph::V(graphs2[[i]])$name[
            igraph::V(graphs2[[i]])$type])
        expect_snapshot(igraph::V(graphs2[[i]])$protOrigin)
        expect_snapshot(igraph::V(graphs3[[i]])$name[
            igraph::V(graphs3[[i]])$type])
        expect_snapshot(igraph::V(graphs3[[i]])$protOrigin)
    }

})
