test_that("generate graphs from fasta",{
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

  res <- bppg::generateGraphsFromFASTA(fasta = fasta,
                                       outpath = temp_dir)

  expect_true(file.exists(file.path(temp_dir, "edgelist_.txt")))
  expect_true(file.exists(file.path(temp_dir, "subgraphs_collprotpept_.rds")))

  expect_snapshot(igraph::as_edgelist(res[[1]]))
  expect_snapshot(igraph::as_edgelist(res[[2]]))
  expect_snapshot(igraph::as_edgelist(res[[3]]))

})

test_that("test generateGraphsFromQuantData", {
    # skip("test generateGraphsFromQuantData")
    temp_dir <- tempfile(pattern = "test_dir")
    dir.create(temp_dir)
    on.exit(unlink(temp_dir, recursive = TRUE))

    # Load fasta
    file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta", package = "bppg")
    fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

    # Load intensity table
    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    D <- readMqPeptideTable(file)

    # Compute function
    graphs <- bppg::generateGraphsFromQuantData(D = D,
                                            fasta = fasta,
                                            outpath = paste0(temp_dir, "/"))

    # Check results
    expect_true(file.exists(file.path(temp_dir, "edgelist_fasta_.xlsx")))
    expect_true(file.exists(file.path(temp_dir, "aggr_peptides_.xlsx")))
    expect_true(file.exists(file.path(temp_dir, "peptide_ratios_.xlsx")))
    expect_true(file.exists(file.path(temp_dir, "edgelist_filtered_.xlsx")))

    for (i in 1:4) {
        for (j in seq_along(graphs[[i]])) {
            expect_snapshot(igraph::as_edgelist(graphs[[i]][[j]]))
            expect_snapshot(igraph::vertex_attr(graphs[[i]][[j]], "pep_logRatio"))
        }
    }

})
