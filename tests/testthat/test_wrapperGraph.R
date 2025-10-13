test_that("generate graphs from fasta",{
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))


  file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

  res <- bppg::generateGraphsFromFASTA(fasta = fasta,
                                       save_intermediate = TRUE,
                                       result_path = paste0(temp_dir, "\\"))

  expect_snapshot(igraph::as_edgelist(res[[1]]))
  expect_snapshot(igraph::as_edgelist(res[[2]]))
  expect_snapshot(igraph::as_edgelist(res[[3]]))

  expect_true(file.exists(file.path(temp_dir, "edgelist_.txt")))
  expect_true(file.exists(file.path(temp_dir, "edgelist_collprotpept_.txt")))
  expect_true(file.exists(file.path(temp_dir, "subgraphs_collprotpept_.rds")))

})

test_that("test generateGraphsFromQuantData", {

  # Create a temporary directory so no permanent files are put on a package users directory
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  # Load fasta
  file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
  fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

  # Create intensity table
  set.seed(4)
  res <- bppg::digestFASTA(fasta)
  peptides <- c()
  for (i in 1:10) {
    peptides <- c(peptides, res[[i]][sample(1:length(res[[i]]), size = round(length(res[[i]])*0.75))])
  }
  peptides <- unique(peptides)
  data_table <- data.frame(Sequence = peptides,
                           sample1_run1 = round(rnorm(598, mean = 20), digits = 4),
                           sample1_run2 = round(rnorm(598, mean = 20), digits = 4),
                           sample2_run1 = round(rnorm(598, mean = 20), digits = 4),
                           sample2_run2 = round(rnorm(598, mean = 20), digits = 4),
                           sample3_run1 = round(rnorm(598, mean = 20), digits = 4),
                           sample3_run2 = round(rnorm(598, mean = 20), digits = 4))
  for (i in 2:7) {
    data_table[sample(1:598, size = 120), i] <- NA # Insert some NAs
  }

  # Compute function
  graphs <- bppg::generateGraphsFromQuantData(D = data_table,
                                            fasta = fasta,
                                            outpath = paste0(temp_dir, "/"))

  # Check results
  expect_true(file.exists(file.path(temp_dir, "edgelist_fasta_.xlsx")))
  expect_true(file.exists(file.path(temp_dir, "aggr_peptides_.xlsx")))
  expect_true(file.exists(file.path(temp_dir, "peptide_ratios_.xlsx")))
  expect_true(file.exists(file.path(temp_dir, "edgelist_filtered_.xlsx")))

  for (i in 1:3) {
    for (j in seq_along(graphs[[i]])) {
      expect_snapshot(igraph::as_edgelist(graphs[[i]][[j]]))
      expect_snapshot(igraph::vertex_attr(graphs[[1]][[1]], "pep_ratio"))
    }
  }


})
