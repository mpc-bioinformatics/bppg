

test_that("test .generateQuantGraphs", {

  # Create a temporary directory so no permanent files are put on a package users directory
  temp_dir <- tempfile(pattern = "test_dir")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  # Create a ratio table and edgelist
  set.seed(8)
  ratio_table <- data.frame(peptides = paste0("pep_", 1:10),
                            ratio_sample1_sample2 = round(runif(10, min = 0.9, max = 1.1), digits = 3),
                            ratio_sample1_sample3 = round(runif(10, min = 0.9, max = 1.1), digits = 3),
                            ratio_sample2_sample3 = round(runif(10, min = 0.9, max = 1.1), digits = 3))
  for (i in 2:4) {
    ratio_table[sample(1:10, size = 2), i] <- NA # Insert some NAs
  }

  proteins <- rep(paste0("prot_", 1:5), times = c(4,2,3,4,4))
  peptides <- c(paste0("pep_", 1:4), paste0("pep_", 3:4), paste0("pep_", 5:7), paste0("pep_", 7:10), paste0("pep_", 7:10))
  edgelist <- data.frame(protein = proteins, peptide = peptides)

  # Compute function
  graphs <- bppg:::.generateQuantGraphs(peptide_ratios = ratio_table,
                                  id_cols = 1,
                                  fasta_edgelist = edgelist,
                                  outpath = temp_dir,
                                  seq_column = "peptides",
                                  collProtNodes = TRUE,
                                  collPeptNodes = TRUE,
                                  suffix = "")

  # Check result attributes
  expect_true(file.exists(paste0(temp_dir, "edgelist_filtered_.xlsx")))
  expect_equal(unname(lapply(graphs, length)), list(3,2,2))
  expect_equal(names(graphs), c("sample1_sample2", "sample1_sample3", "sample2_sample3"))

  for (i in 1:3) {
    for (j in seq_along(graphs[[i]])) {
      expect_snapshot(igraph::as_edgelist(graphs[[i]][[j]]))
    }
  }

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
  res <- digestFASTA(fasta)
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
  graphs <- generateGraphsFromQuantData(D = data_table,
                                            fasta = fasta,
                                            outpath = temp_dir)


  # Check results
  expect_true(file.exists(paste0(temp_dir, "edgelist_fasta_.xlsx")))
  expect_true(file.exists(paste0(temp_dir, "aggr_peptides_.xlsx")))
  expect_true(file.exists(paste0(temp_dir, "peptide_ratios_.xlsx")))

  expect_equal(unname(lapply(graphs[[1]], length)), list(173, 67, 17))
  expect_equal(unname(lapply(graphs[[2]], length)), list(176, 62, 22))
  expect_equal(unname(lapply(graphs[[3]], length)), list(168, 62, 22))

})


