################################################################################
# The tests for the following functions are covered in test-graph_generation_FASTA.R and are skipped here to avoid redundancies.
# bppg::digest_fasta()
# bppg::generate_edgelist()
# bppg::generate_graphs_from_edgelist()
################################################################################



test_that("test aggregate_replicates", {

  # Create fake data (3 samples with 3 runs each)
  df <- list()
  df <- c(df, sequence = list(paste0("pep_", 1:10)))
  for (i in 1:3) {
    for (j in 1:3) {
      set.seed((i+2)^(j+2))
      df[[paste0("sample", i, "_run", j)]] <- runif(10, min = 15, max = 25)
      num_na <- sample(1:10, size = sample(1:3, 1))
      df[[paste0("sample", i, "_run", j)]][num_na] <- NA
    }
  }
  df <- as.data.frame(df)


  D <- aggregate_replicates(D = df,
                            group = factor(rep(1:3, each = 3)))

  # Testing
  expect_equal(ncol(D), 4)
  expect_type(D$sequence, "character")
  expect_type(D[[2]], "double")

  indices_na <- list(c(1,2,6), c(4,9), c(1,3,4,6,7,10))
  control_numbers <- c(20.971, 15.857, 21.479)
  for (i in 1:3) {
    expect_equal(which(is.na(D[[i+1]])), indices_na[[i]])
    expect_equal(round(D[[i+1]][[5]], digits = 3), control_numbers[[i]])
  }

  # Test with laxer missing value limit
  D <- aggregate_replicates(D = df,
                            group = factor(rep(1:3, each = 3)),
                            missing.limit = 0.35)
  indices_na <- list(c(6), c(9), integer(0))
  for (i in 1:3) {
    expect_equal(which(is.na(D[[i+1]])), indices_na[[i]])
  }
  expect_equal(round(D[[2]][[1]], digits = 3), 17.602)

})



test_that("test calculate_peptide_ratios", {

  # Create fake data
  df <- list()
  df <- c(df, sequence = list(paste0("pep_", 1:10)))
  for (i in 1:3) {
    set.seed(i)
    df[[paste0("sample", i)]] <- runif(10, min = 15, max = 25)
    num_na <- sample(1:10, size = sample(0:2, 1))
    df[[paste0("sample", i)]][num_na] <- NA

  }
  df <- as.data.frame(df)


  ratios <- calculate_peptide_ratios(aggr_intensities = df, id_cols = 1)

  res_1 <- c(   NA, 1.176,    NA, 0.693, 1.436, 1.019, 0.666, 1.080, 0.924, 1.313)
  res_2 <- c(   NA, 1.233,    NA, 0.759, 1.235, 0.877,    NA, 0.831, 0.976, 1.364)
  res_3 <- c(0.990, 1.048, 0.909, 1.096, 0.860, 0.861,    NA, 0.769, 1.056, 1.040)

  expect_equal(round(x = ratios[[2]], digits = 3), res_1)
  expect_equal(round(x = ratios[[3]], digits = 3), res_2)
  expect_equal(round(x = ratios[[4]], digits = 3), res_3)


})



test_that("test collapse_edgelist_quant", {

  # Create edgelist (proteins, peptides and pep_ratios and compute the collapsing
  set.seed(4)
  peptides <- c(rep(paste0("pep_", 1:2), each = 2), rep(paste0("pep_", 3:3), each = 4), rep(paste0("pep_", 4:5), each = 3))
  ratios <- round(runif(5, min = 0.9, max = 1.1), digits = 2)
  pep_ratios <- unlist(mapply(rep, ratios, each = c(2, 2, 4, 3, 3)))
  proteins <- c(paste0("prot_", 1:2), paste0("prot_", 1:2), paste0("prot_", 2:5), paste0("prot_", 3:5), paste0("prot_", 3:5))

  edgelist <- data.frame(protein = proteins, peptide = peptides, pep_ratio = pep_ratios)

  collapsed_edgelist <- collapse_edgelist_quant(edgelist = edgelist, collapse_protein_nodes = TRUE, collapse_peptide_nodes = TRUE)


  # The expected result
  res_proteins <- c("prot_1", rep("prot_2", each = 2), rep("prot_3;prot_4;prot_5", each = 2))
  res_peptides <- c(rep("pep_1;pep_2", each = 2), rep("pep_3", each = 2), "pep_4;pep_5")
  res_pep_ratios <- c(rep("0.9;1.02", each = 2), rep("0.96", each = 2), "0.96;1.06")
  res_collapsed_edgelist <- data.frame(protein = res_proteins, peptide = res_peptides, pep_ratio = res_pep_ratios)


  # Check expected vs. actual result
  expect_equal(collapsed_edgelist[["protein"]],
               res_collapsed_edgelist[["protein"]])
  expect_equal(collapsed_edgelist[["peptide"]],
               res_collapsed_edgelist[["peptide"]])
  expect_equal(collapsed_edgelist[["pep_ratio"]],
               res_collapsed_edgelist[["pep_ratio"]])

})




























