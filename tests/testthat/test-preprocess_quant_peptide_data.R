test_that("test read_MQ_peptidetable", {

  file <- system.file("extdata", "peptides.txt", package = "bppg")
  D <- read_MQ_peptidetable(path = file, LFQ = TRUE, remove_contaminants = FALSE)

  expect_equal(nrow(D), 7944)
  expect_equal(ncol(D), 28)
  expect_type(D$Sequence, "character")
  expect_type(D[[2]], "integer")
  expect_type(D[[5]], "double")



  D <- read_MQ_peptidetable(path = file, LFQ = TRUE, remove_contaminants = TRUE)

  expect_equal(nrow(D), 7864)
  expect_equal(ncol(D), 28)
})


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
  control_numbers <- c(20.97100, 15.85685, 21.47879)
  for (i in 1:3) {
    expect_equal(which(is.na(D[[i+1]])), indices_na[[i]])
    expect_equal(round(D[[i+1]][[5]], digits = 5), control_numbers[[i]])
  }

  # Test with laxer missing value limit
  D <- aggregate_replicates(D = df,
                            group = factor(rep(1:3, each = 3)),
                            missing.limit = 0.35)
  indices_na <- list(c(6), c(9), integer(0))
  for (i in 1:3) {
    expect_equal(which(is.na(D[[i+1]])), indices_na[[i]])
  }

})


test_that("test foldChange", {


})


test_that("test calculate_peptide_ratios", {


})
