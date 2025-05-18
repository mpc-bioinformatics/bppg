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



test_that("test foldChange", {

  # Create fake data
  df <- list()
  df <- c(df, sequence = list(paste0("pep_", 1:10)))
  for (i in 1:2) {
    set.seed(i+2)
    df[[paste0("sample", i)]] <- runif(10, min = 15, max = 25)
    num_na <- sample(1:10, size = sample(1:2, 1))
    df[[paste0("sample", i)]][num_na] <- NA

  }
  df <- as.data.frame(df)

  # Calculate fold changes
  fc_1 <- foldChange(D = df,
                     X = "sample1",
                     Y = "sample2",
                     useNA = FALSE)
  fc_2 <- foldChange(D = df,
                     X = "sample1",
                     Y = "sample2",
                     useNA = TRUE)

  res_1 <- c(1.250, 0.654, 0.952, 0.972, 1.101, 0.837,    NA,    NA, 1.179,    NA)
  res_2 <- c(1.250, 0.654, 0.952, 0.972, 1.101, 0.837,   Inf, 0.000, 1.179,   Inf)

  # Testing
  expect_equal(round(fc_1, digits = 3), res_1)
  expect_equal(round(fc_2, digits = 3), res_2)

})



