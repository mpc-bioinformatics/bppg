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
