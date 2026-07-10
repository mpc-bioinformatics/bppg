test_that("test .addUniquenessAttributes", {

  file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
  graphs <- readRDS(file)
  G <- graphs[[1]][[4]]

  G_new <- bppg:::.addUniquenessAttributes(G)

  expect_snapshot(igraph::vertex_attr(G_new))

})

test_that("geometric mean", {

  x <- c(1,5,4.5,2)
  n <- length(x)

  res <- exp(mean(log(x)))
  res_prod <- prod(x)^(1/n)

  expect_equal(res, bppg:::.geomMean(x, useprod = FALSE))
  expect_equal(res_prod, bppg:::.geomMean(x, useprod = TRUE))
})
