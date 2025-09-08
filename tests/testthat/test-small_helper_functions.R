test_that("geometric mean", {

  x <- c(1,5,4.5,2)
  n <- length(x)

  res <- exp(mean(log(x)))
  res_prod <- prod(x)^(1/n)

  expect_equal(res, bppg::.geomMean(x, useprod = FALSE))
  expect_equal(res_prod, bppg::.geomMean(x, useprod = TRUE))
})
