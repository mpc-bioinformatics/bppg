test_that("test .errorEquation", {

  Ri <- c(0.5, 1.3)
  Ci <- c(0.3, 0.7)
  M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
  rj <- c(0.6, 1.2)
  e <- bppg:::.errorEquation(Ri, Ci, M, rj)

  expect_snapshot()
})



test_that("test .minimizeSquaredError", {

  M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = FALSE)
  rj <- c(0.6, 1.2)
  S <- list(X = M, fc = rj)
  res <- bppg:::.minimizeSquaredError(S, verbose = FALSE)

  # test fixed_Ci
  res2 <- bppg:::.minimizeSquaredError(S, verbose = FALSE, fixed.Ci = c(0.3, NA))

  expect_snapshot(res)
  expect_snapshot(res2)
})


