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

################################################################################
# iterate over Ci

test_that("test iterateOverCi", {

    M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = FALSE)
    rj <- c(0.6, 1.2)
    S <- list(X = M, fc = rj)
    res <- iterateOverCi(S, grid.size = 10)
    expect_snapshot(res)
})

test_that("test iterateOverCi with extended grid", {

    M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = FALSE)
    rj <- c(0.6, 1.2)
    S <- list(X = M, fc = rj)

    res2 <- iterateOverCi(S, grid.size = 10,
                          extend_grid_at_borders = TRUE)

    expect_snapshot(res2)
})

################################################################################
# automated analysis

test_that("test automated analysis", {

    M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = FALSE)
    rj <- c(0.6, 1.2)
    S <- list(X = M, fc = rj)

    res <- iterateOverCi(S, grid.size = 10)

    res3 <- automatedAnalysisIteratedCi(S, res)
    expect_snapshot(res3)

})


test_that("test automated analysis with using results from other proteins", {

    M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = FALSE)
    rj <- c(0.6, 1.2)
    S <- list(X = M, fc = rj)

    res <- iterateOverCi(S, grid.size = 10)

    res4 <- automatedAnalysisIteratedCi(S, res,
                                        use_results_from_other_proteins = TRUE)
    expect_snapshot(res4)

})


test_that("test automated analysis with only one protein node", {

    M <- matrix(c(1,1,1), nrow = 3, byrow = FALSE)
    rj <- c(0.6, 1.2, 0.9)
    S <- list(X = M, fc = rj)

    res5 <- automatedAnalysisIteratedCi(S, res = NULL)
    expect_snapshot(res5)
})


test_that("test automated analysis with constant error", {

    M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = FALSE)
    rj <- c(0.6, 1.2)
    S <- list(X = M, fc = rj)

    set.seed(63657)
    res <- data.frame(protein = rep(1:2, each = 9),
                      grid = rep(seq(0.1,0.9, by = 0.1), 2),
                      R1 = rep(0.6, 18),
                      R2 = runif(18, 0,1),
                      C1 = c(seq(0.1,0.9, by = 0.1), rev(seq(0.1, 0.9, by = 0.1))),
                      C2 = c(rev(seq(0.1,0.9, by = 0.1)), seq(0.1, 0.9, by = 0.1)),
                      error = rep(c(1e-5, 1e-10), each = 9))


    res6 <- automatedAnalysisIteratedCi(S, res = res, verbose = TRUE)
    expect_snapshot(res6)
})
