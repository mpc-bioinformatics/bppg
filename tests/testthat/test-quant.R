test_that("test equation", {

  Ri <- c(0.5, 1.3)
  Ci <- c(0.3, 0.7)
  M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
  rj <- c(0.6, 1.2)
  e <- bppg::equation(Ri, Ci, M, rj)

  e[[2]] <- round(e[[2]], digits = 4)
  e[[3]] <- round(e[[3]], digits = 4)

  expected_res <- list(res_Mat = cbind(c(0.50, 0.15), c(0.00, 0.91)),
                       res_equ = c(0.1823, 0.1241),
                       res_squ_err = 0.0486,
                       W = cbind(c(1.0, 0.3), c(0.0, 0.7)))

  expect_equal(e, expected_res)
})



test_that("test minimize_squared_error", {

  M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
  rj <- c(0.6, 1.2)
  S <- list(X = M, fc = rj)
  e <- minimize_squared_error(S)


  e[[1]] <- round(e[[1]], digits = 4)
  e[[2]] <- round(e[[2]], digits = 4)
  expected_res <- list(Ri = c(0.6000, 1.4716),
                       Ci = c(0.3116, 0.6884),
                       outer.iter = 2,
                       convergence = 0)
  expect_equal(e[c(1,2,5,6)], expected_res)


  for (i in 1:4) {
    e[[3]][[i]] <- round(e[[3]][[i]], digits = 4)
  }
  expected_res <- list(res_Mat = cbind(c(0.6, 0.187), c(0.0, 1.013)),
                       res_equ = c(0, 0),
                       res_squ_err = 0,
                       W = cbind(c(1, 0.3116), c(0, 0.6884)))
  expect_equal(e[[3]], expected_res)


  e[[4]] <- round(e[[4]], digits = 4)
  expect_equal(e[[4]], data.frame(iter = c(0, 1),
                                  squ_err = c(0.0828, 0),
                                  R1 = c(-0.737, -0.737),
                                  R2 = c(0.2630, 0.5574),
                                  C1 = c(0.5, 0.3116),
                                  C2 = c(0.5, 0.6884)))

})



test_that("test iterate_over_Ci", {

})



test_that("test automated_analysis_iterated_Ci", {

})




test_that("test protein_elimination", {

})


