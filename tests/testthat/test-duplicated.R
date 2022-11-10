test_that("duplicated.dgCMatrix works column wise", {
  M <- Matrix::Matrix(c(0,0,1,0,0,1,1 ,0, 0,0), byrow = FALSE, nrow = 2)
  expect_equal(duplicated.dgCMatrix(M, MARGIN = 2), c(FALSE, FALSE, FALSE, TRUE, TRUE))
})


test_that("duplicated.dgCMatrix works row wise", {
  M <- Matrix::Matrix(c(0,0,1,0,0,1,1 ,0, 0,0), byrow = TRUE, ncol = 2)
  expect_equal(duplicated.dgCMatrix(M, MARGIN = 1), c(FALSE, FALSE, FALSE, TRUE, TRUE))
})

