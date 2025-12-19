test_that("test .errorEquation", {

  # constructed example
  RiLog <- log2(c(0.5, 1.3))
  Ci <- c(0.3, 0.7)
  M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
  rjLog <- log2(c(0.6, 1.2))
  e <- bppg:::.errorEquation(RiLog = RiLog, Ci = Ci, M = M, rjLog = rjLog)
  expect_snapshot(e)

  # real example
  testfile_path <- file.path(testthat::test_path(), "testfiles")
  graphs <- readRDS(file.path(testfile_path, "quantGraphsForTesting.rds"))

})

test_that("test .minimizeSquaredError", {

  testfile_path <- file.path(testthat::test_path(), "testfiles")
  graphs <- readRDS(file.path(testfile_path, "quantGraphsForTesting.rds"))
  G_N <- graphs[[3]][[1]] # N

  res <- bppg:::.minimizeSquaredError(G_N, verbose = FALSE)

  # test fixed_Ci
  res2 <- bppg:::.minimizeSquaredError(G_N, verbose = FALSE, fixed.Ci = c(0.3, NA))

  testfile_path <- file.path(testthat::test_path(), "testfiles")
  #saveRDS(res, file.path(testfile_path, "test_minimizeSquaredError_file1.rds"))
  #saveRDS(res2, file.path(testfile_path, "test_minimizeSquaredError_file2.rds"))
  res_snap <- readRDS(file.path(testfile_path, "test_minimizeSquaredError_file1.rds"))
  res2_snap <- readRDS(file.path(testfile_path, "test_minimizeSquaredError_file2.rds"))

  expect_equal(res$Ri, res_snap$Ri, tolerance = 1e-05)
  expect_equal(res$Ci, res_snap$Ci, tolerance = 1e-05)
  expect_equal(str(res$Ci), str(res_snap$Ci))

  expect_equal(res2$Ri, res2_snap$Ri, tolerance = 1e-05)
  expect_equal(res2$Ci, res2_snap$Ci, tolerance = 1e-05)
  expect_equal(str(res2$Ci), str(res2_snap$Ci))

})

################################################################################
# iterate over Ci

test_that("test iterateOverCi", {

  testfile_path <- file.path(testthat::test_path(), "testfiles")
  graphs <- readRDS(file.path(testfile_path, "quantGraphsForTesting.rds"))
  G_I <- graphs[[1]][[1]] # I
  G_N <- graphs[[3]][[1]] # N
  G_M <- graphs[[3]][[2]] # M

  res_I <- iterateOverCi(G_I, grid.size = 10)
  #saveRDS(res_I, file.path(testfile_path, "test_iterateOverCi_file_I1.rds"))
  res_I_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_I1.rds"))
  expect_equal(res_I, res_I_snap, tolerance = 1e-05)

  res_N <- iterateOverCi(G_N, grid.size = 10)
  #saveRDS(res_N, file.path(testfile_path, "test_iterateOverCi_file_N1.rds"))
  res_N_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_N1.rds"))
  expect_equal(res_N, res_N_snap, tolerance = 1e-05)

  res_M <- iterateOverCi(G_M, grid.size = 10)
  #saveRDS(res_M, file.path(testfile_path, "test_iterateOverCi_file_M1.rds"))
  res_M_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_M1.rds"))
  expect_equal(res_M, res_M_snap, tolerance = 1e-05)


})

# extend_grid_at_borders = TRUE



test_that("test iterateOverCi with extended grid", {

  testfile_path <- file.path(testthat::test_path(), "testfiles")
  graphs <- readRDS(file.path(testfile_path, "quantGraphsForTesting.rds"))
  G_I <- graphs[[1]][[1]] # I
  G_N <- graphs[[3]][[1]] # N
  G_M <- graphs[[3]][[2]] # M

  res_I <- iterateOverCi(G_I, grid.size = 10, extend_grid_at_borders = TRUE)
  #saveRDS(res_I, file.path(testfile_path, "test_iterateOverCi_file_I2.rds"))
  res_I_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_I2.rds"))
  expect_equal(res_I, res_I_snap, tolerance = 1e-05)

  res_N <- iterateOverCi(G_N, grid.size = 10, extend_grid_at_borders = TRUE)
  #saveRDS(res_N, file.path(testfile_path, "test_iterateOverCi_file_N2.rds"))
  res_N_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_N2.rds"))
  expect_equal(res_N, res_N_snap, tolerance = 1e-05)

  res_M <- iterateOverCi(G_M, grid.size = 10, extend_grid_at_borders = TRUE)
  #saveRDS(res_M, file.path(testfile_path, "test_iterateOverCi_file_M2.rds"))
  res_M_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_M2.rds"))
  expect_equal(res_M, res_M_snap, tolerance = 1e-05)
})



################################################################################
# automated analysis

test_that("test automated analysis", {
  testfile_path <- file.path(testthat::test_path(), "testfiles")
  graphs <- readRDS(file.path(testfile_path, "quantGraphsForTesting.rds"))
  G_I <- graphs[[1]][[1]] # I
  G_N <- graphs[[3]][[1]] # N
  G_M <- graphs[[3]][[2]] # M

  testfile_path <- file.path(testthat::test_path(), "testfiles")
  res_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_I1.rds"))
  res <- automatedAnalysisIteratedCi(G_I, res_snap)
  expect_snapshot(res)

  testfile_path <- file.path(testthat::test_path(), "testfiles")
  res_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_N1.rds"))
  res2 <- automatedAnalysisIteratedCi(G_N, res_snap)
  expect_snapshot(res2)

  testfile_path <- file.path(testthat::test_path(), "testfiles")
  res_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_M1.rds"))
  res3 <- automatedAnalysisIteratedCi(G_M, res_snap)
  expect_snapshot(res3)
})


# ### TODO: anderer usecase wo das wirklich nen Unterschied macht?
# test_that("test automated analysis with using results from other proteins", {
#     M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = FALSE)
#     rj <- c(0.6, 1.2)
#     S <- list(X = M, fc = rj)
#
#     testfile_path <- file.path(testthat::test_path(), "testfiles")
#     res_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file1.rds"))
#
#     res4 <- automatedAnalysisIteratedCi(S, res_snap,
#                                         use_results_from_other_proteins = TRUE)
#     expect_snapshot(res4)
#
# })


# test_that("test automated analysis with only one protein node", {
#
#     M <- matrix(c(1,1,1), nrow = 3, byrow = FALSE)
#     rj <- c(0.6, 1.2, 0.9)
#     S <- list(X = M, fc = rj)
#
#     res5 <- automatedAnalysisIteratedCi(S, res = NULL)
#     expect_snapshot(res5)
# })
#
#
# test_that("test automated analysis with constant error", {
#
#     M <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = FALSE)
#     rj <- c(0.6, 1.2)
#     S <- list(X = M, fc = rj)
#
#     set.seed(63657)
#     res <- data.frame(protein = rep(1:2, each = 9),
#                       grid = rep(seq(0.1,0.9, by = 0.1), 2),
#                       R1 = rep(0.6, 18),
#                       R2 = runif(18, 0,1),
#                       C1 = c(seq(0.1,0.9, by = 0.1), rev(seq(0.1, 0.9, by = 0.1))),
#                       C2 = c(rev(seq(0.1,0.9, by = 0.1)), seq(0.1, 0.9, by = 0.1)),
#                       error = rep(c(1e-5, 1e-10), each = 9))
#
#
#     res6 <- automatedAnalysisIteratedCi(S, res = res, verbose = TRUE)
#     expect_snapshot(res6)
# })
