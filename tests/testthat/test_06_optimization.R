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
    res2 <- bppg:::.minimizeSquaredError(G_N, verbose = FALSE, fixedCi = c(0.3, NA))

    testfile_path <- file.path(testthat::test_path(), "testfiles")
    #saveRDS(res, file.path(testfile_path, "test_minimizeSquaredError_file1.rds"))
    #saveRDS(res2, file.path(testfile_path, "test_minimizeSquaredError_file2.rds"))
    res_snap <- readRDS(file.path(testfile_path, "test_minimizeSquaredError_file1.rds"))
    res2_snap <- readRDS(file.path(testfile_path, "test_minimizeSquaredError_file2.rds"))

    expect_equal(res$Ri, res_snap$Ri, tolerance = 1e-05)
    expect_equal(res$Ci, res_snap$Ci, tolerance = 1e-05)
    # expect_equal(str(res$Ci), str(res_snap$Ci))

    expect_equal(res2$Ri, res2_snap$Ri, tolerance = 1e-05)
    expect_equal(res2$Ci, res2_snap$Ci, tolerance = 1e-05)
    # expect_equal(str(res2$Ci), str(res2_snap$Ci))

})

################################################################################
# iterate over Ci

test_that("test iterateOverCi", {

  testfile_path <- file.path(testthat::test_path(), "testfiles")
  graphs <- readRDS(file.path(testfile_path, "quantGraphsForTesting.rds"))
  G_I <- graphs[[1]][[1]] # I
  G_N <- graphs[[3]][[1]] # N
  G_M <- graphs[[3]][[2]] # M

  res_I <- iterateOverCi(G_I, gridSize = 10)
  #saveRDS(res_I, file.path(testfile_path, "test_iterateOverCi_file_I1.rds"))
  res_I_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_I1.rds"))
  expect_equal(res_I, res_I_snap, tolerance = 1e-05)

  res_N <- iterateOverCi(G_N, gridSize = 10)
  #saveRDS(res_N, file.path(testfile_path, "test_iterateOverCi_file_N1.rds"))
  res_N_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_N1.rds"))
  expect_equal(res_N, res_N_snap, tolerance = 1e-05)

  res_M <- iterateOverCi(G_M, gridSize = 10)
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

  res_I <- iterateOverCi(G_I, gridSize = 10, extend_grid_at_borders = TRUE)
  #saveRDS(res_I, file.path(testfile_path, "test_iterateOverCi_file_I2.rds"))
  res_I_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_I2.rds"))
  expect_equal(res_I, res_I_snap, tolerance = 1e-05)

  res_N <- iterateOverCi(G_N, gridSize = 10, extend_grid_at_borders = TRUE)
  #saveRDS(res_N, file.path(testfile_path, "test_iterateOverCi_file_N2.rds"))
  res_N_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_N2.rds"))
  expect_equal(res_N, res_N_snap, tolerance = 1e-05)

  res_M <- iterateOverCi(G_M, gridSize = 10, extend_grid_at_borders = TRUE)
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

  res_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_I1.rds"))
  res <- automatedAnalysisIteratedCi(G_I, res_snap) # case 3
  expect_snapshot(SummarizedExperiment::assay(res))

  res_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_N1.rds"))
  res2 <- automatedAnalysisIteratedCi(G_N, res_snap) # case 1 and 2
  expect_snapshot(SummarizedExperiment::assay(res2))

  res_snap <- readRDS(file.path(testfile_path, "test_iterateOverCi_file_M1.rds"))
  res3 <- automatedAnalysisIteratedCi(G_M, res_snap) # case 3
  expect_snapshot(SummarizedExperiment::assay(res3))

  res4 <- automatedAnalysisIteratedCi(G_M, res_snap, use_results_from_other_proteins = TRUE) # case 4
  expect_snapshot(SummarizedExperiment::assay(res4))
  # comparison graphID proteinNr  error_min    RiLog RiLog_min RiLog_max Ci Ci_min Ci_max case
  # 1         NA      NA         1 0.00663525 1.026032        NA        NA NA    0.1    0.1    4
  # 2         NA      NA         2 0.00663525 0.983554        NA        NA NA    0.9    0.9    4
  # TODO: it switches to case 4 because for both proteins an almost identical result is achieved

})


