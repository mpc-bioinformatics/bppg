test_that("test .errorEquation", {

    RiLog <- log2(c(0.5, 1.3))
    Ci <- c(0.3, 0.7)

    file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
    graphs <- readRDS(file)

    G <- graphs[[1]][[2]]
    M <- as_biadjacency_matrix(G)
    rjLog <- na.omit(V(G)$pep_logRatio)

    e <- bppg:::.errorEquation(RiLog = RiLog, Ci = Ci, M = M, rjLog = rjLog)
    expect_snapshot(e)

})

test_that("test .minimizeSquaredError", {

    file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
    graphs <- readRDS(file)
    G_N <- graphs[[1]][[2]] # N-shaped graph

    res <- bppg:::.minimizeSquaredError(G_N, verbose = FALSE)

    # test fixed_Ci
    res2 <- bppg:::.minimizeSquaredError(G_N, verbose = FALSE, fixedCi = c(0.3, NA))

    res$RES$res_equ<- round(res$RES$res_equ, 5)
    res2$RES$res_equ<- round(res2$RES$res_equ, 5)

    expect_snapshot(res)
    expect_snapshot(res2)

})

################################################################################
# iterate over Ci

test_that("test iterateOverCi and automated analysis", {

  file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
  graphs <- readRDS(file)
  G_I <- graphs[[1]][[1]] # I
  G_N <- graphs[[2]][[2]] # N
  G_M <- graphs[[1]][[3]] # M

   igraph::V(G_I)$protOrigin[igraph::V(G_I)$type] <- c("test_origin")
  res_I <- iterateOverCi(G_I, gridSize = 10)

  res_I2 <- automatedAnalysisIteratedCi(G_I, res_I,
                                     use_results_from_other_proteins = FALSE)

  res_N <- iterateOverCi(G_N, gridSize = 10)
  res_N2 <- automatedAnalysisIteratedCi(G_N, res_N,
                                        use_results_from_other_proteins = FALSE)

  res_M <- iterateOverCi(G_M, gridSize = 10)
  res_M2 <- automatedAnalysisIteratedCi(G_M, res_M,
                                        use_results_from_other_proteins = FALSE)


  expect_snapshot(round(res_I, 4))
  expect_snapshot(res_I2)
  expect_snapshot(SummarizedExperiment::rowData(res_I))
  expect_snapshot(round(res_N, 4))
  expect_snapshot(res_N2)
  expect_snapshot(round(res_M, 4))
  expect_snapshot(res_M2)

  ## test case when an error term is NaN (which may happen during the optimization,
  ## if a Ci is estimated as 0)
  res_M$error[c(1,5)] <- NaN
  res_M3 <- automatedAnalysisIteratedCi(G_M, res_M)
  expect_snapshot(SummarizedExperiment::assay(res_M3))

})


test_that("test iterateOverCi with extended grid", {

  file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
  graphs <- readRDS(file)
  G_I <- graphs[[1]][[1]] # I
  G_N <- graphs[[2]][[2]] # N
  G_M <- graphs[[1]][[3]] # M

  res_I <- iterateOverCi(G_I, gridSize = 10, extend_grid_at_borders = TRUE)
  res_N <- iterateOverCi(G_N, gridSize = 10, extend_grid_at_borders = TRUE)
  res_M <- iterateOverCi(G_M, gridSize = 10, extend_grid_at_borders = TRUE)

  expect_snapshot(round(res_I, 4))
  expect_snapshot(round(res_N, 4))
  expect_snapshot(round(res_M, 4))
})


