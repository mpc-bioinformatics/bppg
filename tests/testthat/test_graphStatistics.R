test_that("test .calculateProteinNodeInfo", {

  file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
  graphs <- readRDS(file)

  ProtInfo <- bppg:::.calculateProteinNodeInfo(G = graphs, verbose = FALSE)

  expect_snapshot(ProtInfo)
})



test_that("subgraph characteristics table", {

  file <- system.file("extdata", "theoGraphs.rds", package = "bppg")
  graphs <- readRDS(file)

  res <- bppg:::.calculateSubgraphCharacteristics(S = graphs,
                                                  fastalevel = TRUE)


  file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
  graphs <- readRDS(file)
  res2 <- bppg:::.calculateSubgraphCharacteristics(S = graphs,
                                                   fastalevel = FALSE)

  expect_snapshot(res)
  expect_snapshot(res2)

})
