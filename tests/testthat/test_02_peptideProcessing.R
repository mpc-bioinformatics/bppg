test_that("test aggregateReplicates", {

    # Create test data (3 samples with 3 runs each)
    df <- list()
    df <- c(df, sequence = list(paste0("pep_", 1:10)))
    for (i in 1:3) {
        for (j in 1:3) {
            set.seed((i + 2)^(j + 2))
            df[[paste0("sample", i, "_run", j)]] <- runif(10, min = 15, max = 25)
            num_na <- sample(1:10, size = sample(1:3, 1))
            df[[paste0("sample", i, "_run", j)]][num_na] <- NA
        }
    }
    df <- as.data.frame(df)
    rownames(df) <- df$sequence
    df$sequence <- NULL
    D <- SummarizedExperiment::SummarizedExperiment(
            assays = list(intensities = df),
            colData = data.frame(sample = colnames(df)),
            rowData = data.frame(Sequence = rownames(df)))

    D1 <- bppg::aggregateReplicates(D = D, group = factor(rep(1:3, each = 3)))

    D2 <- bppg::aggregateReplicates(D = D,
                                    group = factor(rep(1:3, each = 3)),
                                    missing.limit = 0.35,
                                    method = "median")
    expect_snapshot(SummarizedExperiment::assays(D1)$intensities)
    expect_snapshot(SummarizedExperiment::assays(D2)$intensities)
    expect_snapshot(D1)
    expect_snapshot(D2)
})



test_that("test calculatePeptideRatios", {

    # Create test data
    df <- list()
    df <- c(df, sequence = list(paste0("pep_", 1:10)))
    for (i in 1:3) {
        set.seed(i)
        df[[paste0("sample", i)]] <- runif(10, min = 15, max = 25)
        num_na <- sample(1:10, size = sample(0:2, 1))
        df[[paste0("sample", i)]][num_na] <- NA
    }
    df <- as.data.frame(df)
    rownames(df) <- df$sequence
    df$sequence <- NULL
    D <- SummarizedExperiment::SummarizedExperiment(
        assays = list(intensities = df),
        colData = data.frame(group = colnames(df)),
        rowData = data.frame(Sequence = rownames(df)))

    D1 <- bppg::calculatePeptideRatios(D = D)

    expect_snapshot(SummarizedExperiment::assays(D1)$logRatios)
    expect_snapshot(D1)
})
