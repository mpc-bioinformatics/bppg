

test_that("read MaxQuant Output table", {
    expect_snapshot(bppg::readMqPeptideTable(test_path("testfiles/peptides.txt")))

    expect_snapshot(bppg::readMqPeptideTable(test_path("testfiles/peptides.txt"),
                                             LFQ = TRUE,
                                             further_columns_to_keep = c("Proteins", "Score")))

})


test_that("test aggregateReplicates", {

    # Create test data (3 samples with 3 runs each)
    # TODO: possibility to make generation of test data more simple?
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

    expect_snapshot(bppg::aggregateReplicates(D = df,
                                              group = factor(rep(1:3, each = 3)))
    )

    expect_snapshot(bppg::aggregateReplicates(D = df,
                                              group = factor(rep(1:3, each = 3)),
                                              missing.limit = 0.35,
                                              method = "median")
    )

})



test_that("test calculatePeptideRatios", {

    # Create test data
    # TODO: possibility to make generation of test data more simple?
    df <- list()
    df <- c(df, sequence = list(paste0("pep_", 1:10)))
    for (i in 1:3) {
        set.seed(i)
        df[[paste0("sample", i)]] <- runif(10, min = 15, max = 25)
        num_na <- sample(1:10, size = sample(0:2, 1))
        df[[paste0("sample", i)]][num_na] <- NA
    }
    df <- as.data.frame(df)


    expect_snapshot(bppg::calculatePeptideRatios(aggr_intensities = df, id_cols = 1))

    expect_snapshot(bppg::calculatePeptideRatios(aggr_intensities = df, id_cols = 1, type = "difference"))
})





