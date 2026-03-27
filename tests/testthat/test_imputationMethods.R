test_that("test min 2 impute", {

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

    imputed <- list()
    for (i in 1:3) {
        imputed[[i]] <- t(bppg::min_2_impute(df[, i * c(1:3)], df))
    }

    D_min <- as.data.frame(imputed)


    expect_snapshot(D_min)
})