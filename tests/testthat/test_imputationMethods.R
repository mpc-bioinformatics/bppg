# TODO run exmpales with local data
test_that(".min2impute", {

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
        imputed[[i]] <- t(bppg:::.min2impute(D = df[seq(i * 3 - 2, i * 3)],
                intensities = df))
        names(imputed)[[i]] <- paste0("sample", i)
    }

    D_min <- as.data.frame(imputed)


    expect_snapshot(D_min)
})

test_that(".missForest", {

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

    set.seed(8)

    imputed <- list()
    for (i in 1:3) {
        imputed[[i]] <- t(missForest::missForest(
                        df[seq(i * 3 - 2, i * 3)],
                         verbose = FALSE
                        )$ximp)
        names(imputed)[[i]] <- paste0("sample", i)
    }

    D_missForest <- as.data.frame(imputed)

    expect_snapshot(D_missForest)
})



test_that(".QRILC", {

    df <- list()
    df <- c(df, sequence = list(paste0("pep_", 1:10)))

    for (i in 1:3) {
        for (j in 1:3) {
            set.seed((i + 2)^(j + 2))
            df[[paste0("sample", i, "_run", j)]] <- runif(10, 15, 25)
            num_na <- sample(1:10, size = sample(1:3, 1))
            df[[paste0("sample", i, "_run", j)]][num_na] <- NA
        }
    }

    df <- as.data.frame(df)
    rownames(df) <- df$sequence
    df$sequence <- NULL

    imputed <- list()

    for (i in 1:3) {

        df_subset <- as.matrix(df[, seq(i * 3 - 2, i * 3)])

        res <- bppg:::.QRILC(df_subset)

        imp <- if (is.list(res)) {
            if (!is.null(res$imp)) res$imp else res[[1]]
        } else {
            res
        }

        imputed[[i]] <- t(as.matrix(imp))
    }

    D_QRLIC <- as.data.frame(imputed)

    expect_snapshot(D_QRLIC)
})


test_that(".BPCA", {

    df <- list()
    df <- c(df, sequence = list(paste0("pep_", 1:10)))

    for (i in 1:3) {
        for (j in 1:3) {
            set.seed((i + 2)^(j + 2))

            df[[paste0("sample", i, "_run", j)]] <- runif(10, 15, 25)

            num_na <- sample(1:10, size = sample(1:3, 1))
            df[[paste0("sample", i, "_run", j)]][num_na] <- NA
        }
    }

    df <- as.data.frame(df)
    rownames(df) <- df$sequence
    df$sequence <- NULL

    imputed <- list()

    for (i in 1:3) {

        df_subset <- as.matrix(df[, seq(i * 3 - 2, i * 3)])

        imp <- bppg:::.BPCA(
            D = df_subset,
            intensities = df
        )

        imputed[[i]] <- t(imp)
        names(imputed)[[i]] <- paste0("sample", i)
    }

    D_BPCA <- as.data.frame(imputed)

    expect_snapshot(D_BPCA)
})