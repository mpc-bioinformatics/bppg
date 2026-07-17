# TODO run exmpales with local data
test_that("min2impute", {

    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    D <- bppg::readMqPeptideTable(file)
    D_norm <- bppg::normalizePeptideIntensities(D)

    df <- SummarizedExperiment::assays(D_norm)$intensities_norm
    group <- rep(1:9, each = 3)

    imputed <- list()
    for (g in factor(group)) {
        imputed[[g]] <- t(bppg:::.min2impute(D = df[g == group],
                intensities = df))
    }

    D_min <- as.data.frame(imputed)

    expect_snapshot(D_min)
})

test_that("colMeanImputation", {

    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    D <- bppg::readMqPeptideTable(file)
    D_norm <- bppg::normalizePeptideIntensities(D)

    df <- SummarizedExperiment::assays(D_norm)$intensities_norm
    group <- rep(1:9, each = 3)

    imputed <- list()
    for (g in factor(group)) {
        imputed[[g]] <- bppg:::.colMeanImputation(D = df[g == group],
                intensities = df)
    }

    D_mean <- data.frame(imputed)
    colnames(D_mean) <- colnames(df)

    expect_snapshot(D_mean)
})

test_that("missForest", {

    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    D <- bppg::readMqPeptideTable(file)
    D_norm <- bppg::normalizePeptideIntensities(D)

    df <- SummarizedExperiment::assays(D_norm)$intensities_norm

    set.seed(8)

    D_missForest <- bppg:::.missForest(df)

    expect_snapshot(D_missForest)
})



test_that("QRILC", {
    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    D <- bppg::readMqPeptideTable(file)
    D_norm <- bppg::normalizePeptideIntensities(D)

    df <- SummarizedExperiment::assays(D_norm)$intensities_norm

    set.seed(42)
    imputed <- bppg:::.QRILC(df) # completly random?

    expect_snapshot(imputed)
})


test_that("BPCA", {
    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    D <- bppg::readMqPeptideTable(file)
    D_norm <- bppg::normalizePeptideIntensities(D)

    df <- SummarizedExperiment::assays(D_norm)$intensities_norm

    D_BPCA <- bppg:::.BPCA(df, verbose = FALSE)

    expect_snapshot(D_BPCA)
})
