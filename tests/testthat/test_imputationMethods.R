
test_that("min2impute", {

    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    D <- bppg::readMqPeptideTable(file)
    D_norm <- bppg::normalizePeptideIntensities(D)

    df <- SummarizedExperiment::assays(D_norm)$intensities_norm
    group <- rep(1:9, each = 3)

    imputed <- list()
    for (g in levels(factor(group))) {
        imputed[[g]] <- bppg:::.min2impute(D = df[g == group],
                intensities = df)
    }

    D_min <- data.frame(imputed) # aggregated values
    colnames(D_min) <- levels(factor(group))

    expect_snapshot(D_min)
})

test_that("colMeanImputation", {

    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    D <- bppg::readMqPeptideTable(file)
    D_norm <- bppg::normalizePeptideIntensities(D)

    df <- SummarizedExperiment::assays(D_norm)$intensities_norm
    D_mean <- bppg:::.colMeanImputation(intensities = df)

    expect_snapshot(D_mean)
})

test_that("missForest", {

    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    D <- bppg::readMqPeptideTable(file)
    D_norm <- bppg::normalizePeptideIntensities(D)

    df <- SummarizedExperiment::assays(D_norm)$intensities_norm

    set.seed(8)

    D_missForest <- data.frame(bppg:::.missForest(df))
    colnames(D_missForest) <- colnames(df)

    expect_snapshot(D_missForest)
})



test_that("QRILC", {
    file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
    D <- bppg::readMqPeptideTable(file)
    D_norm <- bppg::normalizePeptideIntensities(D)

    df <- SummarizedExperiment::assays(D_norm)$intensities_norm

    set.seed(42)
    imputed <- bppg:::.QRILC(df)

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
