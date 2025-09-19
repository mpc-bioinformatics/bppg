

test_that("digestion function", {

    P1 <- "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"

    ## no missed cleavages
    digested_P1 <- bppg:::.digest2(P1)
    expect_equal(nrow(digested_P1), 15)
    expect_equal(ncol(digested_P1), 4)
    expect_equal(digested_P1$mc, rep(0, 15))

    ## up to 2 missed cleavages
    digested_P1_mc2 <- bppg:::.digest2(P1, missed = 2)
    expect_equal(nrow(digested_P1_mc2), 42)
    expect_equal(ncol(digested_P1_mc2), 4)
    expect_equal(digested_P1_mc2$mc, c(rep(0, 15), rep(1, 14), rep(2, 13)))

    ## digest a protein with K or R as last amino acid
    P2 <- "MVHLTPEEKSAVTALWGKVNVDEVGGEALGR"
    digested_P2 <- bppg:::.digest2(P2)

    ## test if undefined enzyme is being used
    expect_error(bppg:::.digest2(P2, enzyme = "LysC"))

    ## test trypsin.strict enzyme
    P3 <- "MVHLTPEEKSAVTALWGKPVNVDEVGGEALGR"
    digested_P3 <- bppg:::.digest2(P3, enzyme = "trypsin")
    digested_P3_strict <- bppg:::.digest2(P3, enzyme = "trypsin.strict")
    expect_in("SAVTALWGKPVNVDEVGGEALGR", digested_P3$sequence)
    expect_in("SAVTALWGK", digested_P3_strict$sequence)

    expect_warning(bppg:::.digest2("ABCD")) # warning for no cleavage site
    expect_warning(bppg:::.digest2("ABCRD", missed = 2)) # warning for not possible 2 missed cleavage

})





test_that("digestion of a FASTA file", {
    digested_proteins <- readRDS(testthat::test_path("testfiles/digested_proteins_test.rds"))

    file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
    fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
    names(fasta) <- limma::strsplit2(names(fasta), "\\|")[,2]
    res <- digestFASTA(fasta)

    expect_equal(res, digested_proteins)
})



