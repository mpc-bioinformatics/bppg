test_that("test .digest2", {
    file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
    fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

    digested_proteins1 <- bppg:::.digest2(sequence = fasta[[1]],
                                            enzyme = "trypsin",
                                            missed = 0,
                                            remove_initial_M = TRUE)

    # test sequence that ends on R or K
    digested_proteins2 <- bppg:::.digest2(sequence = paste0(fasta[[1]], "R"),
                                            enzyme = "trypsin",
                                            missed = 2,
                                            remove_initial_M = FALSE)
    ## test if undefined enzyme is being used
    expect_error(bppg:::.digest2(fasta[[1]], enzyme = "LysC")) 
     # warning for no cleavage site
    expect_warning(bppg:::.digest2("ABCD"))
    # warning for not possible 2 missed cleavage
    expect_warning(bppg:::.digest2("ABCRD", missed = 2)) 

    expect_snapshot(digested_proteins1)
    expect_snapshot(digested_proteins2)
})


test_that("digestion of a FASTA file", {
    file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
    fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
    names(fasta) <- limma::strsplit2(names(fasta), "\\|")[,2]
    res <- digestFASTA(fasta)

    protOrigin <- as.list(c(rep("yeast", 4), rep("spike_in", 3)))
    res2 <- digestFASTA(fasta, protOrigin = protOrigin)

    expect_snapshot(res)
    expect_snapshot(res2)
})


