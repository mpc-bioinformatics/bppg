test_that("test .digest2", {
    file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta",
                        package = "bppg")
    fasta <- Biostrings::readAAStringSet(file)

    digested_proteins1 <- bppg:::.digest2(sequence = as.character(fasta[[1]]),
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
    file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta",
                        package = "bppg")
    fasta <- Biostrings::readAAStringSet(file) 
    names(fasta) <- vapply(strsplit(names(fasta), "|", fixed = TRUE), "[", 
        character(1), 2)
    res <- digestFASTA(fasta)

    # fake protOrigin just for testing and combiend input
    fasta_tmp <- Biostrings::readAAStringSet(file)
    names(fasta_tmp) <- vapply(strsplit(names(fasta_tmp), "|", fixed = TRUE), "[", 
        character(1), 2)
    fasta <- c(fasta, fasta_tmp)
    protOrigin <- as.list(c(rep("yeast", 7), rep("spike_in", 2), rep("neither", 9)))
    res2 <- digestFASTA(fasta, protOrigin = protOrigin)

    # test list input
    fasta2 <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)      
    names(fasta2) <- vapply(strsplit(names(fasta2), "|", fixed = TRUE), "[", 
        character(1), 2)
    res3 <- digestFASTA(fasta2)

    expect_snapshot(res)
    expect_snapshot(res2)
    expect_snapshot(res3)
})


