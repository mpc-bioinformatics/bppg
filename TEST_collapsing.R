library(seqinr)
file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
graphs <- bppg::generate_graphs_from_FASTA(fasta)


digested_proteins <- digest_fasta(fasta)
edgelist <- generate_edgelist(digested_proteins)


collapse_protein_nodes = TRUE
collapse_peptide_nodes = TRUE

X <- collapse_edgelist(edgelist)


Y <- collapse_edgelist(edgelist, collapse_peptide_nodes = FALSE)


Z <- collapse_edgelist(edgelist, collapse_protein_nodes = FALSE)
