
# Test-files for examples:

# original FASTA file: Uniprot S. cerevisiae proteome, downloaded on 2025-08-20, version 2025-03

# 5 Proteins are selected, which will lead to 3 graphs with different shapes
# The 5 proteins are:
# P09938
# P39708
# P07262
# P40212
# Q12690

library(seqinr)
file <- system.file("extdata", "uniprotkb_proteome_Scerevisiae_UP000002311_20250820_v202503.fasta", package = "bppg")
fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

proteins <- c("P09938", "P39708", "P07262", "P40212", "Q12690")

# get all fasta entries for the respective 5 proteins
ind <- NULL
for(i in seq_along(proteins)) {
    ind <- c(ind, grep(proteins[[i]], names(fasta)))
}

fasta_filtered <- fasta[ind]
seqinr::write.fasta(sequences = fasta_filtered, names = names(fasta_filtered),
                    file.out = "inst/extdata/uniprot_proteome_Scerevisiae_filtered.fasta")

# TODO: explain the origin of .raw files and quantification via MaxQuant


file <- system.file("extdata", "peptides.txt", package = "bppg")
D <- read.table(file, sep = "\t", header = TRUE)

# find rows that belong to proteins associated with the 5 proteins
ind <- NULL
for(i in seq_along(proteins)) {
    ind <- c(ind, grep(proteins[[i]], D$Proteins))
}

D_filtered <- D[unique(ind),]
write.table(D_filtered, file = "inst/extdata/peptides_filtered.txt", sep = "\t", row.names = FALSE)



################################################################################
### generate quant graphs

library(seqinr)
file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta", package = "bppg")
fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
edgelist <- digestFASTA(fasta)

file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
group <- factor(rep(1:9, each = 3))
D <- readMqPeptideTable(path = file, group = group, LFQ = FALSE, remove_contaminants = FALSE)

dAgg <- aggregateReplicates(D, group = group)
exp_peptide_ratios <- calculatePeptideRatios(dAgg)

# graphs with collapsed protein nodes and NOT collapsed peptide nodes (should be used for optimization)
res <- generateQuantGraphs(exp_peptide_ratios, edgelist)
saveRDS(res, file = "inst/extdata/quantGraphs.rds")


# graphs with collapsed protein nodes and collapsed peptide nodes (should be used for visualization)
res <- generateQuantGraphs(exp_peptide_ratios, edgelist, collPeptNodes = TRUE)
saveRDS(res, file = "inst/extdata/quantGraphs_collpept.rds")

plotBipartiteGraph(res[[1]][[3]], legend = FALSE)





