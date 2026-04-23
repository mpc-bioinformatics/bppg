
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
res <- generateQuantGraphs(exp_peptide_ratios, edgelist)
saveRDS(res, file = "inst/extdata/quantGraphs_collpept.rds")




file <- system.file("extdata", "quantGraphs_collpept.rds", package = "bppg")
graphs <- readRDS(file)
G <- graphs$"1_2"[[3]]

plotBipartiteGraph(G, three_shapes = TRUE, useCanonicalPermutation = TRUE, legend.x = 0)











################################################################################
#### testing around

M1 <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
G1 <- igraph::graph_from_biadjacency_matrix(M1)
plot(G1, layout = igraph::layout_as_bipartite)

M2 <- matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE)
G2 <- igraph::graph_from_biadjacency_matrix(M2)
plot(G2, layout = igraph::layout_as_bipartite)


bppg:::.isomorphicBipartite(G1, G2)

#' M1 <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
#' G1 <- igraph::graph_from_biadjacency_matrix(M1)
#'
#' M2 <- matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE)
#' G2 <- igraph::graph_from_biadjacency_matrix(M2)
#'
#' .isomorphicBipartite(G1, G2)
#'

