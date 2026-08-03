
# Test-files for examples and tests:

# original FASTA file: Uniprot S. cerevisiae proteome, downloaded on 2025-08-20, version 2025-03

# 9 Proteins are selected, which will lead to 4 graphs with different shapes

library(seqinr)
file <- "inst/original_data/uniprotkb_proteome_Scerevisiae_UP000002311_20250820_v202503.fasta"
fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)

proteins <- c("P09938",
              "P39708", "P07262",
              "P40212", "Q12690",
              "P00330", "P00331", "P07246", "P38113")

# get all fasta entries for the respective 9 proteins
ind <- NULL
for(i in seq_along(proteins)) {
    ind <- c(ind, grep(proteins[[i]], names(fasta)))
}

fasta_filtered <- fasta[ind]
seqinr::write.fasta(sequences = fasta_filtered, names = names(fasta_filtered),
                    file.out = "inst/extdata/uniprot_proteome_Scerevisiae_filtered.fasta")



# The file "peptides.txt" is the peptide quantity table generated via
# MaxQuant. The raw data were downloaded from the PRIDE repository (PXD001819)
# and processed with MaxQuant version 2.7.3.0.
# It was then filtered for rows belonging to peptides that are present in the
# 9 selected proteins.


file <- "inst/original_data/peptides.txt"
D <- read.table(file, sep = "\t", header = TRUE)

# find rows that belong to proteins associated with the 9 proteins
ind <- NULL
for(i in seq_along(proteins)) {
    ind <- c(ind, grep(proteins[[i]], D$Proteins))
}

D_filtered <- D[unique(ind),]
write.table(D_filtered, file = "inst/extdata/peptides_filtered.txt", sep = "\t", row.names = FALSE)


################################################################################
### generate theoretical graphs

file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta", package = "bppg")
fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
edgelist <- digestFASTA(fasta)
graphs <- generateGraphsFromEdgelist(edgelist, collProtNodes = FALSE, collPeptNodes = FALSE)
graphs_collPeptProt <- generateGraphsFromEdgelist(edgelist, collProtNodes = TRUE, collPeptNodes = TRUE)
graphs_collProt <- generateGraphsFromEdgelist(edgelist, collProtNodes = TRUE, collPeptNodes = FALSE)

saveRDS(graphs, file = "inst/extdata/theoGraphs.rds")
saveRDS(graphs_collPeptProt, file = "inst/extdata/theoGraphs_collpeptprot.rds")
saveRDS(graphs_collProt, file = "inst/extdata/theoGraphs_collprot.rds")


################################################################################
### generate quant graphs

library(seqinr)
file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta", package = "bppg")
fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
edgelist <- digestFASTA(fasta)

file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
group <- factor(rep(1:9, each = 3))
D <- readMqPeptideTable(path = file, group = group, LFQ = FALSE, remove_contaminants = FALSE)
D_norm <- normalizePeptideIntensities(D)
dAgg <- aggregateReplicates(D_norm, group = group)
exp_peptide_ratios <- calculatePeptideRatios(dAgg)

# graphs with collapsed protein nodes and NOT collapsed peptide nodes (should be used for optimization)
res <- generateQuantGraphs(exp_peptide_ratios, edgelist)
saveRDS(res, file = "inst/extdata/quantGraphs.rds")


# graphs with collapsed protein nodes and collapsed peptide nodes (should be used for visualization)
res <- generateQuantGraphs(exp_peptide_ratios, edgelist, collPeptNodes = TRUE)
saveRDS(res, file = "inst/extdata/quantGraphs_collpept.rds")


################################################################################
### generate result objects for testing

file <- system.file("extdata", "quantGraphs.rds", package = "bppg")
graphs <- readRDS(file)

RES <- list()
for (comp in seq_along(graphs)) {
    comparison <- names(graphs)[comp]
    graphs_tmp <- graphs[[comp]]
    RES_tmp <- NULL
    for (i in 1:3) {  # only first  graphs per comparison to speed up testing
        G <- graphs_tmp[[i]]
        res <- iterateOverCi(G, gridSize = 100)
        if (is.null(RES_tmp)) {
            RES_tmp <- automatedAnalysisIteratedCi(G, res)
        } else {
            RES_tmp <- rbind(RES_tmp, automatedAnalysisIteratedCi(G, res))
        }
    }
    RES <- c(RES, RES_tmp)
}
names(RES) <- names(graphs)

resultsList <- RES

saveRDS(resultsList, file = "inst/extdata/resultsList.rds")




