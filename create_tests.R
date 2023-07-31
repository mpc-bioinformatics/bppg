library(bppg)

library(seqinr)
file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
names(fasta) <- limma::strsplit2(names(fasta), "\\|")[,2]
res <- digest_fasta(fasta)

#saveRDS(res, file = "inst/extdata/digested_proteins_test.rds")
saveRDS(res, file = "tests/testthat/digested_proteins_test.rds")


edgelist <- generate_edgelist(res)
#saveRDS(edgelist, file = "inst/extdata/edgelist_test.rds")
saveRDS(edgelist, file = "tests/testthat/edgelist_test.rds")



edgelist_coll <- bppg::collapse_edgelist(edgelist,
                                         collapse_protein_nodes = TRUE,
                                         collapse_peptide_nodes = TRUE)
#saveRDS(edgelist_coll, file = "inst/extdata/edgelist_coll_pept_prot_test.rds")
saveRDS(edgelist_coll, file = "tests/testthat/edgelist_coll_pept_prot_test.rds")

graphs <- bppg::generate_graphs_from_edgelist(edgelist_coll)
#saveRDS(graphs, file = "inst/extdata/graphs_coll_pept_prot_test.rds")
saveRDS(graphs, file = "tests/testthat/graphs_coll_pept_prot_test.rds")



edgelist_coll2 <- bppg::collapse_edgelist(edgelist,
                                         collapse_protein_nodes = TRUE,
                                         collapse_peptide_nodes = FALSE)
#saveRDS(edgelist_coll2, file = "inst/extdata/edgelist_coll_prot_test.rds")
saveRDS(edgelist_coll2, file = "tests/testthat/edgelist_coll_prot_test.rds")

graphs2 <- bppg::generate_graphs_from_edgelist(edgelist_coll2)
#saveRDS(graphs2, file = "inst/extdata/graphs_coll_prot_test.rds")
saveRDS(graphs2, file = "tests/testthat/graphs_coll_prot_test.rds")


################################################################################

library(testthat)
library(igraph)
library(bppg)

G <- readRDS(test_path("testfiles/graphs_coll_pept_prot_test.rds"))





res <- bppg::calculate_subgraph_characteristics(S = G, #S2, S3,
                                               fastalevel = TRUE,
                                               #comparison = NULL,
                                               file = NULL)
res


expected <- data.frame(
  graph_ID = c(1,2),
  nr_protein_nodes = c(7,1),
  nr_peptide_nodes = c(13,1),
  nr_unique_peptid_nodes = c(7,1),
  nr_shared_peptid_nodes = c(6,0),
  nr_edges = c(22,1),
  nr_protein_accession = c(7,1),
  nr_peptide_sequences = c(476, 204),
  comparison = c(1,1)
)




