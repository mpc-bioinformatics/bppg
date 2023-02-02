# library(BBmisc)    # for collapse()
# library(seqinr)    # for reading in fasta files
# library(limma)     # for strsplit2()
# library(pbapply)   # for progress bars for apply functions
# library(igraph)    # for graph functionality
# library(openxlsx)
# library(tidyverse) # long format + ggplot2

library(bppg)

################################################################################
################################################################################
### prepare edgelist from FASTA file

#source("R Scripts/bipartite_graph_generation_and_analysis/helper_functions/Digest2.R")

# C:/Users/schorkka/UNI/Promotion/promotion_project/

fasta1 <- seqinr::read.fasta(file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_fasta/20220722_uniprot_proteome_UP000005640_HomoSapiens_v2022_02.fasta",
                     seqtype = "AA", as.string = TRUE)
names(fasta1) <- limma::strsplit2(names(fasta1), "\\|")[,2]
fasta2 <- seqinr::read.fasta(file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_fasta/20220722_uniprot_proteome_UP000000589_MusMusculus_v2022_02.fasta",
                     seqtype = "AA", as.string = TRUE)
names(fasta2) <- limma::strsplit2(names(fasta2), "\\|")[,2]
fasta3 <- seqinr::read.fasta(file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_fasta/MaxQuant_contaminants_v2_3_1.fasta",
                     seqtype = "AA", as.string = TRUE)
names(fasta3) <- paste0("CON_", names(fasta3))

fasta <- c(fasta1, fasta2, fasta3)
prot_origin <- c(rep("human", length(fasta1)), rep("mouse", length(fasta2)), rep("contaminant", length(fasta3)))



digested_proteins <- bppg::digest_fasta(fasta, missed_cleavages = 2, min_aa = 7, max_aa = 50)
# 4min

protein_list <- names(digested_proteins)

edgelist <- generate_edgelist(digested_proteins)
# 36s

write.table(edgelist, "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_fasta/edgelist.txt", sep = "\t", row.names = FALSE)



################################################################################
################################################################################
## read in quantitative peptide data

DATA <- bppg::read_MQ_peptidetable("C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/peptides.txt",
                                   LFQ = TRUE, remove_contaminants = FALSE)

validvalues <- rowSums(!is.na(DATA[,-1]))
barplot(table(validvalues))

### only keep peptides with at least one valid value

DATA_filtered <- DATA[validvalues >= 1, ]





################################################################################
### map peptides to proteins (by using the edgelist from the fasta file!)

edgelist <- read.table("C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_fasta/edgelist.txt",
                       sep = "\t", header = TRUE)
peptides <- DATA_filtered$Sequence

edgelist_filtered <- edgelist[edgelist[,2] %in% peptides, ]



#proteins <- character(length(peptides))

proteins <- pbapply::pblapply(peptides, function(x, edgelist, prot_origin) {
  prot <- edgelist[,1][edgelist[,2] == x]

  ### TODO: mixed könnte auch gemischt zwischen contaminant und human sein?

  ind <- which(protein_list %in% prot)
  prot_origin_tmp <- unique(prot_origin[ind])
  if (length(prot_origin_tmp) > 1) prot_origin_tmp <- "mixed"
  if (length(prot_origin_tmp) == 0) prot_origin_tmp <- ""


  return(c(proteins = BBmisc::collapse(prot, "/"), prot_origin = prot_origin_tmp))
}, edgelist = edgelist_filtered, prot_origin = prot_origin)
### 14min

proteins2 <- BBmisc::convertListOfRowsToDataFrame(proteins)


#peptides[which(proteins == "")]  ### diese Peptide sind >50 AA lang!
D <- cbind(peptides = DATA_filtered$Sequence, proteins = proteins2$proteins,
           prot_origin = proteins2$prot_origin, DATA_filtered[,-1])
D <- D[D$proteins != "",]


write.table(D, "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/preprocessed_peptide_data_D4.txt",
            sep = "\t", row.names = FALSE)




################################################################################
### calculate peptide ratios

#source("R Scripts/bipartite_graph_generation_and_analysis/helper_functions/calculate_peptide_ratios.R")

D <- read.table(file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/preprocessed_peptide_data_D4.txt", header = TRUE)
group <- factor(limma::strsplit2(colnames(D)[-c(1:3)], "_")[,1])

D_aggr <- bppg::aggregate_replicates(D, method = "mean", missing.limit = 0.4,
                                        group = group, id_cols = 1:3)

write.table(D_aggr, file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/aggr_pept_int_D4.txt",
            row.names = FALSE, sep = "\t")
save(group, file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/group.RData")

D_aggr <- read.table("C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/aggr_pept_int_D4.txt",
                     sep = "\t", header = TRUE)
groups  <- levels(group)



peptide_ratios <- bppg::calculate_peptide_ratios(aggr_intensities = D_aggr, id_cols = 1:3, group_levels = groups)
write.table(peptide_ratios, file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/peptide_ratios.txt",
            row.names = FALSE, sep = "\t")



################################################################################
### generate graphs

peptide_ratios <- read.table(file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/peptide_ratios.txt",
                 header = TRUE, sep = "\t")

edgelist <- read.table("C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_fasta/edgelist.txt", sep = "\t", header = TRUE)

subgraphs <- bppg::generate_quant_graphs(peptide_ratios = peptide_ratios, id_cols = 1:3, fasta_edgelist=edgelist)
saveRDS(subgraphs, "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/subgraphs.rds")



################################################################################
### collapse peptide and protein ratios
subgraphs <- readRDS("C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/subgraphs.rds")

subgraphs_coll_prot <- lapply(subgraphs, bppg::collapse_protein_nodes, sparse = FALSE, fast = FALSE, fc = TRUE)

saveRDS(subgraphs_coll_prot, "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/subgraphs_coll_prot.rds")


subgraphs_coll_prot_pep <- lapply(subgraphs_coll_prot, bppg::collapse_peptide_nodes, sparse = FALSE, fc = TRUE, fast = TRUE)

saveRDS(subgraphs_coll_prot_pep, "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/subgraphs_coll_prot_pep.rds")



################################################################################
################################################################################
################################################################################
#### calculate isomorph lists separately per comparison

#library(igraph)
#source("R Scripts/bipartite_graph_generation_and_analysis/helper_functions/isomorph_classes_calculation_and_plotting_functions.R")

S <- readRDS("C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/subgraphs_coll_prot_pep.rds")
#S2 <- readRDS("C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/Subgraphs_with_pep_ratios_mergedProteins.rds")
#names(S) <- names(S2)

isomorph_list <- list()

for (i in 1:21) {
  print(i)
  S_tmp <- S[[i]]
  isomorph_list[[i]] <- bppg::calculateIsomorphList(S_tmp, matrix = FALSE)
}
### TODO: relativ langsam

names(isomorph_list) <- names(S)

save(isomorph_list, file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/isomorph_classes/isomorph_classes_coll_prot_pep_D4.RData")


#### TODO: hier weitermachen



################################################################################
#### calculate isomorph lists together over all comparisons

lengths(S)
# [1] 4498 4424 3545 3594 3114 2812 4622 3796 3825 3347 3040 4263 4436 3997 3734 4351 3940 3665 4803 4549 4595
S_all <- unlist(S, recursive = FALSE)

isomorph_all <- bppg::calculateIsomorphList(S_all, matrix = FALSE)
save(isomorph_all, file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/isomorph_classes/isomorph_classes_all.RData")






################################################################################

source("R Scripts/bipartite_graph_generation/helper_functions/calculate_subgraph_characteristics_2.R")

G <- readRDS("data/D4_without_isoforms/D4_quant/preprocessed/Subgraphs_mergedPeptides.rds")

RES <- NULL

for(i in 1:21) {

  G_tmp <- G[[i]]
  D_tmp <- calculate_subgraph_characteristics_2(G_tmp)
  D_tmp$comparison <- rep(names(G)[[i]], nrow(D_tmp))

  RES <- rbind(RES, D_tmp)
}


write.xlsx(RES, "data/D4_without_isoforms/D4_quant/table_subgraph_characteristics_D4_quant.xlsx")


#################################################################################
#################################################################################
### plot prototypes of the isomorphism classes

prototypes_D4 <- readRDS("data/D4_without_isoforms/D4_quant/isomorph_classes/prototypes_D4.rds")

pdf("data/D4_without_isoforms/D4_quant/isomorph_classes/prototype_plots.pdf", height = 10, width = 10)

for (i in 1:length(prototypes_D4)) {
  print(i)
  bppg::plotBipartiteGraph(prototypes_D4[[i]], legend = FALSE, vertex.size = 28, vertex.label.cex = 4,
                           vertex.size2=28, edge.width = 2, three_shapes = TRUE)
}

dev.off()








################################################################################
### quality control of the peptide data

source("R Scripts/QC_and_normalization/automatedNormalization_v1_3.R")
source("R Scripts/QC_and_normalization/MA_Plots_v1_3.R")
source("R Scripts/QC_and_normalization/PCA_plot_v1_3.R")
source("R Scripts/QC_and_normalization/ValidValue_Plot_v1_3.R")


LFQ_values_filtered2 <- LFQ_values_filtered
colnames(LFQ_values_filtered2) <- substr(colnames(LFQ_values_filtered2), 15, 100)
group <- factor(substr(colnames(LFQ_values_filtered2), 1, 7))


LFQ_values_filtered2_long <- pivot_longer(LFQ_values_filtered2, 1:20)
LFQ_values_filtered2_long$group <- factor(substr(LFQ_values_filtered2_long$name, 1, 7))

Boxplots(X = LFQ_values_filtered2_long, groupvar_name = "Group", sample_filter = NULL,
         plot_device = "png", suffix = "LFQ_peptides",
         plot_height = 10, plot_width = 15, plot_dpi = 300,
         log_data = TRUE, log_base = 2, group_colours = NULL,
         method = "boxplot", base_size = 20,
         output_path = "data/D4_without_isoforms/D4_quant/QC_and_Norm/")



ValidValuePlot(X = LFQ_values_filtered2_long, groupvar_name = "Group", sample_filter=NULL,
               suffix = "LFQ_peptides", plot_device = "png",
               group_colours = NULL, plot_height = 10,
               plot_width = 15, plot_dpi = 300, ylim = NULL, title = NULL,
               output_path = "data/D4_without_isoforms/D4_quant/QC_and_Norm/")


### TODO: Einfärbungen der verschiedenen Peptid-Arten (human, mouse, shared!)
MAPlots(X = LFQ_values_filtered2, log = TRUE, alpha = TRUE, suffix="LFQ_peptides",
        plot_height=15, plot_width=15, output_path = "data/D4_without_isoforms/D4_quant/QC_and_Norm/")


### TODO: hier wäre wahrscheinlich eine LTS-Normalisierung angebracht, weil sich
###       die Proben durch die verschiedenen Mischungen stark ändern?


PCA_Plot(X = LFQ_values_filtered2, id = NULL, log_data = TRUE, log_base = 2,
         impute = TRUE, impute_method = "mean", propNA = 0.3,
         scale. = TRUE,
         groupvar1 = group, groupvar2 = NULL, groupvar1_name = "Group", groupvar2_name = NULL,
         point.size = 4, base_size = 11,
         group_colours = NULL, returnPCA = FALSE, title = NULL,
         plot_device = "png", plot_height = 10, plot_width = 10, plot_dpi = 300,
         output_path = "data/D4_without_isoforms/D4_quant/QC_and_Norm/", suffix = "LFQ_peptides_imputed", xlim = NULL, ylim = NULL,
         label = FALSE, PCx = 1, PCy = 2, label_size = 4)

PCA_Plot(X = LFQ_values_filtered2, id = NULL, log_data = TRUE, log_base = 2,
         impute = FALSE, impute_method = "mean", propNA = 0.5,
         scale. = TRUE,
         groupvar1 = group, groupvar2 = NULL, groupvar1_name = "Group", groupvar2_name = NULL,
         point.size = 4, base_size = 11,
         group_colours = NULL, returnPCA = FALSE, title = NULL,
         plot_device = "png", plot_height = 20, plot_width = 20, plot_dpi = 300,
         output_path = "data/D4_without_isoforms/D4_quant/QC_and_Norm/", suffix = "LFQ_peptides", xlim = NULL, ylim = NULL,
         label = TRUE, PCx = 1, PCy = 2, label_size = 4)

