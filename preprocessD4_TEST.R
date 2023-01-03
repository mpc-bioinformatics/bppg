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
names(fasta1) <- strsplit2(names(fasta1), "\\|")[,2]
fasta2 <- seqinr::read.fasta(file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_fasta/20220722_uniprot_proteome_UP000000589_MusMusculus_v2022_02.fasta",
                     seqtype = "AA", as.string = TRUE)
names(fasta2) <- strsplit2(names(fasta2), "\\|")[,2]
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

write.table(edgelist, "data/D4_without_isoforms/D4_fasta/edgelist.txt", sep = "\t", row.names = FALSE)

################################################################################
################################################################################
## read in quantitative peptide data

# library(BBmisc)    # for collapse()
# library(seqinr)    # for reading in fasta files
# library(limma)     # for strsplit2()
# library(pbapply)   # for progress bars for apply functions
# library(igraph)


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

D_mean_NA_mind2 <- aggregate_replicates(D, method = "mean", use0 = FALSE, missing.limit = 0.4,
                                        group = group, accession.cols = 1:2)

write.table(D_mean_NA_mind2, file = "data/D4_without_isoforms/D4_quant/preprocessed/D4_mean_NA_mind2.txt",
            row.names = FALSE, sep = "\t")


D_aggr <- read.table("data/D4_without_isoforms/D4_quant/preprocessed/D4_mean_NA_mind2.txt", sep = "\t", header = TRUE)
groups  <- levels(group)

FC <- NULL

for (i in 1:(length(groups)-1)) {
  lev1 <- groups[[i]]
  for(j in (i+1):length(groups)){
    lev2 <- groups[[j]]
    name <- paste0("FC_", lev1, "_", lev2)
    FC <- cbind(FC, foldChange(D_aggr, lev1, lev2))
    colnames(FC)[ncol(FC)] <- name
  }

}



FC <- cbind(D_aggr[, 1:2], FC)
write.table(data.frame(FC, prot_origin), file = "data/D4_without_isoforms/D4_quant/preprocessed/FC_NA_mind2.txt", row.names = FALSE, sep = "\t")




### TODO: hier weitermachen!!


################################################################################

library(igraph)

FC <- read.table(file = "data/D4_without_isoforms/D4_quant/preprocessed/FC_NA_mind2.txt", header = TRUE, sep = "\t")

edgelist <- read.table("data/D4_without_isoforms/D4_fasta/edgelist.txt", sep = "\t")
edgelist_filtered <- edgelist[edgelist[,2] %in% FC$peptides, ]


graphs <- list()
subgraphs <- list()
#i = 3

for (i in 3:23) {

  comparison <- substr(colnames(FC)[i], 4 , 100)


  fc <- FC[,i]
  peptides_tmp <- FC$peptides[!is.na(fc)]  ## peptides that are quantified in this specific comparison
  fc <- na.omit(fc)
  edgelist_filtered2 <- edgelist_filtered[edgelist_filtered[,2] %in% peptides_tmp, ]

  G_tmp <- igraph::graph_from_edgelist(as.matrix(edgelist_filtered2[, 1:2]), directed = FALSE)
  V(G_tmp)[name %in% edgelist_filtered2[,1]]$type <- TRUE
  V(G_tmp)[name %in% edgelist_filtered2[,2]]$type <- FALSE

  vertex_attr(G_tmp, "pep_ratio", index = V(G_tmp)[!V(G_tmp)$type]) <- fc[match(V(G_tmp)$name[!V(G_tmp)$type], peptides_tmp)]

  subgraphs_tmp <- igraph::decompose(G_tmp)

  graphs[[i-2]] <- G_tmp
  names(graphs)[[i-2]] <- comparison

  subgraphs[[i-2]] <- subgraphs_tmp
  names(subgraphs)[[i-2]] <- comparison

}


saveRDS(graphs, file = "data/D4_without_isoforms/D4_quant/preprocessed/Graphs_with_pep_ratios.rds")
saveRDS(subgraphs, file = "data/D4_without_isoforms/D4_quant/preprocessed/Subgraphs_with_pep_ratios.rds")

x <- sapply(Subgraphs_with_pep_ratios[[1]], gorder)
## größter Graph ist 250!



################################################################################
################################################################################
### alte Version der Protein und Peptid-Zusammenfassung

source("R Scripts/bipartite_graph_generation/helper_functions/calculate_biadjacency_matrix_and_submatrices.R")
source("R Scripts/bipartite_graph_generation/helper_functions/duplicated_for_sparse_matrices.R")

library(igraph)


subgraphs <- readRDS("data/D4_without_isoforms/D4_quant/preprocessed/Subgraphs_with_pep_ratios.rds")


#subgraphs2 <- form_proteingroups(subgraphs, )


subgraphs2 <- lapply(subgraphs, form_proteingroups, sparse = FALSE, fast = TRUE, matrix = FALSE)

subgraphs3 <- lapply(subgraphs2, merge_Peptides, sparse = FALSE, fast = TRUE, matrix = FALSE, fc = FALSE)


saveRDS(subgraphs2, file = "data/D4_without_isoforms/D4_quant/preprocessed/Subgraphs_mergedProteins.rds")
saveRDS(subgraphs3, file = "data/D4_without_isoforms/D4_quant/preprocessed/Subgraphs_mergedPeptides.rds")


-####################### neue Version
  ################################################################################
################################################################################
#### TODO: proteine müssen zusammengefasst werden, Peptide aber nicht!

library(pbapply)

subgraphs <- readRDS("data/D4_without_isoforms/D4_quant/preprocessed/Subgraphs_with_pep_ratios.rds")

subgraphs_merged_proteins <- list()


for (l in 1:21) {

  print(l)


  subgraphs_tmp <- subgraphs[[l]]

  number_of_graphs <- length(subgraphs_tmp)

  subgraphs_with_merged_proteins <- list()
  merged_subgraphs <- vector(mode = 'list', length = number_of_graphs)

  #add progress bar to loop
  number_of_iterations <- number_of_graphs
  pb <- startpb(0, number_of_graphs)
  on.exit(closepb(pb))


  merged_subgraphs_tmp <- vector(mode = 'list', length = number_of_graphs)

  for(k in 1:(number_of_graphs)){

    #print(k)


    G <- subgraphs_tmp[[k]]

    peptide_ids <- as.numeric(V(G)[!V(G)$type])  # ids of peptide nodes (type = FALSE)
    protein_ids <- as.numeric(V(G)[V(G)$type])   # ids of protein nodes (type = TRUE)

    #assign initial number attribute to all vertices

    for(i in 1:gorder(G)){
      vertex_attr(G, "number", index = i) <- 1
    }

    # duplicated_peptides <- list() #list of peptide ids to be deleted due to merging
    #
    # if(length(peptide_ids) > 1){
    #   for(i in 1:(length(peptide_ids) - 1)){ # for all peptides i until the penultimate
    #     for(j in (i+1):length(peptide_ids)){ # for all peptides j>i
    #       vertex_id_1 <- peptide_ids[[i]] #id of peptide i
    #       vertex_id_2 <- peptide_ids[[j]] #id of peptide j
    #
    #       # if the vertex j has not been marked for deletion yet
    #       if(vertex_attr(G, "number", index = vertex_id_2) != 0){
    #
    #         # if vertices i and j have the same amount of edges
    #         if(length(G[[vertex_id_1]][[1]]) == length(G[[vertex_id_2]][[1]])){
    #
    #           # if all edges of i and j are identical
    #           if(!FALSE %in% unlist(unname(G[[vertex_id_1]][[1]] == G[[vertex_id_2]][[1]]))){
    #
    #             # raise number attribute of vertex i by 1 as it will be merged with vertex j
    #             vertex_attr(G, "number", index = vertex_id_1) <- vertex_attr(G, "number", index = vertex_id_1) + 1
    #
    #             # set number attribute of vertex j to 0 as it has been merged and does not need to be checked again
    #             vertex_attr(G, "number", index = vertex_id_2) <- 0
    #
    #             #append id of vertex j to the list of peptides to be deleted
    #             list_length <- length(duplicated_peptides)
    #             duplicated_peptides[[(list_length + 1)]] <- vertex_id_2
    #           }
    #         }
    #       }
    #     }
    #   }
    # }

    duplicated_proteins <- list() #list of protein ids to be deleted due to merging

    if(length(protein_ids) > 1){
      for(i in 1:(length(protein_ids) - 1)){ # for all peptides i until the penultimate
        for(j in (i+1):length(protein_ids)){ # for all peptides j>i
          vertex_id_1 <- protein_ids[[i]] #id of protein i
          vertex_id_2 <- protein_ids[[j]] #id of protein j

          # if the vertex j has not been marked for deletion yet
          if(vertex_attr(G, "number", index = vertex_id_2) != 0){

            # if vertices i and j have the same amount of edges
            if(length(G[[vertex_id_1]][[1]]) == length(G[[vertex_id_2]][[1]])){

              # if all edges of i and j are identical
              if(!FALSE %in% unlist(unname(G[[vertex_id_1]][[1]] == G[[vertex_id_2]][[1]]))){

                subgraph_list_length <- length(subgraphs_with_merged_proteins)
                subgraphs_with_merged_proteins[[subgraph_list_length + 1]] <- k

                # raise number attribute of vertex i by 1 as it will be merged with vertex j
                vertex_attr(G, "number", index = vertex_id_1) <- vertex_attr(G, "number", index = vertex_id_1) + 1

                #append name of second protein to first protein's node

                vertex_attr(G, "name", index = vertex_id_1) <- paste(vertex_attr(G, "name", index = vertex_id_1),  vertex_attr(G, "name", index = vertex_id_2), sep = ", ")

                # set number attribute of vertex j to 0 as it has been merged and does not need to be checked again
                vertex_attr(G, "number", index = vertex_id_2) <- 0

                #append id of vertex j to the list of proteins to be deleted
                list_length <- length(duplicated_proteins)
                duplicated_proteins[[(list_length + 1)]] <- vertex_id_2
              }
            }
          }
        }
      }
    }

    duplicated_vertices <- c(duplicated_proteins)

    #build new graph with merged peptides and proteins
    G2 <- delete_vertices(G, duplicated_vertices)


    setpb(pb, (k))
    merged_subgraphs_tmp[[k]] <- G2
  }
  # 1min30
  subgraphs_merged_proteins[[l]] <- merged_subgraphs_tmp
}


saveRDS(subgraphs_merged_proteins, "data/D4_without_isoforms/D4_quant/preprocessed/Subgraphs_with_pep_ratios_mergedProteins.rds")


### TODO: Namen der Liste = Comparisons! Das fehlt hier noch!




################################################################################
#### Fasse Protein und Peptid-Knoten zusammen

library(pbapply)

subgraphs <- readRDS("data/D4_without_isoforms/D4_quant/preprocessed/Subgraphs_with_pep_ratios.rds")

subgraphs_merged_proteins <- list()


for (l in 1:21) {

  print(l)


  subgraphs_tmp <- subgraphs[[l]]

  number_of_graphs <- length(subgraphs_tmp)

  subgraphs_with_merged_proteins <- list()
  merged_subgraphs <- vector(mode = 'list', length = number_of_graphs)

  #add progress bar to loop
  number_of_iterations <- number_of_graphs
  pb <- startpb(0, number_of_graphs)
  on.exit(closepb(pb))


  merged_subgraphs_tmp <- vector(mode = 'list', length = number_of_graphs)

  for(k in 1:(number_of_graphs)){

    #print(k)


    G <- subgraphs_tmp[[k]]

    peptide_ids <- as.numeric(V(G)[!V(G)$type])  # ids of peptide nodes (type = FALSE)
    protein_ids <- as.numeric(V(G)[V(G)$type])   # ids of protein nodes (type = TRUE)

    #assign initial number attribute to all vertices

    for(i in 1:gorder(G)){
      vertex_attr(G, "number", index = i) <- 1
    }

    duplicated_peptides <- list() #list of peptide ids to be deleted due to merging

    if(length(peptide_ids) > 1){
      for(i in 1:(length(peptide_ids) - 1)){ # for all peptides i until the penultimate
        for(j in (i+1):length(peptide_ids)){ # for all peptides j>i
          vertex_id_1 <- peptide_ids[[i]] #id of peptide i
          vertex_id_2 <- peptide_ids[[j]] #id of peptide j

          # if the vertex j has not been marked for deletion yet
          if(vertex_attr(G, "number", index = vertex_id_2) != 0){

            # if vertices i and j have the same amount of edges
            if(length(G[[vertex_id_1]][[1]]) == length(G[[vertex_id_2]][[1]])){

              # if all edges of i and j are identical
              if(!FALSE %in% unlist(unname(G[[vertex_id_1]][[1]] == G[[vertex_id_2]][[1]]))){

                # raise number attribute of vertex i by 1 as it will be merged with vertex j
                vertex_attr(G, "number", index = vertex_id_1) <- vertex_attr(G, "number", index = vertex_id_1) + 1

                # set number attribute of vertex j to 0 as it has been merged and does not need to be checked again
                vertex_attr(G, "number", index = vertex_id_2) <- 0

                #append id of vertex j to the list of peptides to be deleted
                list_length <- length(duplicated_peptides)
                duplicated_peptides[[(list_length + 1)]] <- vertex_id_2
              }
            }
          }
        }
      }
    }

    duplicated_proteins <- list() #list of protein ids to be deleted due to merging

    if(length(protein_ids) > 1){
      for(i in 1:(length(protein_ids) - 1)){ # for all peptides i until the penultimate
        for(j in (i+1):length(protein_ids)){ # for all peptides j>i
          vertex_id_1 <- protein_ids[[i]] #id of protein i
          vertex_id_2 <- protein_ids[[j]] #id of protein j

          # if the vertex j has not been marked for deletion yet
          if(vertex_attr(G, "number", index = vertex_id_2) != 0){

            # if vertices i and j have the same amount of edges
            if(length(G[[vertex_id_1]][[1]]) == length(G[[vertex_id_2]][[1]])){

              # if all edges of i and j are identical
              if(!FALSE %in% unlist(unname(G[[vertex_id_1]][[1]] == G[[vertex_id_2]][[1]]))){

                subgraph_list_length <- length(subgraphs_with_merged_proteins)
                subgraphs_with_merged_proteins[[subgraph_list_length + 1]] <- k

                # raise number attribute of vertex i by 1 as it will be merged with vertex j
                vertex_attr(G, "number", index = vertex_id_1) <- vertex_attr(G, "number", index = vertex_id_1) + 1

                #append name of second protein to first protein's node

                vertex_attr(G, "name", index = vertex_id_1) <- paste(vertex_attr(G, "name", index = vertex_id_1),  vertex_attr(G, "name", index = vertex_id_2), sep = ", ")

                # set number attribute of vertex j to 0 as it has been merged and does not need to be checked again
                vertex_attr(G, "number", index = vertex_id_2) <- 0

                #append id of vertex j to the list of proteins to be deleted
                list_length <- length(duplicated_proteins)
                duplicated_proteins[[(list_length + 1)]] <- vertex_id_2
              }
            }
          }
        }
      }
    }

    duplicated_vertices <- c(duplicated_proteins)

    #build new graph with merged peptides and proteins
    G2 <- delete_vertices(G, duplicated_vertices)


    setpb(pb, (k))
    merged_subgraphs_tmp[[k]] <- G2
  }
  # 1min30
  subgraphs_merged_proteins[[l]] <- merged_subgraphs_tmp
}


saveRDS(subgraphs_merged_proteins, "data/D4_without_isoforms/D4_quant/preprocessed/Subgraphs_with_pep_ratios_mergedProteins_and_Peptides.rds")


### TODO: Namen der Liste = Comparisons! Das fehlt hier noch!




















################################################################################
################################################################################
################################################################################
#### calculate isomorph lists separately per comparison

library(igraph)
source("R Scripts/bipartite_graph_generation_and_analysis/helper_functions/isomorph_classes_calculation_and_plotting_functions.R")

S <- readRDS("data/D4_without_isoforms/D4_quant/preprocessed/Subgraphs_with_pep_ratios_mergedProteins.rds")

S2 <- readRDS("data/D4_without_isoforms/D4_quant/preprocessed/Subgraphs_with_pep_ratios_mergedProteins.rds")
names(S) <- names(S2)

isomorph_list <- list()

for (i in 1:21) {
  print(i)
  S_tmp <- S[[i]]
  isomorph_list[[i]] <- calculateIsomorphList(S_tmp, matrix = FALSE)
}

names(isomorph_list) <- names(S)

save(isomorph_list, file = "data/D4_without_isoforms/D4_quant/isomorph_classes/isomorph_classes_merged_Proteins_D4.RData")




################################################################################
#### calculate isomorph lists together over all comparisons

lengths(S)
# [1] 4498 4424 3545 3594 3114 2812 4622 3796 3825 3347 3040 4263 4436 3997 3734 4351 3940 3665 4803 4549 4595
S_all <- unlist(S, recursive = FALSE)

isomorph_all <- calculateIsomorphList(S_all, matrix = FALSE)
save(isomorph_all, file = "data/D4_without_isoforms/D4_quant/isomorph_classes/isomorph_classes_all.RData")






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

