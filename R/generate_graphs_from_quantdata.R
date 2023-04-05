

# D: data frame with peptide-level data (1st column peptide Sequence, then Intensities)
# fasta: fasta file (to assign peptides to proteins)

# fasta1 <- seqinr::read.fasta(file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_fasta/20220722_uniprot_proteome_UP000005640_HomoSapiens_v2022_02.fasta",
#                              seqtype = "AA", as.string = TRUE)
# names(fasta1) <- limma::strsplit2(names(fasta1), "\\|")[,2]
# fasta2 <- seqinr::read.fasta(file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_fasta/20220722_uniprot_proteome_UP000000589_MusMusculus_v2022_02.fasta",
#                              seqtype = "AA", as.string = TRUE)
# names(fasta2) <- limma::strsplit2(names(fasta2), "\\|")[,2]
# fasta3 <- seqinr::read.fasta(file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_fasta/MaxQuant_contaminants_v2_3_1.fasta",
#                              seqtype = "AA", as.string = TRUE)
# names(fasta3) <- paste0("CON_", names(fasta3))
# fasta <- c(fasta1, fasta2, fasta3)
#
# DATA <- bppg::read_MQ_peptidetable("C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/peptides.txt",
#                                    LFQ = TRUE, remove_contaminants = FALSE)
#
# DATA2 <- bppg::read_MQ_peptidetable("C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/peptides.txt",
#                                    LFQ = FALSE, remove_contaminants = FALSE)
#


### TODO: save end and intermediate results

#' Generate graphs from quantitative peptide-level data
#'
#' @param D data set with peptide sequence as first column and peptide intensities in subsequent columns
#' (e.g. output from bppg::read_MQ_peptidetable)
#' @param fasta fasta file used for identification of peptides in D
#'
#' @return list of list of graphs
#' @export
#'
#' @examples
generate_graphs_from_quant_data <- function(D, fasta, ...) {
  message("Digesting FASTA file...")
  digested_proteins <- bppg::digest_fasta(fasta, ...)#, ...)
  message("Generating edgelist ...")
  edgelist <- bppg::generate_edgelist(digested_proteins)
  #### TODO: Möglichkeit, eine vorberechnete edgelist anzugeben!
  #edgelist_filtered <- edgelist[edgelist[,2] %in% D$Sequence, ]


  #prepare for normalization
  #id <- D[,1]
  intensities <- D[,-1]
  #normalize Intensities
  # TODO: auch Median, Quantilsnormalisierung und LFQ-Normalisierung von MaxQuant erlauben?
  intensities <- 2^limma::normalizeBetweenArrays(log2(intensities), method = "cyclicloess")


  #D <- read.table(file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/preprocessed_peptide_data_D4.txt", header = TRUE)
  group <- factor(limma::strsplit2(colnames(intensities), "_")[,1])
  D_aggr <- bppg::aggregate_replicates(D, method = "mean", missing.limit = 0.4,
                                       group = group, id_cols = 1)
  #write.table(D_aggr, file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/aggr_pept_int_D4.txt",
  #            row.names = FALSE, sep = "\t")
  #save(group, file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/group.RData")
#  D_aggr <- read.table("C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/aggr_pept_int_D4.txt",
 #                        sep = "\t", header = TRUE)
  groups  <- levels(group)

  peptide_ratios <- bppg::calculate_peptide_ratios(aggr_intensities = D_aggr, id_cols = 1, group_levels = groups)
  #write.table(peptide_ratios, file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/D4_without_isoforms/D4_quant/preprocessed/peptide_ratios.txt",
  #            row.names = FALSE, sep = "\t")

  ## generierung der Graphen (man braucht peptide_ratios und fast_edgelist!)

  graphs <- bppg::generate_quant_graphs(peptide_ratios = peptide_ratios, id_cols = 1, fasta_edgelist = edgelist)
  return(graphs)

}




#' Generate graphs from peptide ratio table, using an edgelist calculated on the fasta file
#'
#' @param peptide_ratios table with peptide ratios
#' @param id_cols columns with ids, e.g. peptide sequences (everything except the peptide ratios)
#' @param fasta_edgelist Edgelist created from the corresponding FASTA file
#'
#' @return list of list of subgraphs
#' @export
#'
#' @examples
generate_quant_graphs <- function(peptide_ratios, id_cols = 1, fasta_edgelist) {

  ### broad filtering for edgelist for only quantifies peptides
  edgelist_filtered <- fasta_edgelist[fasta_edgelist[,2] %in% peptide_ratios$peptides, ]

  id <- peptide_ratios[, id_cols, drop = FALSE]
  peptide_ratios <- peptide_ratios[, -(id_cols), drop = FALSE]

  colnames_split <- limma::strsplit2(colnames(peptide_ratios), "_")
  comparisons <- paste(colnames_split[,2], colnames_split[,3], sep = "_")


 # graphs <- list()
  subgraphs <- list()
  #i = 3

  for (i in 1:ncol(peptide_ratios)) {

    comparison <- comparisons[i]
      #substr(colnames(peptide_ratios)[i], 4 , 100)
    fc <- peptide_ratios[,i]
    peptides_tmp <- id$peptides[!is.na(fc)]  ## peptides that are quantified in this specific comparison
    fc <- stats::na.omit(fc)
    edgelist_filtered2 <- edgelist_filtered[edgelist_filtered[,2] %in% peptides_tmp, ]

    ## generate whole bipartite graph
    G_tmp <- igraph::graph_from_edgelist(as.matrix(edgelist_filtered2[, 1:2]), directed = FALSE)
    igraph::V(G_tmp)[igraph::V(G_tmp)$name %in% edgelist_filtered2[,1]]$type <- TRUE
    igraph::V(G_tmp)[igraph::V(G_tmp)$name %in% edgelist_filtered2[,2]]$type <- FALSE

    ## set peptide ratios as vertex attributes
    igraph::vertex_attr(G_tmp, "pep_ratio", index = igraph::V(G_tmp)[!igraph::V(G_tmp)$type]) <- fc[match(igraph::V(G_tmp)$name[!igraph::V(G_tmp)$type], peptides_tmp)]

    ## calculate connected components
    subgraphs_tmp <- igraph::decompose(G_tmp)

    #   graphs[[i-2]] <- G_tmp
    #   names(graphs)[[i-2]] <- comparison

    subgraphs[[i]] <- subgraphs_tmp
    names(subgraphs)[[i]] <- comparison

  }

  return(subgraphs)

}






