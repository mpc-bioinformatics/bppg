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
generate_quant_graphs <- function(peptide_ratios, id_cols = 1, fasta_edgelist, outpath = NULL, seq_column = "Sequence",
                                  collapse_protein_nodes = TRUE, collapse_peptide_nodes = FALSE, suffix = "") {

  ### broad filtering for edgelist for only quantifies peptides
  edgelist_filtered <- fasta_edgelist[fasta_edgelist[,2] %in% peptide_ratios[,seq_column], ]


  if (!is.null(outpath)) {
    openxlsx::write.xlsx(edgelist_filtered, file = paste0(outpath, "edgelist_filtered_", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)
  }


  id <- peptide_ratios[, id_cols, drop = FALSE]
  peptide_ratios <- peptide_ratios[, -(id_cols), drop = FALSE]

  colnames_split <- limma::strsplit2(colnames(peptide_ratios), "_")
  comparisons <- paste(colnames_split[,2], colnames_split[,3], sep = "_")

  subgraphs <- list()

  for (i in 1:ncol(peptide_ratios)) {
    comparison <- comparisons[i]

    fc <- peptide_ratios[,i]
    peptides_tmp <- id[,seq_column][!is.na(fc)]  ## peptides that are quantified in this specific comparison
    fc <- stats::na.omit(fc)
    edgelist_filtered2 <- edgelist_filtered[edgelist_filtered[,2] %in% peptides_tmp, ]


    ## add peptide ratios
    edgelist_filtered2$pep_ratio <- peptide_ratios[,i][match(edgelist_filtered2$peptide, id[, seq_column])]


    ## generate whole bipartite graph
    edgelist_coll <- bppg::collapse_edgelist_quant(edgelist_filtered2,
                                             collapse_protein_nodes = collapse_protein_nodes,
                                             collapse_peptide_nodes = collapse_peptide_nodes)

    G <- bppg::generate_graphs_from_edgelist(edgelist_coll[, 1:2])

    ### set peptide ratios as vertex attributes
    for (j in 1:length(G)){

      G[[j]] <- igraph::set_vertex_attr(graph = G[[j]], name = "pep_ratio",
                                        index = igraph::V(G[[j]])[!igraph::V(G[[j]])$type],
                                        value = edgelist_coll$pep_ratio[match(igraph::V(G[[j]])$name[!igraph::V(G[[j]])$type], edgelist_coll$peptide)])
    }

    subgraphs[[i]] <- G
    names(subgraphs)[[i]] <- comparison

  }

  return(subgraphs)

}


#' Generate graphs from quantitative peptide-level data
#'
#' @param D data set with peptide sequence as first column and peptide intensities in subsequent columns
#' (e.g. output from bppg::read_MQ_peptidetable)
#' @param fasta fasta file used for identification of peptides in D
#' @param outpath bla
#' @param normalize currently only loess normalization possible
#' @param missed_cleavages bla
#' @param min_aa bla
#' @param max_aa bla
#' @param ... currently not in use
#'
#' @return list of list of graphs
#' @export
#'
#' @examples
generate_graphs_from_quant_data <- function(D, fasta, outpath = NULL, normalize = FALSE,
                                            missed_cleavages = 2, min_aa = 6, max_aa = 50,
                                            id_columns = 1, seq_column = "Sequence",
                                            collapse_protein_nodes = TRUE, collapse_peptide_nodes = FALSE,
                                            suffix = "",
                                            ...) {

  message("Digesting FASTA file...")
  digested_proteins <- bppg::digest_fasta(fasta, missed_cleavages = missed_cleavages,
                                          min_aa = min_aa, max_aa = max_aa)#, ...)
  message("Generating edgelist ...")
  edgelist <- bppg::generate_edgelist(digested_proteins)

  if (!is.null(outpath)) {
    openxlsx::write.xlsx(edgelist, file = paste0(outpath, "edgelist_fasta_", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)
  }


  # remove peptides outside the desired length range
  D <- D[nchar(D[, seq_column]) >= min_aa & nchar(D[, seq_column]) <= max_aa,]


  #normalize Intensities
  intensities <- D[,-id_columns]

  ### aggregate replicates by calculating the mean
  group <- factor(limma::strsplit2(colnames(intensities), "_")[,1])
  D_aggr <- bppg::aggregate_replicates(D, method = "mean", missing.limit = 0.4,
                                       group = group, id_cols = id_columns)
  if (!is.null(outpath)) {
    openxlsx::write.xlsx(D_aggr, file = paste0(outpath, "aggr_peptides_", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)
  }

  ### calculate the peptide ratio table
  groups  <- levels(group)
  peptide_ratios <- bppg::calculate_peptide_ratios(aggr_intensities = D_aggr, id_cols = id_columns, group_levels = groups)
  if (!is.null(outpath)) {
    openxlsx::write.xlsx(peptide_ratios, file = paste0(outpath, "peptide_ratios_", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)
  }


  ## Generierung der Graphen (man braucht peptide_ratios und fast_edgelist!)
  graphs <- bppg::generate_quant_graphs(peptide_ratios = peptide_ratios, id_cols = id_columns, fasta_edgelist = edgelist,
                                        outpath = outpath, seq_column = seq_column,
                                        collapse_protein_nodes = collapse_protein_nodes,
                                        collapse_peptide_nodes = collapse_peptide_nodes,
                                        suffix = suffix)
  return(graphs)

}







