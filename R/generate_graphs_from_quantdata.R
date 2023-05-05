

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
#' ### TODO: Einstellbar, ob Peptid-Knoten auch gemergt werden sollen (dann mit geom. Mittel als peptid-ratio).
generate_quant_graphs <- function(peptide_ratios, id_cols = 1, fasta_edgelist, outpath = NULL, seq_column = "Sequence",
                                  collapse_protein_nodes = TRUE, collapse_peptide_nodes = FALSE) {

  ### broad filtering for edgelist for only quantifies peptides

  edgelist_filtered <- fasta_edgelist[fasta_edgelist[,2] %in% peptide_ratios[,seq_column], ]

  if (!is.null(outpath)) {
    openxlsx::write.xlsx(edgelist_filtered, file = paste0(outpath, "edgelist_filtered.xlsx"), overwrite = TRUE, keepNA = TRUE)
  }


  id <- peptide_ratios[, id_cols, drop = FALSE]
  peptide_ratios <- peptide_ratios[, -(id_cols), drop = FALSE]

  colnames_split <- limma::strsplit2(colnames(peptide_ratios), "_")
  comparisons <- paste(colnames_split[,2], colnames_split[,3], sep = "_")


  # graphs <- list()
  subgraphs <- list()
  #i = 3

  ### TODO: progress bar!
  for (i in 1:ncol(peptide_ratios)) { #

    #print(i)

    comparison <- comparisons[i]
    #substr(colnames(peptide_ratios)[i], 4 , 100)
    fc <- peptide_ratios[,i]
    peptides_tmp <- id[,seq_column][!is.na(fc)]  ## peptides that are quantified in this specific comparison
    fc <- stats::na.omit(fc)
    edgelist_filtered2 <- edgelist_filtered[edgelist_filtered[,2] %in% peptides_tmp, ]

    ## generate whole bipartite graph
    #G_tmp <- igraph::graph_from_edgelist(as.matrix(edgelist_filtered2[, 1:2]), directed = FALSE)
    edgelist_coll <- bppg::collapse_edgelist(edgelist_filtered2,
                                             collapse_protein_nodes = TRUE,
                                             collapse_peptide_nodes = FALSE)
    #igraph::V(G_tmp)[igraph::V(G_tmp)$name %in% edgelist_filtered2[,1]]$type <- TRUE
    #igraph::V(G_tmp)[igraph::V(G_tmp)$name %in% edgelist_filtered2[,2]]$type <- FALSE

    G <- bppg::generate_graphs_from_edgelist(edgelist_coll)

    ### set peptide ratios as vertex attributes
    for (j in 1:length(G)){
      #G[[i]] <- set_vertex_attr(G, "pep_ratio", value = pep_ratio)

      #  print(igraph::V(G[[i]])[!igraph::V(G[[i]])$type])
      # print(peptide_ratios[match(igraph::V(G[[i]])$name[!igraph::V(G[[i]])$type], id$Sequence), i])
      #
      # print(igraph::V(G[[i]])$name[!igraph::V(G[[i]])$type])
      ## print(head(id))
      #print(j)

      G[[j]] <- igraph::set_vertex_attr(graph = G[[j]], name = "pep_ratio",
                                        index = igraph::V(G[[j]])[!igraph::V(G[[j]])$type],
                                        value = peptide_ratios[match(igraph::V(G[[j]])$name[!igraph::V(G[[j]])$type], id[, seq_column]), i])

      #igraph::set_vertex_attr(G[[i]], "pep_ratio", index = igraph::V(G[[i]])[!igraph::V(G[[i]])$type]) <- peptide_ratios[match(igraph::V(G[[i]])$name[!igraph::V(G[[i]])$type], id$Sequence), i]
    }

    #  ## set peptide ratios as vertex attributes
    #    igraph::vertex_attr(G_tmp, "pep_ratio", index = igraph::V(G_tmp)[!igraph::V(G_tmp)$type]) <- fc[match(igraph::V(G_tmp)$name[!igraph::V(G_tmp)$type], peptides_tmp)]

    ## calculate connected components
    #   subgraphs_tmp <- igraph::decompose(G_tmp)

    #   graphs[[i-2]] <- G_tmp
    #   names(graphs)[[i-2]] <- comparison

    subgraphs[[i]] <- G
    names(subgraphs)[[i]] <- comparison

  }

  return(subgraphs)

}



### TODO: save end and intermediate results

#' Generate graphs from quantitative peptide-level data
#'
#' @param D data set with peptide sequence as first column and peptide intensities in subsequent columns
#' (e.g. output from bppg::read_MQ_peptidetable)
#' @param fasta fasta file used for identification of peptides in D
#' @param outpath
#' @param normalize currently only loess normalization possible
#' @param missed_cleavages
#' @param min_aa
#' @param max_aa
#' @param ... currently not in use
#'
#' @return list of list of graphs
#' @export
#'
#' @examples # TODO
generate_graphs_from_quant_data <- function(D, fasta, outpath = NULL, normalize = FALSE,
                                            missed_cleavages = 2, min_aa = 6, max_aa = 50,
                                            id_columns = 1, seq_column = "Sequence",
                                            ...) {

  message("Digesting FASTA file...")
  digested_proteins <- bppg::digest_fasta(fasta, missed_cleavages = missed_cleavages,
                                          min_aa = min_aa, max_aa = max_aa)#, ...)
  message("Generating edgelist ...")
  edgelist <- bppg::generate_edgelist(digested_proteins)
  #### TODO: Möglichkeit, eine vorberechnete edgelist anzugeben!
  #edgelist_filtered <- edgelist[edgelist[,2] %in% D$Sequence, ]

  if (!is.null(outpath)) {
    openxlsx::write.xlsx(edgelist, file = paste0(outpath, "edgelist_fasta.xlsx"), overwrite = TRUE, keepNA = TRUE)
  }


  # remove peptides outside the desired length range
  # TODO: remove also peptides with too many missed cleavages
  D <- D[nchar(D[, seq_column]) >= min_aa & nchar(D[, seq_column]) <= max_aa,]


  #prepare for normalization
  #id <- D[,1]
  intensities <- D[,-id_columns]
  #normalize Intensities
  # TODO: auch Median, Quantilsnormalisierung und LFQ-Normalisierung von MaxQuant erlauben?
  if (normalize) {
    intensities <- 2^limma::normalizeBetweenArrays(log2(intensities), method = "cyclicloess")
  }


  ### aggregate replicates by calculating the mean
  group <- factor(limma::strsplit2(colnames(intensities), "_")[,1])
  D_aggr <- bppg::aggregate_replicates(D, method = "mean", missing.limit = 0.4,
                                       group = group, id_cols = id_columns)
  if (!is.null(outpath)) {
    openxlsx::write.xlsx(D_aggr, file = paste0(outpath, "aggr_peptides.xlsx"), overwrite = TRUE, keepNA = TRUE)
  }

  ### calculate the peptide ratios
  groups  <- levels(group)
  peptide_ratios <- bppg::calculate_peptide_ratios(aggr_intensities = D_aggr, id_cols = id_columns, group_levels = groups)
  if (!is.null(outpath)) {
    openxlsx::write.xlsx(peptide_ratios, file = paste0(outpath, "peptide_ratios.xlsx"), overwrite = TRUE, keepNA = TRUE)
  }


  ## Generierung der Graphen (man braucht peptide_ratios und fast_edgelist!)
  graphs <- bppg::generate_quant_graphs(peptide_ratios = peptide_ratios, id_cols = id_columns, fasta_edgelist = edgelist,
                                        outpath = outpath, seq_column = seq_column)
  return(graphs)

}


