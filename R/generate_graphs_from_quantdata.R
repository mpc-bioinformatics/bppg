#' Generate graphs from peptide ratio table, using an edgelist calculated on the fasta file.
#'
#' @param peptide_ratios           \strong{data.frame} \cr
#'                                 A table with peptide ratios.
#' @param id_cols                  \strong{integer vector} \cr
#'                                 The columns with ids, e.g. peptide sequences (everything except the peptide ratios)
#' @param fasta_edgelist           \strong{data.frame} \cr
#'                                 An edgelist created from the corresponding FASTA file, eg. created with [bppg::generate_edgelist()].
#' @param outpath                  \strong{character} \cr
#'                                 The output path for the results.
#' @param seq_column               \strong{character} \cr
#'                                 The column name of the peptide sequence.
#' @param collapse_protein_nodes   \strong{logical} \cr
#'                                 If \code{TRUE}, the protein nodes will be collapsed.
#' @param collapse_peptide_nodes   \strong{logical} \cr
#'                                 If \code{TRUE}, the peptide nodes will be collapsed.
#' @param suffix                   \strong{character} \cr
#'                                 The suffix for saving results.
#'
#' @return A list of list of subgraphs
#' @export
#'
#' @seealso [bppg::generate_edgelist()]
#'
#' @examples
#'

generate_quant_graphs <- function(peptide_ratios,
                                  id_cols = 1,
                                  fasta_edgelist,
                                  outpath = NULL,
                                  seq_column = "Sequence",
                                  collapse_protein_nodes = TRUE,
                                  collapse_peptide_nodes = FALSE,
                                  suffix = "") {

  # filter out na, leave valid rows only
  peptide_ratios <- stats::na.omit(peptide_ratios)

  # filter edgelist so only valid peptides remain
  edgelist_filtered <- fasta_edgelist[fasta_edgelist[, 2] %in% peptide_ratios[, seq_column], ]


  if (!is.null(outpath)) {
    openxlsx::write.xlsx(edgelist_filtered, file = paste0(outpath, "edgelist_filtered_", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)
  }


  id <- peptide_ratios[, id_cols, drop = FALSE]
  fc <- peptide_ratios[, -(id_cols), drop = FALSE]

  colnames_split <- limma::strsplit2(colnames(fc)[1], "_")
  comparison <- paste(colnames_split[, 2], colnames_split[, 3], sep = "_")

  subgraphs <- list()

  ## add peptide ratios
  edgelist_filtered$pep_ratio <- fc[match(edgelist_filtered$peptide, id[, seq_column]), 1]
  edgelist_filtered$imputed <- fc[match(edgelist_filtered$peptide, id[, seq_column]), 2]


  ## generate whole bipartite graph
  edgelist_coll <- bppg::collapse_edgelist_quant(edgelist_filtered,
                                            collapse_protein_nodes = collapse_protein_nodes,
                                            collapse_peptide_nodes = collapse_peptide_nodes)

  edgelist_coll_pep <- bppg::collapse_edgelist_quant(edgelist_filtered,
                                            collapse_protein_nodes = collapse_protein_nodes,
                                            collapse_peptide_nodes = TRUE)

  # create graphs and return decomposed graph list
  G <- bppg::generate_graphs_from_edgelist(edgelist_coll[, 1:2])

  ### set peptide ratios as vertex attributes
  for (j in 1:length(G)){

    G[[j]] <- igraph::set_vertex_attr(graph = G[[j]], name = "pep_ratio",
                                      index = igraph::V(G[[j]])[!igraph::V(G[[j]])$type],
                                      value = edgelist_coll$pep_ratio[match(igraph::V(G[[j]])$name[!igraph::V(G[[j]])$type], edgelist_coll$peptide)])
    G[[j]] <- igraph::set_vertex_attr(graph = G[[j]], name = "imputed",
                                      index = igraph::V(G[[j]])[!igraph::V(G[[j]])$type],
                                      value = edgelist_coll$imputed[match(igraph::V(G[[j]])$name[!igraph::V(G[[j]])$type], edgelist_coll$peptide)])
  }

  subgraphs <- G
  names(subgraphs) <- comparison

  return(subgraphs)

}


#' Generate graphs from quantitative peptide-level data
#'
#' @param D                        \strong{data.frame} \cr
#'                                 A data set with peptide sequence as first column
#'                                 and peptide intensities in subsequent columns
#'                                 ,e.g. created with [bppg::read_MQ_peptidetable()].
#' @param fasta                    \strong{list of vector of characters} \cr
#'                                 A fasta file used for identification of peptides in D,
#'                                 already read into R by [seqinr::read.fasta()].
#' @param outpath                  \strong{character} \cr
#'                                 The output path for the results.
#' @param missed_cleavages         \strong{integer} \cr
#'                                 The number of allowed missed cleavages in a peptide.
#' @param min_aa                   \strong{integer} \cr
#'                                 The minimum number of amino acids in a peptide.
#' @param max_aa                   \strong{integer} \cr
#'                                 The maximum number of amino acids in a peptide.
#' @param id_columns               \strong{integer vector} \cr
#'                                 The columns of D that contain ID information (the rest should contain only peptide intensities, properly normalized).
#' @param seq_column               \strong{character} \cr
#'                                 The column name of the column with the peptide sequences.
#' @param collapse_protein_nodes   \strong{logical} \cr
#'                                 If \code{TRUE}, the protein nodes will be collapsed.
#' @param collapse_peptide_nodes   \strong{logical} \cr
#'                                 If \code{TRUE}, the peptide nodes will be collapsed.
#' @param suffix                   \strong{character} \cr
#'                                 The suffix for output files.
#' @param ...                      currently not in use
#'
#' @return A list of list of graphs
#' @export
#'
#' @seealso [bppg::read_MQ_peptidetable()], [seqinr::read.fasta()],
#'          [bppg::generate_quant_graphs()], [bppg::generate_graphs_from_FASTA()]
#'
#' @examples
#'

generate_graphs_from_quant_data <- function(D,
                                            fasta,
                                            outpath = NULL,
                                            #normalize = FALSE,
                                            missed_cleavages = 2,
                                            min_aa = 6,
                                            max_aa = 50,
                                            id_columns = 1,
                                            seq_column = "Sequence",
                                            collapse_protein_nodes = TRUE,
                                            collapse_peptide_nodes = FALSE,
                                            suffix = "",
                                            ...) {

  message("Digesting FASTA file...")
  digested_proteins <- bppg::digest_fasta(fasta, missed_cleavages = missed_cleavages,
                                          min_aa = min_aa, max_aa = max_aa)
  message("Generating edgelist ...")
  edgelist <- bppg::generate_edgelist(digested_proteins)

  if (!is.null(outpath)) {
    openxlsx::write.xlsx(edgelist, file = paste0(outpath, "edgelist_fasta_", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)
  }

  # remove peptides outside the desired length range
  D <- D[nchar(D[, seq_column]) >= min_aa & nchar(D[, seq_column]) <= max_aa,]

  intensities <- D[,-id_columns]

  ### aggregate replicates by calculating the mean
  group <- factor(limma::strsplit2(colnames(intensities), "_")[,1])
  D_aggr <- bppg::aggregate_replicates(D, method = "mean", missing.limit = 0.4,
                                       group = group, id_cols = id_columns)
  if (!is.null(outpath)) {
    openxlsx::write.xlsx(D_aggr$agg, file = paste0(outpath, "aggr_peptides_", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)
    openxlsx::write.xlsx(D_aggr$imp, file = paste0(outpath, "aggr_peptides_MV_", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)
  }

  ### calculate the peptide ratio (list of data.frame per comparison)
  groups  <- levels(group)
  peptide_ratios <- bppg::calculate_peptide_ratios(aggr_intensities = D_aggr, id_cols = id_columns, group_levels = groups)
  if (!is.null(outpath)) { 
    # TODO: this will not work with the new data structure, per comparison
    openxlsx::write.xlsx(peptide_ratios, file = paste0(outpath, "peptide_ratios_", suffix, ".xlsx"), overwrite = TRUE, keepNA = TRUE)
  }

  #TODO: loop per comparison, will ich dann list of graphs returnen?
  ## Generierung der Graphen (man braucht peptide_ratios und fast_edgelist!)
  graphs <- bppg::generate_quant_graphs(peptide_ratios = peptide_ratios, id_cols = id_columns, fasta_edgelist = edgelist,
                                        outpath = outpath, seq_column = seq_column,
                                        collapse_protein_nodes = collapse_protein_nodes,
                                        collapse_peptide_nodes = collapse_peptide_nodes,
                                        suffix = suffix)
  return(graphs)

}







