
#' Tries to recursively remove protein nodes from graph while keeping the error term in the optimization step low
#'
#' @param G igraph object
#' @param threshold threshold for increase of error term, default 1.05 refers to 5% increase
#' @param iter iteration (this is a recursive function)
#' @param min_error_ref minimal error term with all available protein nodes
#' @param min_error_current current minimal error term
#' @param protein_nodes_list current list of protein nodes
#' @param combination_list list of node combinations
#' @param error_list list of error terms
#' @param comb_current current combination of protein nodes
#' @param G_current current graph (with removed protein nodes)
#' @param n_comb_current current number of protein nodes
#'
#' @return list
#' @export
#'
#' @examples # TODO
protein_elimination <- function(G,
                                threshold = 1.05,
                                iter = 0,
                                min_error_ref = NULL,
                                min_error_current = NULL,
                                protein_nodes_list = igraph::V(G)[igraph::V(G)$type],
                                combination_list = NULL,
                                error_list = NULL,
                                comb_current = NULL,
                                G_current = NULL,
                                n_comb_current = NA) {

  #### 0: calculate reference error
  if (iter == 0) {
    G <- bppg::add_uniqueness_attributes(G)
    proteinnodes <- igraph::V(G)[igraph::V(G)$type]
    nr_unique_peptides <- igraph::V(G)$nr_unique_peptides[igraph::V(G)$type]

    X <- igraph::as_incidence_matrix(G)
    fc <- stats::na.omit(igraph::V(G)$pep_ratio)
    S <- list(X = X, fc = fc)
    n <- ncol(S$X) # number of protein nodes

    opti <- minimize_squared_error(S, #error.type = "multiplicative",
                                   fixed.Ci = NULL,
                                   verbose = FALSE, #error.trans = "square",
                                   reciprocal = FALSE, log_level = TRUE, control = list(trace = 0, delta = 1e-9))
    min_error_ref <- opti$RES$res_squ_err

    # initialization
    combination_list <- paste(proteinnodes, collapse = ",")
    comb_current <- paste(proteinnodes, collapse = ",")
    G_current <- G
    error_list <- min_error_ref
    n_comb_current <- length(proteinnodes)
  }


  ### try to remove every protein node, if the error is sill small enough, try to remove the next protein
  proteinnodes <- igraph::V(G)[igraph::V(G)$type]
  nr_unique_peptides <- igraph::V(G)$nr_unique_peptides[igraph::V(G)$type]

  for (i in 1:length(proteinnodes)) {

    proteinnodes_tmp <- protein_nodes_list[-i]
    combination <- paste(proteinnodes_tmp, collapse = ",")

    ### this combination has already been tested
    if (combination %in% combination_list) {
      next
    }

    combination_list <- c(combination_list, combination)

    error_list <- c(error_list, NA)


    # skip if protein has unique peptides
    if (nr_unique_peptides[i] > 0) {
      next
    }

    # skip if it has more proteins than the current combination
    if (length(proteinnodes_tmp) > n_comb_current) {
      next
    }


    # this will delete a protein node and all associated edges
    G_tmp <- igraph::delete_vertices(G, proteinnodes[i])
    # different results are possible
    # 1) the graph is still connected and all peptide nodes are still covered
    # 2) at least one peptide node is not connected anymore (this has to be skipped then!)
    # 3) all peptide nodes are covered but the graph is not connected anymore (this has to be skipped then!)

    # check if all peptide nodes are still connected to at least one protein node
    if (any(igraph::ego_size(G_tmp, order = 1, nodes = igraph::V(G_tmp)[!igraph::V(G_tmp)$type], mindist = 1) == 0)) {
      next
    }

    ## decompose into connected components if possible
    G_CC <- igraph::decompose(G_tmp)

    min_error_tmp <- 0
    for (i in 1:length(G_CC)) {

      X <- igraph::as_incidence_matrix(G_CC[[i]])
      fc <- stats::na.omit(igraph::V(G_CC[[i]])$pep_ratio)
      S <- list(X = X, fc = fc)
      n <- ncol(S$X) # number of protein groups


      opti <- minimize_squared_error(S, #error.type = "multiplicative",
                                     fixed.Ci = NULL,
                                     verbose = FALSE, #error.trans = "square",
                                     reciprocal = FALSE, log_level = TRUE, control = list(trace = 0, delta = 1e-9))
      min_error_tmp <- min_error_tmp + opti$RES$res_squ_err
    }
    error_list[length(error_list)] <- min_error_tmp

    ## falls min_error_tmp NA sein sollte, skippen (seltene Fälle wenn ein Ci 0 ist)
    if (is.na(min_error_tmp)) {
      next
    }

    ## check if error is small enough compared to reference error
    if (min_error_tmp > min_error_ref * threshold) {
      next
    }

    min_error_current <- min_error_tmp
    comb_current <- combination
    G_current <- G_CC
    n_comb_current <- length(proteinnodes_tmp)


    RES <- protein_elimination(G = G_tmp, threshold = threshold, iter = 1, min_error_ref = min_error_ref,
                                 min_error_current = min_error_current, protein_nodes_list = proteinnodes_tmp, error_list = error_list,
                                 combination_list = combination_list, comb_current = comb_current, G_current = G_current,
                                 n_comb_current = n_comb_current)
    min_error_current <- RES$min_error_current
    combination_list <- RES$combination_list
    error_list <- RES$error_list
    comb_current <- RES$comb_current
    G_current <- RES$G_current
    n_comb_current <- RES$n_comb_current
  }



  return(list(min_error_ref = min_error_ref, min_error_current = min_error_current, protein_nodes_list = protein_nodes_list,
              combination_list = combination_list, error_list = error_list,
              comb_current = comb_current, G_current = G_current, n_comb_current = n_comb_current))
}





