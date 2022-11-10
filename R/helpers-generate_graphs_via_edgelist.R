

#' calculate edgelist from a list of digested proteins
#'
#' @param digested_proteins Output from digest_fasta (List of vectors of peptide sequences)
#'
#' @return edgelist
#' @export
#'
#' @examples
#' library(seqinr)
#' #file <- system.file("extdata", "2020_01_31_proteome_S_cerevisae.fasta", package = "bppg")
#' #fasta <- seqinr::read.fasta(file = file)
#' #digested_proteins <- digest_fasta(fasta)
#' #res <- calculate_edgelist(digested_proteins)
#'
calculate_edgelist <- function(digested_proteins) {

  #calculate necessary number of edges by counting the peptides belonging to each protein
  mat_length <- sum(lengths(digested_proteins))

  #generate empty edge matrix of size (#edges)x2
  edgelist <- matrix(nrow = mat_length, ncol = 2)

  #add progress bar to loop
  number_of_iterations <- length(digested_proteins)
  pb <- pbapply::startpb(0, length(digested_proteins))
  on.exit(pbapply::closepb(pb))

  #add an entry to the edge matrix for each peptide-protein relation in the digested_proteins matrix
  current_row <- 1
  for (i in 1:length(digested_proteins)){
    if(length(digested_proteins[[i]]) != 0){
      for (j in 1:length(digested_proteins[[i]])){
        edgelist[current_row, 1] <- names(digested_proteins)[[i]]
        edgelist[current_row, 2] <- digested_proteins[[i]][[j]]
        current_row <- current_row + 1
      }
      pbapply::setpb(pb, i)
    }
  }

  #progress bar command
  invisible(NULL)

  #find and remove duplicate rows that would lead to duplicate edges
  duplicate_rows <- duplicated(edgelist, margin = 1)
  edgelist <- edgelist[-which(duplicate_rows),]

  return(edgelist)

}




#filter_edgelist <- function(edgelist, )




#' Generate bipartite peptide-protein graphs from a list of digested proteins via an edgelist
#'
#' @param edgelist Output from calculate_edgelist
#'
#' @return List of subgraphs as igraph objects.
#' @export
#'
#' @examples
#' library(seqinr)
#' #file <- system.file("extdata", "2020_01_31_proteome_S_cerevisae.fasta", package = "bppg")
#' #fasta <- seqinr::read.fasta(file = file)
#' #digested_proteins <- digest_fasta(fasta)
#' #edgelist <- calculate_edgelist(digested_proteins)
#' #res <- generate_graphs_via_edgelist(digested_proteins)
#'
generate_graphs_from_edgelist <- function(edgelist) {

  #generate graph from edge matrix
  G <- igraph::graph_from_edgelist(edgelist, directed = FALSE)

  #assign vertex types to proteins and peptides for the graph to be bipartite
  igraph::V(G)[igraph::V(G)$name %in% names(edgelist[,1])]$type <- TRUE
  igraph::V(G)[igraph::V(G)$name %in% names(edgelist[,2])]$type <- FALSE
  ### TODO: export G

  #decompose graph into connected components
  subgraphs <- igraph::decompose(G)
  return(subgraphs)

}
