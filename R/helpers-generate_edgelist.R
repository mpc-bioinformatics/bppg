#' Generate edgelist from list of in silico digested proteins.
#'
#' @param digested_proteins   \strong{list of vector of characters} \cr
#'                            The output from [digest_fasta()] (List of vectors of peptide sequences)
#' @param prot_origin         \strong{vector of characters} \cr
#'                            origin of the protein (e.g. organism, spike-in/background etc)
#'
#' @return An edgelist.
#' @export
#'
#' @seealso [digest_fasta()]
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' digested_proteins <- digest_fasta(fasta)
#' edgelist <- generate_edgelist(digested_proteins)
#'
#'

generate_edgelist <- function(digested_proteins, prot_origin = NULL) {
  #calculate necessary number of edges by counting the peptides belonging to each protein
  mat_length <- sum(lengths(digested_proteins))

  #generate empty edge matrix of size (#edges)x2
  if (is.null(prot_origin)) {
    edgelist <- matrix(nrow = mat_length, ncol = 2)
  } else {
    edgelist <- matrix(nrow = mat_length, ncol = 3)
  }


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

        if (!is.null(prot_origin)) {
          edgelist[current_row, 3] <- prot_origin[[i]]
        }
        current_row <- current_row + 1
      }
      pbapply::setpb(pb, i)
    }
  }

  #progress bar command
  invisible(NULL)

  #find and remove duplicate rows that would lead to duplicate edges
  duplicate_rows <- duplicated(edgelist, margin = 1)
  edgelist <- edgelist[!duplicate_rows,]

  edgelist <- as.data.frame(edgelist)
  if(is.null(prot_origin)) {
    colnames(edgelist) <- c("protein", "peptide")
  } else {
    colnames(edgelist) <- c("protein", "peptide", "prot_origin")
  }


  return(edgelist)
}
