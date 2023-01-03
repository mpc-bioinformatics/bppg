
#' Assign protein accessions to a list of peptides, depending on a FASTA file.
#'
#' @param sequence vector of peptide sequences
#' @param fasta_vec names vector of protein sequences from fasta file(s)
#'
#' @return list of assigned proteins
#' @export
#'
#' @examples
#' ### TODO
assign_protein_accessions <- function(sequence, fasta_vec) {

  protein_accessions <- names(fasta_vec)

  assigned_proteins <- pbapply::pblapply(sequence, function(x) {
    ind <- which(grepl(x, fasta_vec))
    proteins <- protein_accessions[ind]

    ### test if protein is really tryptic
    is_tryptic <- rep(TRUE, length(ind))
    for (i in 1:length(ind)) {
      dig <- Digest2(fasta_vec[ind[i]], missed = 2)
      if (!(x %in% dig$sequence)) {
        if(any(dig$start == 1)) {
          dig2 <- dig[dig$start == 1, ]  # only peptides at protein N-terminus, where initial M could possible have been cut off
          if (!(paste0("M", x) %in% dig2$sequence)) {
            is_tryptic[i] <- FALSE
            next()
          }
        } else {
          is_tryptic[i] <- FALSE
        }
      }
    }
    proteins <- proteins[is_tryptic]
    proteins <- sort(proteins)

    proteins <- BBmisc::collapse(proteins, sep = "/")
    return(proteins)
  })

  return(unlist(assigned_proteins))
}

### TODO: wird abgelöst durch Grapherstellung über die Edgelist (Studienprojekt WS 22/23)

# sequence, fasta_vec

# protein_accessions <- names(fasta_vec)
#
# ### 1000 -> 12min
#
# x <- pblapply(sequence, function(x) {
#   ind <- grep(paste0("(?:^M|K|R|^)(", x, ")(?=[^P]|$)"), fasta_vec, perl = TRUE)
#   ### TODO: passt noch nicht ganz 100%ig, weil
#   #proteins <- protein_accessions[ind]
#   #proteins <- BBmisc::collapse(proteins, sep = "/")
#   return(ind)
# }
# )
# #
# #
# #
# # table(x[[1]])
#
# # (?:^M|K|R|^)(WLSPEEVL)(?=[^P]|$)
#
# sequence[1:10]
# fasta_vec[x[[1]]]
