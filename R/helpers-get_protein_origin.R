

#' For each protein node, get information about protein origin (contaminant, spike-in
#' or a specific species).
#'
#' @param accessions character vector with protein accessions (inside a protein node), separated by ";"
#' @param contaminants character vector of contaminant protein accessions
#' @param spike_ins character vector of spike-in protein accessions
#' @param organisms names list with character vectors of protein accessions for each organism
#'
#' @return character vector with protein origin for each protein accession
#' @export
#'
#' @examples # TODO
get_protein_origin <- function(accessions,
                               contaminants = NULL,
                               spike_ins = NULL,
                               organisms = NULL) {


  ## x = string with accessions, separated by ";"
  get_origin <- function(x, contaminants, spike_ins, organisms) {

    accessions <- unlist(strsplit(x, ";"))

    origin_tmp <- character(length(accessions))
    for (i in 1:length(accessions))  {
      if (accessions[i] %in% contaminants) {
        origin_tmp[i] <- "Contaminant"
        next
      } else if (accessions[i] %in% spike_ins) {
        origin_tmp[i] <- "Spike-in"
        next
      } else {
        for (org in names(organisms)) {
          if (accessions[i] %in% organisms[[org]]) {
            origin_tmp[i] <- org
            next
          }
        }
      }
    }


    ### if there are multiple accessions, the origin has to be summarized
    if (length(accessions) > 1) {
      if (length(unique(accessions)) == 1) {
        origin <- origin_tmp[1]
      } else {
        origin <- paste0(sort(unique(origin_tmp)), collapse = ";")
      }
    } else {
      origin <- origin_tmp
    }

    return(origin)
  }

  origin <- pbapply::pbsapply(accessions, get_origin, contaminants = contaminants,
                              spike_ins = spike_ins, organisms = organisms)
  return(origin)
}
