### aus Studienprojekt WS 22/23
#' Collapsing of peptide and protein nodes of an edgelist.
#'
#' @param edgelist edgelist
#' @param collapse_protein_nodes
#' @param collapse_peptide_nodes
#'
#' @return Edgelist with collapsed protein and peptide nodes
#' @export
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' digested_proteins <- bppg::digest_fasta(fasta)
#' edgelist <- bppg::generate_edgelist(digested_proteins)
#' edgelist_collapsed <- bppg::collapse_edgelist(edgelist)
#'

collapse_edgelist <- function(edgelist,
                              collapse_protein_nodes = TRUE,
                              collapse_peptide_nodes = TRUE) {

  if (!collapse_protein_nodes & !collapse_peptide_nodes) {
    return(edgelist)
  }


  if (collapse_protein_nodes) {
    # function(x) paste(sort(unique(x)), collapse = ";")
    ### aggregate peptide sequences that belong to the same protein accession (1 row per protein accession)
    protEdges <- aggregate(data = edgelist, x = peptide ~ protein, function(x) paste(x, collapse = ";"))
    ### aggregate proteins with the same set of peptides (-> protein nodes)
    protNodes <- aggregate(data = protEdges, x = protein ~ peptide, function(x) paste(sort(unique(x)), collapse = ";"))
  } else {
    protEdges <- aggregate(data = edgelist, x = peptide ~ protein, function(x) paste(sort(unique(x)), collapse = ";"))
    protNodes <- protEdges
  }


  if (collapse_peptide_nodes) {
    ### aggregate protein accessions belonging to the same peptide sequences (1 row per peptide sequence)
    pepEdges <- aggregate(data = edgelist, x = protein ~ peptide, function(x) paste(sort(unique(x)), collapse = ";"))
    ### aggregate peptides with the same set of proteins (-> peptide nodes)
    pepNodes <- aggregate(data = pepEdges, x = peptide ~ protein, function(x) paste(sort(unique(x)), collapse = ";"))
  } else {
    pepEdges <- aggregate(data = edgelist, x = protein ~ peptide, function(x) paste(sort(unique(x)), collapse = ";"))
    pepNodes <- pepEdges
  }


  edgelist2 <- edgelist
  #keep <- logical(nrow(edgelist2))

  pepNodes2 <- pepNodes
  pepNodes2$peptide <- limma::strsplit2(pepNodes2$peptide, ";")[,1]  # erstes Peptid aus Liste!
  edgelist2 <- edgelist[edgelist$peptide %in% pepNodes2$peptide,]

  protNodes2 <- protNodes
  protNodes2$protein <- limma::strsplit2(protNodes2$protein, ";")[,1]  # erstes Peptid aus Liste!
  edgelist3 <- edgelist2[edgelist2$protein %in% protNodes2$protein,]

  edgelist4 <- edgelist3
  edgelist4$protein <- protNodes$protein[match(edgelist3$protein, protNodes2$protein)]
  edgelist4$peptide <- pepNodes$peptide[match(edgelist3$peptide, pepNodes2$peptide)]

  invisible(NULL)

  return(edgelist4)
}


