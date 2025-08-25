#' Collapsing of peptide and protein nodes of an edgelist.
#'
#' @param edgelist                 \strong{data.frame} \cr
#'                                 An edgelist  with peptide ratios eg. created with [generateEdgelist()].
#' @param collapse_protein_nodes   \strong{logical} \cr
#'                                 If \code{TRUE}, the protein nodes will be collapsed.
#' @param collapse_peptide_nodes   \strong{logical} \cr
#'                                 If \code{TRUE}, the peptide nodes will be collapsed.
#'
#' @return An edgelist with collapsed protein and/or peptide nodes.
#' @export
#'
#' @seealso For edgelists without peptide ratios: [.collapseEdgelist()] \cr
#'          [generateGraphsFromFASTA()], [.generateQuantGraphs()], [generateEdgelist()]
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' digested_proteins <- bppg::digestFASTA(fasta)
#' edgelist <- bppg::generateEdgelist(digested_proteins)
#' edgelist_collapsed <- bppg::.collapseEdgelist(edgelist)
#'


.collapseEdgelistQuant <- function(edgelist,
                                    collapse_protein_nodes = TRUE,
                                    collapse_peptide_nodes = TRUE) {

  if (!collapse_protein_nodes & !collapse_peptide_nodes) {
    return(edgelist)
  }

  ### Calculate list if protein nodes
  if (collapse_protein_nodes) {
    ### aggregate peptide sequences that belong to the same protein accession (1 row per protein accession)
    protEdges <- stats::aggregate(data = edgelist, x = cbind(peptide, pep_ratio) ~ protein, function(x) paste(sort(unique(x)), collapse = ";"))
    ### aggregate proteins with the same set of peptides (-> protein nodes)
    protNodes <- stats::aggregate(data = protEdges, x = protein ~ peptide+pep_ratio, function(x) paste(sort(unique(x)), collapse = ";"))
  } else {
    protEdges <- stats::aggregate(data = edgelist, x = peptide ~ protein, function(x) paste(sort(unique(x)), collapse = ";"))
    protNodes <- protEdges
  }


  ### calculate list of peptide nodes
  if (collapse_peptide_nodes) {
    ### aggregate protein accessions belonging to the same peptide sequences (1 row per peptide sequence)
    pepEdges <- stats::aggregate(data = edgelist, x = protein ~ peptide + pep_ratio, function(x) paste(sort(unique(x)), collapse = ";"))
    ### aggregate peptides with the same set of proteins (-> peptide nodes)
    pepNodes <- stats::aggregate(data = pepEdges, x = cbind(peptide, pep_ratio) ~ protein, function(x) paste(sort(unique(x)), collapse = ";"))
  } else {
    pepEdges <- stats::aggregate(data = edgelist, x = protein ~ peptide + pep_ratio, function(x) paste(sort(unique(x)), collapse = ";"), simplify = FALSE)
    pepNodes <- pepEdges
  }


  edgelist2 <- edgelist

  pepNodes2 <- pepNodes
  pepNodes2$peptide <- limma::strsplit2(pepNodes2$peptide, ";")[,1]  # first peptide from list
  edgelist2 <- edgelist[edgelist$peptide %in% pepNodes2$peptide,]

  protNodes2 <- protNodes
  protNodes2$protein <- limma::strsplit2(protNodes2$protein, ";")[,1]  # first peptide from list
  edgelist3 <- edgelist2[edgelist2$protein %in% protNodes2$protein,]

  edgelist4 <- edgelist3
  edgelist4$protein <- protNodes$protein[match(edgelist3$protein, protNodes2$protein)]
  edgelist4$peptide <- pepNodes$peptide[match(edgelist3$peptide, pepNodes2$peptide)]
  edgelist4$pep_ratio <- pepNodes$pep_ratio[match(edgelist3$peptide, pepNodes2$peptide)]

  invisible(NULL)

  return(edgelist4)
}


