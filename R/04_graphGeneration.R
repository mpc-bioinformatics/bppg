#' Functions in this file:
#' .collapseEdgelist()
#' .collapseEdgelistQuant()
#' .generateGraphsFromEdgelist()
#' .generateQuantGraphs()



#' Collapsing of peptide and protein nodes of an edgelist.
#'
#' @param edgelist                 \strong{data.frame} \cr
#'                                 An edgelist eg. created with
#'                                 [digestFASTA()].
#' @param collProtNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the protein nodes
#'                                 will be collapsed.
#' @param collPeptNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the peptide nodes
#'                                 will be collapsed.
#'
#' @return An edgelist with collapsed protein and/or peptide nodes.
#'
#'
#' @seealso For edgelists with peptide ratios: [.collapseEdgelistQuant()] \cr
#'          [generateGraphsFromFASTA()], [.generateQuantGraphs()],
#'          [digestFASTA()]
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' edgelist <- bppg::digestFASTA(fasta)
#' edgelist_collapsed <- bppg:::.collapseEdgelist(edgelist)
#'

.collapseEdgelist <- function(edgelist,
    collProtNodes = TRUE,
    collPeptNodes = TRUE) {

    if (!collProtNodes && !collPeptNodes) {
        return(edgelist)
    }

    ## Calculate list if protein nodes
    if (collProtNodes) {
        ## aggregate peptide sequences that belong to the same protein accession
        ## (1 row per protein accession)
        protEdges <- stats::aggregate(data = edgelist, x = peptide ~ protein,
            function(x) paste(sort(unique(x)), collapse = ";"))
        ## aggregate proteins with the same set of peptides (-> protein nodes)
        protNodes <- stats::aggregate(data = protEdges, x = protein ~ peptide,
            function(x) paste(sort(unique(x)), collapse = ";"))
    } else {
        protEdges <- stats::aggregate(data = edgelist, x = peptide ~ protein,
            function(x) paste(sort(unique(x)), collapse = ";"))
        protNodes <- protEdges
    }

    ## calculate list of peptide nodes
    if (collPeptNodes) {
        ## aggregate protein accessions belonging to the same peptide sequences
        ## (1 row per peptide sequence)
        pepEdges <- stats::aggregate(data = edgelist, x = protein ~ peptide,
            function(x) paste(sort(unique(x)), collapse = ";"))
        ## aggregate peptides with the same set of proteins (-> peptide nodes)
        pepNodes <- stats::aggregate(data = pepEdges, x = peptide ~ protein,
            function(x) paste(sort(unique(x)), collapse = ";"))
    } else {
        pepEdges <- stats::aggregate(data = edgelist, x = protein ~ peptide,
            function(x) paste(sort(unique(x)), collapse = ";"))
        pepNodes <- pepEdges
    }

    edgelist2 <- edgelist

    pepNodes2 <- pepNodes
    ## first peptide from list
    pepNodes2$peptide <- limma::strsplit2(pepNodes2$peptide, ";")[, 1]
    edgelist2 <- edgelist[edgelist$peptide %in% pepNodes2$peptide, ]

    protNodes2 <- protNodes
    ## first peptide from list
    protNodes2$protein <- limma::strsplit2(protNodes2$protein, ";")[, 1]
    edgelist3 <- edgelist2[edgelist2$protein %in% protNodes2$protein, ]

    edgelist4 <- edgelist3
    edgelist4$protein <- protNodes$protein[match(edgelist3$protein,
            protNodes2$protein)]
    edgelist4$peptide <- pepNodes$peptide[match(edgelist3$peptide,
            pepNodes2$peptide)]
    invisible(NULL)
    return(edgelist4)
}

#' Collapsing of peptide and protein nodes of an edgelist.
#'
#' @param edgelist                 \strong{data.frame} \cr
#'                                 An edgelist  with peptide ratios eg. created with [digestFASTA()].
#' @param collProtNodes   \strong{logical} \cr
#'                                 If \code{TRUE}, the protein nodes will be collapsed.
#' @param collPeptNodes   \strong{logical} \cr
#'                                 If \code{TRUE}, the peptide nodes will be collapsed.
#'
#' @return An edgelist with collapsed protein and/or peptide nodes.
#'
#'
#' @seealso For edgelists without peptide ratios: [.collapseEdgelist()] \cr
#'          [generateGraphsFromFASTA()], [.generateQuantGraphs()], [digestFASTA()]
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' edgelist <- bppg::digestFASTA(fasta)
#' edgelist_collapsed <- bppg:::.collapseEdgelist(edgelist)
#'

#### TODO TO LONG simplified
.collapseEdgelistQuant <- function(edgelist,
    collProtNodes = TRUE,
    collPeptNodes = TRUE) {

    if (!collProtNodes && !collPeptNodes) {
        return(edgelist)
    }

    ## Calculate list if protein nodes
    if (collProtNodes) {
        ## aggregate peptide sequences that belong to the same protein accession
        ## (1 row per protein accession)
        protEdges <- stats::aggregate(data = edgelist,
            x = cbind(peptide, pep_ratio) ~ protein,
            function(x) paste(sort(unique(x)), collapse = ";"))
        ## aggregate proteins with the same set of peptides (-> protein nodes)
        protNodes <- stats::aggregate(data = protEdges,
            x = protein ~ peptide + pep_ratio,
            function(x) paste(sort(unique(x)), collapse = ";"))
    } else {
        protEdges <- stats::aggregate(data = edgelist,
            x = peptide ~ protein,
            function(x) paste(sort(unique(x)), collapse = ";"))
        protNodes <- protEdges
    }

    ## calculate list of peptide nodes
    if (collPeptNodes) {
        ## aggregate protein accessions belonging to the same peptide sequences
        ## (1 row per peptide sequence)
        pepEdges <- stats::aggregate(data = edgelist,
            x = protein ~ peptide + pep_ratio,
            function(x) paste(sort(unique(x)), collapse = ";"))
        ## aggregate peptides with the same set of proteins (-> peptide nodes)
        pepNodes <- stats::aggregate(data = pepEdges,
            x = cbind(peptide, pep_ratio) ~ protein,
            function(x) paste(sort(unique(x)), collapse = ";"))
    } else {
        pepEdges <- stats::aggregate(data = edgelist,
            x = protein ~ peptide + pep_ratio,
            function(x) paste(sort(unique(x)), collapse = ";"),
            simplify = FALSE)
        pepNodes <- pepEdges
    }

    edgelist2 <- edgelist
    pepNodes2 <- pepNodes
    ## first peptide from list
    pepNodes2$peptide <- limma::strsplit2(pepNodes2$peptide, ";")[, 1]
    edgelist2 <- edgelist[edgelist$peptide %in% pepNodes2$peptide, ]

    protNodes2 <- protNodes
    protNodes2$protein <- limma::strsplit2(protNodes2$protein, ";")[, 1]
    ## first peptide from list
    edgelist3 <- edgelist2[edgelist2$protein %in% protNodes2$protein, ]

    edgelist4 <- edgelist3
    edgelist4$protein <- protNodes$protein[match(edgelist3$protein,
            protNodes2$protein)]
    edgelist4$peptide <- pepNodes$peptide[match(edgelist3$peptide,
            pepNodes2$peptide)]
    edgelist4$pep_ratio <- pepNodes$pep_ratio[match(edgelist3$peptide,
            pepNodes2$peptide)]

    invisible(NULL)
    return(edgelist4)
}

#' Generate bipartite peptide-protein graphs from a list of digested proteins
#' via an edgelist.
#'
#' @param edgelist   \strong{data.frame} \cr
#'                   An edgelist, output from [digestFASTA()].
#'
#' @return A list of subgraphs as igraph objects.
#'
#'
#' @seealso [digestFASTA()]
#'
#' @examples
#' ## TODO: example takes longer than 5s
#' library(seqinr)
#' file <- system.file("extdata", "2020_01_31_proteome_S_cerevisae.fasta",
#'  package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' edgelist <- digestFASTA(fasta)
#' res <- bppg:::.generateGraphsFromEdgelist(edgelist)
#'

.generateGraphsFromEdgelist <- function(edgelist) {
    #generate graph from edge matrix
    G <- igraph::graph_from_edgelist(as.matrix(edgelist[,1:2]),
        directed = FALSE)

    #assign vertex types to proteins and peptides for the graph to be bipartite
    igraph::V(G)[igraph::V(G)$name %in% edgelist[, 1]]$type <- TRUE
    igraph::V(G)[igraph::V(G)$name %in% edgelist[, 2]]$type <- FALSE

    #decompose graph into connected components
    igraph::decompose(G)
}



#' Generate graphs from peptide ratio table, using an edgelist calculated
#' on the fasta file.
#'
#' @param peptide_ratios           \strong{data.frame} \cr
#'                                 A table with peptide ratios.
#' @param id_cols                  \strong{integer vector} \cr
#'                                 The columns with ids, e.g. peptide sequences
#'                                 (everything except the peptide ratios)
#' @param fasta_edgelist           \strong{data.frame} \cr
#'                                 An edgelist created from the corresponding
#'                                 FASTA file, eg. created with
#'                                 [bppg::digestFASTA()].
#' @param outpath                  \strong{character} \cr
#'                                 The output path for the results.
#' @param seq_column               \strong{character} \cr
#'                                 The column name of the peptide sequence.
#' @param collProtNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the protein nodes will
#'                                 be collapsed.
#' @param collPeptNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the peptide nodes will
#'                                 be collapsed.
#' @param suffix                   \strong{character} \cr
#'                                 The suffix for saving results.
#'
#' @return A list of list of subgraphs
#' @export
#'
#' @seealso [bppg::digestFASTA()]
#'
#' @examples
#'

.generateQuantGraphs <- function(peptide_ratios,
                                  id_cols = 1,
                                  fasta_edgelist,
                                  outpath = NULL,
                                  seq_column = "Sequence",
                                  collProtNodes = TRUE,
                                  collPeptNodes = FALSE,
                                  suffix = "") {
    ## broad filtering for edgelist for only quantifies peptides
    edgelist_filtered <- fasta_edgelist[fasta_edgelist[, 2]
        %in% peptide_ratios[, seq_column], ]

    if (!is.null(outpath)) {
        openxlsx::write.xlsx(edgelist_filtered,
            file = file.path(outpath, paste0("edgelist_filtered_", suffix, ".xlsx")),
            overwrite = TRUE, keepNA = TRUE)
    }

    id <- peptide_ratios[, id_cols, drop = FALSE]
    peptide_ratios <- peptide_ratios[, -(id_cols), drop = FALSE]
    colnames_split <- limma::strsplit2(colnames(peptide_ratios), "_")
    comparisons <- paste(colnames_split[,2], colnames_split[,3], sep = "_")

    subgraphs <- list()
    for (i in 1:ncol(peptide_ratios)) {
        comparison <- comparisons[i]
        fc <- peptide_ratios[,i]
        ## peptides that are quantified in this specific comparison
        peptides_tmp <- id[, seq_column][!is.na(fc)]
        fc <- stats::na.omit(fc)
        edgelist_filtered2 <- edgelist_filtered[edgelist_filtered[, 2]
            %in% peptides_tmp, ]

        ## add peptide ratios
        edgelist_filtered2$pep_ratio <- peptide_ratios[, i][
            match(edgelist_filtered2$peptide, id[, seq_column])]

        ## generate whole bipartite graph
        edgelist_coll <- .collapseEdgelistQuant(edgelist_filtered2,
            collProtNodes = collProtNodes, collPeptNodes = collPeptNodes)

        G <- .generateGraphsFromEdgelist(edgelist_coll[, 1:2])
        ## set peptide ratios as vertex attributes
        ## TODO APPlY
        for (j in 1:length(G)){
            G[[j]] <- igraph::set_vertex_attr(graph = G[[j]],
                name = "pep_ratio",
                index = igraph::V(G[[j]])[!igraph::V(G[[j]])$type],
                value = edgelist_coll$pep_ratio[
                    match(igraph::V(G[[j]])$name[!igraph::V(G[[j]])$type],
                        edgelist_coll$peptide)])
        }
        subgraphs[[i]] <- G
        names(subgraphs)[[i]] <- comparison
    }
    return(subgraphs)
}
