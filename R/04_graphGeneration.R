# Functions in this file:
# .getContractMapping
# .contractGraph
# generateGraphsFromEdgelist()
# generateQuantGraphs()


#' Create Mapping signature for igraph::contract function.
#'
#' @inheritParams generateQuantGraphs
#' @param edgelist                 \strong{data.frame} \cr
#'                                 An edgelist eg. created with
#'                                 [digestFASTA()].
#'
#' @return A list with two dataframes, one for peptide and one for protein
#'         signatures.
#'
#' @seealso For collapsing graphs: \cr
#'          [generateGraphsFromFASTA()], [generateQuantGraphs()]
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' edgelist <- digestFASTA(fasta)
#' res <- bppg:::.getContractMapping(edgelist)
#'
#' @importFrom stats aggregate

.getContractMapping <- function(edgelist,
    collProtNodes = TRUE,
    collPeptNodes = FALSE) {
    if (collProtNodes) {
        protSignature <- stats::aggregate(data = edgelist,
        x = peptide ~ protein,
        function(x) paste(sort(unique(x)), collapse = ";"))
    } else {
        # map back to themselves, not to peptide set
        # this then mixes with  peptide groups
        protSignature <- edgelist
        protSignature$peptide <- paste0("prot_", protSignature$protein)
    }
    if (collPeptNodes) {
        peptSignature <- stats::aggregate(data = edgelist,
        x = protein ~ peptide,
        function(x) paste(sort(unique(x)), collapse = ";"))
    } else {
        peptSignature <- edgelist
        peptSignature$protein <- paste0("pept_", peptSignature$peptide)
    }
    return(list(peptides = peptSignature, proteins = protSignature))
}

#' Contracting of peptide and protein nodes.
#'
#' @inheritParams generateQuantGraphs
#' @param G                        \strong{igraph graph} \cr
#'                                 Bipartrite peptide protein graph.
#' @param vMapping                 \strong{list} \cr
#'                                 A list with two dataframes from
#'                                 [.getContractMapping()].
#'
#' @return An edgelist with collapsed protein and/or peptide nodes.
#'
#' @seealso For graph collapsing: [.getContractMapping()] \cr
#'          [generateGraphsFromEdgelist()]
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' edgelist <- digestFASTA(fasta)
#' vMapping <- bppg:::.getContractMapping(edgelist)
#' G <- igraph::graph_from_edgelist(as.matrix(edgelist), directed = FALSE)
#' igraph::V(G)[igraph::V(G)$name %in% edgelist[, 1]]$type <- TRUE
#' igraph::V(G)[igraph::V(G)$name %in% edgelist[, 2]]$type <- FALSE
#' res <- bppg:::.contractGraph(G, vMapping)
#'
#' @importFrom igraph contract set_vertex_attr simplify V
#' @importFrom stats na.omit

.contractGraph <- function(G, vMapping,
    collProtNodes,
    collPeptNodes) {
    G <- igraph::set_vertex_attr(graph = G,
        name = "collSignature",
        index = igraph::V(G)[igraph::V(G)$type],
        value = vMapping$proteins$peptide[
            match(igraph::V(G)$name[igraph::V(G)$type],
                vMapping$proteins$protein)])

    G <- igraph::set_vertex_attr(graph = G,
        name = "collSignature",
        index = igraph::V(G)[!igraph::V(G)$type],
        value = vMapping$peptides$protein[
            match(igraph::V(G)$name[!igraph::V(G)$type],
                vMapping$peptides$peptide)])

    gColl <- igraph::contract(G,
        factor(stats::na.omit(igraph::V(G)$collSignature)),
        vertex.attr.comb = c)

    # remove duplicate edges
    gColl  <- igraph::simplify(gColl)

    # reset attributes
    igraph::V(gColl)$type <- vapply(igraph::V(gColl)$type, "[", 1,
        FUN.VALUE = logical(1))
    igraph::V(gColl)$name <- vapply(igraph::V(gColl)$name, paste, collapse=";",
        FUN.VALUE = character(1))
    # this is not ordered - > same ratio order


    if (!is.null(igraph::V(gColl)$pep_logRatio)) {
        if (collPeptNodes) {
            igraph::V(gColl)$pep_ratio_mean[!igraph::V(gColl)$type] <-
                vapply(igraph::V(gColl)$pep_logRatio[!igraph::V(gColl)$type],
                    mean, FUN.VALUE = numeric(1))
        } else {
            igraph::V(gColl)$pep_logRatio <- vapply(
                igraph::V(gColl)$pep_logRatio,  "[", 1, FUN.VALUE = numeric(1))
        } }

    if (!is.null(igraph::V(gColl)$protOrigin)) {
        if (collProtNodes) {
            igraph::V(gColl)$protOrigin[igraph::V(gColl)$type] <- vapply(
                    igraph::V(gColl)$protOrigin[igraph::V(gColl)$type],
                    function(x) { paste(unique(x), collapse = ";")
                    }, FUN.VALUE = character(1))
        }
        # igraph::V(gColl)$protOrigin[!igraph::V(gColl)$type] <- NA
        igraph::V(gColl)$protOrigin <- vapply(igraph::V(gColl)$protOrigin, "[",
        1, FUN.VALUE = character(1))
    }

    return(igraph::delete_vertex_attr(gColl, "collSignature"))
}


#' Generate bipartite peptide-protein graphs from a list of digested proteins
#' via an edgelist. Peptide and protein nodes can be contracted.
#' @inheritParams generateQuantGraphs
#' @param edgelist                 \strong{data.frame} \cr
#'                                 An edgelist, output from [digestFASTA()].
#'                                 Quant data needs to be in the column
#'                                 \strong{$pep_logRatio}.
#' @return A list of subgraphs as igraph objects.
#' @export
#'
#' @seealso [digestFASTA()]
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' edgelist <- digestFASTA(fasta)
#' res <- bppg::generateGraphsFromEdgelist(edgelist, collProtNodes = TRUE,
#'     collPeptNodes = TRUE)
#'
#' @importFrom igraph graph_from_edgelist set_vertex_attr V
#'
generateGraphsFromEdgelist <- function(edgelist,
    collProtNodes = FALSE,
    collPeptNodes = FALSE) {
    checkmate::assertDataFrame(edgelist)
    checkmate::assertFlag(collProtNodes)
    checkmate::assertFlag(collPeptNodes)

    if (collProtNodes || collPeptNodes) {
        vertexMapping <- .getContractMapping(
            edgelist[, c("protein", "peptide")], collProtNodes, collPeptNodes)
    }

    #generate graph from edge matrix
    G <- igraph::graph_from_edgelist(as.matrix(edgelist[, c(1,2)]),
        directed = FALSE)

    #assign vertex types to proteins and peptides for the graph to be bipartite
    igraph::V(G)[igraph::V(G)$name %in% edgelist[, 1]]$type <- TRUE
    igraph::V(G)[igraph::V(G)$name %in% edgelist[, 2]]$type <- FALSE

    if (!is.null(edgelist$pep_logRatio)) {
        G <- igraph::set_vertex_attr(graph = G,
            name = "pep_logRatio",
            index = igraph::V(G)[!igraph::V(G)$type],
            value = edgelist$pep_logRatio[
                match(igraph::V(G)$name[!igraph::V(G)$type],
                    edgelist$peptide)])
    }

    if (!is.null(edgelist$protOrigin)) {
        protOriginDF <- edgelist[, c("protein", "protOrigin")]
        protOriginDF <- protOriginDF[!duplicated(protOriginDF), ]
        G <- igraph::set_vertex_attr(graph = G,
            name = "protOrigin",
            index = igraph::V(G)[igraph::V(G)$type],
            value = protOriginDF$protOrigin[
                match(igraph::V(G)$name[igraph::V(G)$type],
                    protOriginDF$protein)])
    }

    if (collProtNodes || collPeptNodes) {
        G <- .contractGraph(G, vertexMapping, collProtNodes, collPeptNodes)
    }
    return(igraph::decompose(G))
}



#' Generate graphs from peptide ratio table, using an edgelist calculated
#' on the fasta file.
#'
#' @param exp_peptide_ratios       \strong{SummarizedExperiment} \cr
#'                                 A SummarizedExperiment from
#'                                 [bppg::calculatePeptideRatios] with peptide
#'                                 ratios.
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
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_proteome_Scerevisiae_filtered.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' edgelist <- digestFASTA(fasta)
#'
#' file <- system.file("extdata", "peptides_filtered.txt", package = "bppg")
#' group <- factor(rep(1:9, each = 3))
#' D <- readMqPeptideTable(path = file, group = group, LFQ = TRUE, remove_contaminants = FALSE)
#' D_norm <- bppg::normalizePeptideIntensities(D)
#' dAgg <- aggregateReplicates(D_norm)
#' exp_peptide_ratios <- calculatePeptideRatios(dAgg)
#'
#' res <- generateQuantGraphs(exp_peptide_ratios, edgelist)

generateQuantGraphs <- function(exp_peptide_ratios,
    fasta_edgelist,
    seq_column = "Sequence",
    outpath = NULL,
    collProtNodes = TRUE,
    collPeptNodes = FALSE,
    suffix = "") {
    checkmate::assertClass(exp_peptide_ratios, "SummarizedExperiment")
    checkmate::assertDataFrame(SummarizedExperiment::assays(
        exp_peptide_ratios)$logRatios, all.missing=FALSE)
    checkmate::assertDataFrame(fasta_edgelist)
    checkmate::assertFlag(collProtNodes)
    checkmate::assertFlag(collPeptNodes)
    checkmate::assertCharacter(suffix)

    ## broad filtering for edgelist for only quantifies peptides
    edgelist_filtered <- fasta_edgelist[fasta_edgelist[, 2]
        %in% SummarizedExperiment::rowData(exp_peptide_ratios)[, seq_column], ]

    if (!is.null(outpath)) {
        checkmate::assertPathForOutput(outpath, overwrite = TRUE)
        openxlsx::write.xlsx(edgelist_filtered,
            file = file.path(outpath, paste0("edgelist_filtered_", suffix,
                ".xlsx")), overwrite = TRUE, keepNA = TRUE)
    }

    colnames_split <- limma::strsplit2(colnames(exp_peptide_ratios), "_")
    comparisons <- paste(colnames_split[, 2], colnames_split[, 3], sep = "_")

    subgraphs <- lapply(seq_len(ncol(exp_peptide_ratios)),
        function(i) {
            compRatio <- SummarizedExperiment::assays(
                exp_peptide_ratios)$logRatios[, i, drop=FALSE]
            compRatio <- compRatio[!is.na(compRatio), 1, drop=FALSE]

            compEdgelist <- edgelist_filtered[edgelist_filtered$peptide
                %in% rownames(compRatio), ]
            compEdgelist$pep_logRatio <- compRatio[match(compEdgelist$peptide,
                rownames(compRatio)), 1]
            generateGraphsFromEdgelist(compEdgelist, collProtNodes,
                collPeptNodes)
        })

    names(subgraphs) <- comparisons
    return(subgraphs)
}
