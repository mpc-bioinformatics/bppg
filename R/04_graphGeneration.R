# Functions in this file:
# .imputationFilter
# .getContractMapping
# .contractGraph
# generateGraphsFromEdgelist
# generateQuantGraphs

#' Filter peptide ratios to exclude in peptide nodes contradicting
#' imputed values.
#' @param G     graph on which the filtering to apply
#' @return  Graph which has filtered some obviouse wrong imputations

.imputationFilter <- function(G) {
    print("currently not used")
}


#' Create Mapping signature for igraph::contract function.
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
#' @return A list with two dataframes, one for peptide and one for protein
#'         signatures.
#'
#' @seealso For collapsing graphs: \cr
#'          [generateGraphsFromFASTA()], [generateQuantGraphs()]
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' edgelist <- digestFASTA(fasta)
#' res <- bppg:::.getContractMapping(edgelist)
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
#' @param G                        \strong{igraph graph} \cr
#'                                 Bipartrite peptide protein graph.
#' @param vMapping                 \strong{list} \cr
#'                                 A list with two dataframes from
#'                                 [.getContractMapping()].
#' @param collProtNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the protein nodes
#'                                 will be collapsed.
#' @param collPeptNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the peptide nodes
#'                                 will be collapsed.
#'
#' @return An edgelist with collapsed protein and/or peptide nodes.
#'
#' @seealso For graph collapsing: [.getContractMapping()] \cr
#'          [generateGraphsFromEdgelist()]
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' edgelist <- digestFASTA(fasta)
#' vMapping <- bppg:::.getContractMapping(edgelist)
#' G <- igraph::graph_from_edgelist(as.matrix(edgelist), directed = FALSE)
#' igraph::V(G)[igraph::V(G)$name %in% edgelist[, 1]]$type <- TRUE
#' igraph::V(G)[igraph::V(G)$name %in% edgelist[, 2]]$type <- FALSE
#' res <- bppg:::.contractGraph(G, vMapping)

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

    gCollapsed <- igraph::contract(G, 
        factor(stats::na.omit(igraph::V(G)$collSignature)),
        vertex.attr.comb = c)

    # remove duplicate edges
    gCollapsed  <- igraph::simplify(gCollapsed)

    # reset general attributes
    igraph::V(gCollapsed)$type <- vapply(igraph::V(gCollapsed)$type, "[", 1,
        FUN.VALUE = logical(1))
    igraph::V(gCollapsed)$name <- vapply(igraph::V(gCollapsed)$name,
        paste, collapse=";", FUN.VALUE = character(1)) 
    # this is not ordered - > same ratio order


    if (!is.null(igraph::V(gCollapsed)$pep_logRatio)) {
        pepMask <- !igraph::V(gCollapsed)$type
        if (collPeptNodes) {
            igraph::V(gCollapsed)$pep_ratio_mean[pepMask] <-
                vapply(igraph::V(gCollapsed)$pep_logRatio[pepMask],
                    mean, FUN.VALUE = numeric(1))
            if (!is.null(igraph::V(gCollapsed)$imputed)) {
# how to combine this one?
            }
        } else {
            igraph::V(gCollapsed)$pep_logRatio <- vapply(
                igraph::V(gCollapsed)$pep_logRatio,  "[", 1,
                FUN.VALUE = numeric(1))
        } 
    }

    if (!is.null(igraph::V(gCollapsed)$protOrigin) && collProtNodes) {
        igraph::V(gCollapsed)$protOrigin[igraph::V(gCollapsed)$type] <-
            vapply(igraph::V(gCollapsed)$protOrigin[igraph::V(gCollapsed)$type],
                unique, FUN.VALUE = character(1))
    }

    return(igraph::delete_vertex_attr(gCollapsed, "collSignature"))
}


#' Generate bipartite peptide-protein graphs from a list of digested proteins
#' via an edgelist. Peptide and protein nodes can be contracted.
#'
#' @param edgelist                 \strong{data.frame} \cr
#'                                 An edgelist, output from [digestFASTA()].
#'                                 Quant data needs to be in the column
#'                                 \strong{$pep_logRatio}.
#' @param collProtNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the protein nodes
#'                                 will be contracted.
#' @param collPeptNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the peptide nodes
#'                                 will be contracted.
#' @return A list of subgraphs as igraph objects.
#' @export
#'
#' @seealso [digestFASTA()]
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' edgelist <- digestFASTA(fasta)
#' res <- bppg::generateGraphsFromEdgelist(edgelist)
#'
generateGraphsFromEdgelist <- function(edgelist,
    collProtNodes = FALSE,
    collPeptNodes = FALSE) {
    checkmate::assertDataFrame(edgelist)
    checkmate::assertFlag(collProtNodes)
    checkmate::assertFlag(collPeptNodes)

    if (collProtNodes || collPeptNodes) {
        vertexMapping <- .getContractMapping(edgelist, collProtNodes,
            collPeptNodes)
    }

    #generate graph from edge matrix
    G <- igraph::graph_from_edgelist(as.matrix(edgelist[, 1:2]),
        directed = FALSE)

    #assign vertex types to proteins and peptides for the graph to be bipartite
    igraph::V(G)[igraph::V(G)$name %in% edgelist[, 1]]$type <- TRUE
    igraph::V(G)[igraph::V(G)$name %in% edgelist[, 2]]$type <- FALSE

    if (!is.null(edgelist$pep_logRatio)) {
        matchIndex <- match(igraph::V(G)$name[!igraph::V(G)$type],
            edgelist$peptide)
        G <- igraph::set_vertex_attr(graph = G,
            name = "pep_logRatio",
            index = igraph::V(G)[!igraph::V(G)$type],
            value = edgelist$pep_logRatio[matchIndex])
        if (!is.null(edgelist$imputed)) {
            G <- igraph::set_vertex_attr(graph = G,
                name = "imputed",
                index = igraph::V(G)[!igraph::V(G)$type],
                value = edgelist$imputed[matchIndex])
            impProt <- stats::aggregate(imputed ~ protein, edgelist, any)
            G <- igraph::set_vertex_attr(graph = G,
                name = "imputed",
                index = igraph::V(G)[igraph::V(G)$type],
                value = impProt$imputed[
                    match(igraph::V(G)$name[igraph::V(G)$type],
                        impProt$protein)])
        } 
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
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' edgelist <- digestFASTA(fasta)
#'
#' file <- system.file("extdata", "peptides.txt", package = "bppg")
#' D <- readMqPeptideTable(path = file, LFQ = TRUE, remove_contaminants = FALSE)
#' group <- factor(rep(1:9, each = 3))
#' dAgg <- aggregateReplicates(D, group = group)
#' exp_peptide_ratios <- calculatePeptideRatios(dAgg)
#'
#' res <- generateQuantGraphs(exp_peptide_ratios, edgelist)

generateQuantGraphs <- function(exp_peptide_ratios,
    fasta_edgelist,
    seq_column = "Sequence", ## How to assert? could be int
    outpath = NULL,
    collProtNodes = TRUE,
    collPeptNodes = FALSE,
    suffix = "") {
    checkmate::assertClass(exp_peptide_ratios, "SummarizedExperiment")
    checkmate::assertDataFrame(SummarizedExperiment::assays(
        exp_peptide_ratios)$logRatios, all.missing=FALSE)
    checkmate::assert_logical(S4Vectors::metadata(exp_peptide_ratios)$imputed)
    checkmate::assertDataFrame(fasta_edgelist)
    checkmate::assertFlag(collProtNodes)
    checkmate::assertFlag(collPeptNodes)
    checkmate::assertCharacter(suffix)

    ## broad filtering for FASTA edgelist for only quantified peptides
    edgelist_filtered <- fasta_edgelist[fasta_edgelist[, 2]
        %in% SummarizedExperiment::rowData(exp_peptide_ratios)[, seq_column], ]

    if (!is.null(outpath)) {
        checkmate::assertPathForOutput(outpath, overwrite = TRUE)
        openxlsx::write.xlsx(edgelist_filtered,
            file = file.path(outpath,
                paste0("edgelist_filtered_", suffix, ".xlsx")),
            overwrite = TRUE, keepNA = TRUE)
    }
    colnames_split <- limma::strsplit2(colnames(exp_peptide_ratios), "_")
    comparisons <- paste(colnames_split[, 2], colnames_split[, 3], sep = "_")
    # first built graphs and try to identify missing type after
    subgraphs <- lapply(seq_len(ncol(exp_peptide_ratios)),
        function(i) {
            # only the data from this comparison
            compSE <- exp_peptide_ratios[, i, drop = FALSE]
            # remove peptides with missing values
            compSE <- compSE[!is.na(SummarizedExperiment::assays(
                compSE)$logRatios), drop = FALSE]
            # get data
            compRatio <- SummarizedExperiment::assays(compSE)$logRatios

            compEdgelist <- edgelist_filtered[edgelist_filtered$peptide
                %in% rownames(compRatio), ]
            matchIndex <- match(compEdgelist$peptide, 
                rownames(compRatio))
            compEdgelist$pep_logRatio <- compRatio[matchIndex, 1]

            if (S4Vectors::metadata(compSE)$imputed) {
                compEdgelist$imputed <- SummarizedExperiment::assays(
                    compSE)$maskImputed[matchIndex, 1]
            } 
            # TODO does generate Graph has to be adapted too?
            generateGraphsFromEdgelist(compEdgelist,
                collProtNodes, collPeptNodes)
        })

    names(subgraphs) <- comparisons
    return(subgraphs)
}
