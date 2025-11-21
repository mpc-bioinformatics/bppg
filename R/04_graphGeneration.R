#' Functions in this file:
#' .getContractMapping
#' .contractGraph
#' generateGraphsFromEdgelist()
#' generateQuantGraphs()


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
#' @return A list with two dataframes, one for peptides one for proteins.
#'
#' @seealso For collapsing graphs: \cr
#'          [generateGraphsFromFASTA()], [generateQuantGraphs()]
#'
#' @examples
#' TODO
.getContractMapping <- function(edgelist,
                                  collProtNodes = TRUE,
                                  collPeptNodes = FALSE) {
    if (collProtNodes) {
        protSignature <- stats::aggregate(data = edgelist, x = peptide ~ protein,
        function(x) paste(sort(unique(x)), collapse = ";"))
    } else {
        # map back to themselves, not to peptide set
        # this then mixes with orginal 
        protSignature <- edgelist
        protSignature$peptide <- paste0("prot_", protSignature$protein)
    }
    if (collPeptNodes) {
        peptSignature <- stats::aggregate(data = edgelist, x = protein ~ peptide,
        function(x) paste(sort(unique(x)), collapse = ";"))
    } else {
        peptSignature <- edgelist
        peptSignature$protein <- paste0("pept_", peptSignature$peptide)
    }   
    list(peptides = peptSignature, proteins = protSignature)
}

#' Contracting of peptide and protein nodes of in a graph.
#'
#' @param G                        \strong{igraph graph} \cr
#'                                 bipartrite peptide protein graph
#' @param vMapping                 \strong{list} \cr
#'                                 A list with two dataframes from
#'                                 [.getContractMapping()]
#'
#' @return An edgelist with collapsed protein and/or peptide nodes.
#'
#' @seealso For graph collapsing: [.getContractMapping()] \cr
#'          [generateGraphsFromEdgelist()]
#'
#' @examples
#' TODO
#'
.contractGraph <- function(G, vMapping){
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

    message("contracting")
    gCollapsed <- igraph::contract(G, factor(na.omit(V(G)$collSignature)), 
        vertex.attr.comb = c)

    # remove duplicate edges 
    gCollapsed  <- igraph::simplify(gCollapsed)

    # reset attributes
    # message("resetting attributes")
    igraph::V(gCollapsed)$collSignature <- sapply(
        igraph::V(gCollapsed)$collSignature, "[", 1)
    igraph::V(gCollapsed)$type <- sapply(igraph::V(gCollapsed)$type, "[", 1)
    igraph::V(gCollapsed)$name <- sapply(igraph::V(gCollapsed)$name, 
        paste, collapse=";") # this is not ordered - > aber gleiche order wie ratio

    if (!is.null(igraph::V(gCollapsed)$pep_ratio)) {
        igraph::V(gCollapsed)$pep_ratio_mean[!igraph::V(gCollapsed)$type] <- # TODO geht das überhaupt????
            sapply(igraph::V(gCollapsed)$pep_ratio[!igraph::V(gCollapsed)$type], mean)
    }

    if (!is.null(igraph::V(gCollapsed)$protOrigin)) {
        igraph::V(gCollapsed)$protOrigin[igraph::V(gCollapsed)$type] <- 
            sapply(igraph::V(gCollapsed)$protOrigin[igraph::V(gCollapsed)$type], "[", 1)
    }

    # TODO: remove signature flag
    igraph::delete_vertex_attr(gCollapsed, "collSignature")
}


#' Generate bipartite peptide-protein graphs from a list of digested proteins
#' via an edgelist. Peptide and Proteins can be contracted.
#'
#' @param edgelist                 \strong{data.frame} \cr
#'                                 An edgelist, output from [digestFASTA()].
#'                                 For quant data it needs to be in the column
#'                                 $pep_ratio.
#' @param collProtNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the protein nodes
#'                                 will be collapsed.
#' @param collPeptNodes            \strong{logical} \cr
#'                                 If \code{TRUE}, the peptide nodes
#'                                 will be collapsed.
#' @return A list of subgraphs as igraph objects.
#' @export
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
#' res <- bppg::generateGraphsFromEdgelist(edgelist)
#'
generateGraphsFromEdgelist <- function(edgelist,
                                  collProtNodes = FALSE,
                                  collPeptNodes = FALSE) {
    checkmate::checkDataFrame(edgelist)
    checkmate::checkFlag(collProtNodes)
    checkmate::checkFlag(collPeptNodes)

    if(collProtNodes || collPeptNodes) {                                
        vertexMapping <- .getContractMapping(edgelist, collProtNodes,
        collPeptNodes)  
    } 

    #generate graph from edge matrix
    G <- igraph::graph_from_edgelist(as.matrix(edgelist[, 1:2]),
        directed = FALSE)

    #assign vertex types to proteins and peptides for the graph to be bipartite
    igraph::V(G)[igraph::V(G)$name %in% edgelist[, 1]]$type <- TRUE
    igraph::V(G)[igraph::V(G)$name %in% edgelist[, 2]]$type <- FALSE
    
    if (!is.null(edgelist$pep_ratio)){
        G <- igraph::set_vertex_attr(graph = G,
            name = "pep_ratio",
            index = igraph::V(G)[!igraph::V(G)$type],
            value = edgelist$pep_ratio[
                match(igraph::V(G)$name[!igraph::V(G)$type],
                    edgelist$peptide)])
    }

    if(collProtNodes || collPeptNodes) {                                
        G <- .contractGraph(G, vertexMapping)
    } 
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
#' TODO
generateQuantGraphs <- function(peptide_ratios,
                                  id_cols = 1,
                                  fasta_edgelist,
                                  outpath = NULL,
                                  seq_column = "Sequence",
                                  collProtNodes = TRUE,
                                  collPeptNodes = FALSE,
                                  suffix = "") {
    checkmate::checkDataFrame(peptide_ratios)
    checkmate::checkInteger(id_cols)
    checkmate::checkDataFrame(fasta_edgelist)
    
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

        G <- generateGraphsFromEdgelist(edgelist_filtered2, 
            collProtNodes, collPeptNodes)
        ## set peptide ratios as vertex attributes
        subgraphs[[i]] <- G
        names(subgraphs)[[i]] <- comparison
    }
    return(subgraphs)
}
