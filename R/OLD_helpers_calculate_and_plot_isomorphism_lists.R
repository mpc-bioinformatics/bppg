
### TODO: neue Version davon, die Repräsentanten ausrechnet und diese dann zuordnet?
### TODO: use isomorphic function that considers node tpes

#' Calculation of isomorphism lists from submatrizes or subgraphs
#'
#' @param G list of subgraphs
#' #'
#' @return isomorph_list is a list of indizes that belong in the same isomorphism class
#'         Graphs are graph representatives
#' @export
#'
#' @examples
#' ### TODO
calculateIsomorphList <- function(G) {


  isomorph_list <- list()
  k <- 1

  ### TODO: progress bar

  for (i in 1:length(G)) {
    print(i)
    G_tmp <- G[[i]]
    if (k == 1) {isomorph_list[[k]] <- i; k <- k + 1; next}

    for (j in 1:(k-1)) {
      iso <- igraph::isomorphic(G_tmp, G[[isomorph_list[[j]][1]]]) # compare graph with 1st element if each isomorphism class
      if (iso) {
        cG <- igraph::canonical_permutation(G)
        cG <- igraph::permute(G_tmp, cG$labeling)
        cG2 <- igraph::canonical_permutation(Graphs[[isomorph_list[[j]][1]]])
        cG2 <- igraph::permute(Graphs[[isomorph_list[[j]][1]]], cG2$labeling)

        iso2 <- all(igraph::V(cG)$type == igraph::V(cG2)$type)   # test if graphs are really the same (considering the node types)

        if(iso2) {
          isomorph_list[[j]] <- c(isomorph_list[[j]], i); break
        }

      }


      #  same_nr_of_prot_and_pep <- sum(V(G)$type) == sum(V(Graphs[[isomorph_list[[j]][1]]])$type) # compare number of type 1 nodes
      #  if(iso&same_nr_of_prot_and_pep){isomorph_list[[j]] <- c(isomorph_list[[j]], i); break}  # add graph to group of isomorphic graphs
    }

    if (!iso) {isomorph_list[[k]] <- i; k <- k + 1; next}                      # if it is not isomorphic to an existing class, start a new one

  }
  return(list(isomorph_list = isomorph_list, Graphs = Graphs))
}






#### function that plots list of isomorph classes, sorted by number of occurences
## isomorph_list: result of calculateIsomorphList
## Graphs: list of graphs
## path: path to save plots
## title: if TRUE, title is added with number of occurence and percentage value
## pdf: if TRUE, plot is saved in a single odf, if FALSE as multiple pngs
## cex.title: size of title
## which_graphs: ranks of classes that should be plottet (e.g. 1:10 for top 10 classes)
## mfrow: mfrow agrument of par function to arrance graphs in one plot
## save: if TRUE, graphs will be saved as pdf or png
## ...: further arguments to plotBipartiteGraph

#' function that plots list of isomorph classes, sorted by number of occurences
#'
#' @param isomorph_list result of calculateIsomorphList
#' @param Graphs list of graphs
#' @param path path to save plots
#' @param title if TRUE, title is added with number of occurence and percentage value
#' @param pdf if TRUE, plot is saved in a single odf, if FALSE as multiple pngs
#' @param cex.title size of title
#' @param which_graphs ranks of classes that should be plottet (e.g. 1:10 for top 10 classes)
#' @param mfrow mfrow agrument of par function to arrance graphs in one plot
#' @param save if TRUE, graphs will be saved as pdf or png
#' @param title_format "times+percent" or "percent"
#' @param ... further arguments to plotBipartiteGraph
#' @param height height of plot
#' @param width width of plot
#'
#' @return plots saved as a pdf or multiple png files
#' @export
#'
#' @examples
#' ### TODO
plotIsomorphList <- function(isomorph_list, Graphs, path, title = TRUE, pdf = TRUE,
                             cex.title = 1, which_graphs = NULL, mfrow = c(1,1), save = TRUE,
                             title_format = "times+percent", ..., height = 15, width = 15) {

  ord_le_iso <- order(lengths(isomorph_list), decreasing = TRUE) # order by number of occurrences
  if (!is.null(which_graphs)) ord_le_iso <- ord_le_iso[which_graphs]
  le_iso <- lengths(isomorph_list)  # sizes of isomorph classes
  le_iso_total <- sum(le_iso)       # total number of graphs
  percentages <- round(le_iso/le_iso_total * 100, 2) # percentages (proportion of all graphs)



  if(pdf & save) grDevices::pdf(paste0(path, ".pdf"))
  graphics::par(mfrow = mfrow)
  graphics::par(mai = c(0.1, 0.5, 0.4, 0), oma = c(0,0,0,0), mfrow = mfrow)
  j <- 1
  for (i in ord_le_iso) {
    if (!pdf & save) grDevices::png(paste0(path, "_", j, ".png"), res = 500, units = "cm", height = height, width = width)
    ind <- isomorph_list[[i]][1]  # plot first element for each isomorph group
    G <- Graphs[[ind]]
    types <- igraph::vertex_attr(G)$type # type = 0 peptides, type = 1 proteins
    G <- igraph::set_vertex_attr(G, name = "name", value = c(1:sum(!types), LETTERS[1:sum(types)]))

    bppg::plotBipartiteGraph(G, vertex.label.dist = 0, ...)
    if(title & title_format == "times+percent") title(paste0(le_iso[i], " times (",
                                                             formatC(percentages[i], digits = 2, format = "f"), "%)"), cex.main = cex.title, line = 0.5)
    if(title & title_format == "percent") title(paste0(formatC(percentages[i], digits = 2, format = "f"), "%"), cex.main = cex.title, line = 0.5)
    if(!pdf & save) grDevices::dev.off()
    j <- j + 1
  }
  graphics::par(mfrow = c(1,1))
  if(pdf & save) grDevices::dev.off()


}

