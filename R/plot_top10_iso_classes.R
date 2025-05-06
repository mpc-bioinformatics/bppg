

#' Plots the 10 most common isomorphism classes from a prototypelist
#'
#' @param prototypelist prototype list, result of generatePrototypeList()
#' @param out_file output file
#' @param plot_height heigh of plot (in inches)
#' @param plot_width width of plot (in inches)
#'
#' @return nothing, save plot as pdf
#' @export
#'
#' @examples # TODO
plot_top10_iso_classes <- function(prototypelist, out_file, plot_height = 10, plot_width = 4) {

  ### sort prototypelist decreasing by counter
  ord <- order(prototypelist$counter, decreasing = TRUE)
  graphs_tmp <- prototypelist$graphs[ord]
  counter_tmp <- prototypelist$counter[ord]

  grDevices::pdf(out_file, height = plot_height, width = plot_width)#, res = 300, units = "cm")
  graphics::layout(matrix(c(1:11,11), 6, 2, byrow = TRUE), widths=c(1,1), heights=c(1,1,1,1,1,1))
  graphics::par(mai = c(0.1, 0.5, 0.4, 0), oma = c(0,0,0,0), xpd=NA)
  for (i in 1:10) {
    title = paste0(i, ") ", counter_tmp[i], " (", formatC(counter_tmp[i]/sum(counter_tmp)*100, digits = 2, format = "f"), "%)")


    plotBipartiteGraph(graphs_tmp[[i]], vertex.label.dist = 0, legend = FALSE,
                       vertex.color = c("mediumseagreen", "cadetblue2", "coral1"),
                       vertex.size = 35, vertex.label.cex = 1.5, edge.width = 3, vertex.size2=35,
                       useCanonicalPermutation = TRUE, three_shapes = TRUE,
                       node_labels_proteins = "letters",
                       node_labels_peptides = "numbers",
                       round_digits = 2)
    title(main = title, cex.main = 1.2, line = 0.5)
  }
  vertex.color = c("mediumseagreen", "cadetblue2", "coral1")
  plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=c(0, 0.4), ylim=c(0, 0.1),
       oma = c(0,0,0,0), mar = c(0,0,0,0), xpd=NA)
  graphics::legend(0.2, 0.12, legend = c("protein", "shared peptide", "unique peptide"),
         col = vertex.color, pch = c(19, 15, 18), cex = 2, pt.cex = c(2.2, 2.2, 2.4), bty = "n",
         xjust = 0.5, yjust = 1)
  grDevices::dev.off()

}
