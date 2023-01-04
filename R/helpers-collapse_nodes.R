
### TODO: umschreiben, wenn Studienprojekt WS22/23 fertig ist, für effizientere Lösung,
###       die auch die Namen der Proteine behält.
###       Diese Version müsste dann auch Peptid-Knoten collapsen können, sodass man nur eine
###       Funktion braucht.


#' collapse protein nodes (merges duplicated columns of biadjacency submatrices)
#'
#' @param subgraphs subgraphs list
#' @param sparse TRUE, if sparse matrix type is used
#' @param fast TRUE, if combining of column names should be skipped (makes calculation faster if there are large subgraphs)
#' @param fc Are peptide ratios included in the graphs?
#'
#' @return submatrizes or subgraphs, depending on matrix argument
#' @export
#'
#' @examples
#' ### TODO
#'
collapse_protein_nodes <- function(subgraphs, sparse = FALSE, fast = FALSE, fc = FALSE) {

  submatrix <- lapply(subgraphs, igraph::as_incidence_matrix)

  if(sparse) submatrix <- lapply(submatrix, function(x) as(x, "dgCMatrix"))

  if (fc) {
  peptide_ratios <- lapply(subgraphs, function(x) igraph::V(x)$pep_ratio)
  peptide_ratios <- lapply(peptide_ratios, stats::na.omit)
  }

  for (i in seq_along(submatrix)) { # for each submatrix in submatrix
   # print(i)
    tmp <- submatrix[[i]]




    if(ncol(tmp) == 1) next

    if (sparse) {
      ind <- duplicated(tmp, MARGIN = 2) ## duplicated.dgCMatrix
      if (!any(ind)) next
      tmp2 <- tmp[, !ind, drop = FALSE]
      tmp3 <- tmp[, ind, drop = FALSE]

      if (!fast) {
        for (j in 1:ncol(tmp2)) {
          for (k in 1:ncol(tmp3)) {

            if (all(tmp3[,k] == tmp2[,j]))  {
              groupname <- BBmisc::collapse(c(colnames(tmp3)[k], colnames(tmp2)[j]), sep = ";")
              colnames(tmp2)[j] <- groupname
            }
          }
        }
      }

    } else {
      ind <- duplicated(tmp, MARGIN = 2)
      tmp2 <- tmp[, !ind, drop = FALSE]

      if (!fast) {
        for (j in 1:ncol(tmp2)) {
          ind <- apply(tmp, 2, function(x) all(x == tmp2[,j]))
          groupname <- BBmisc::collapse(colnames(tmp)[ind], sep = ";")
          colnames(tmp2)[j] <- groupname
        }
      }

    }
    submatrix[[i]] <- tmp2
  }

  subgraphs <- lapply(submatrix, igraph::graph_from_incidence_matrix)

  if (fc) {
    for (i in 1:length(subgraphs)) {
      G_tmp <- subgraphs[[i]]
      peptide_ratios_tmp <- peptide_ratios[[i]]

      igraph::vertex_attr(G_tmp, "pep_ratio", index = igraph::V(G_tmp)[!igraph::V(G_tmp)$type]) <- peptide_ratios_tmp

      subgraphs[[i]] <- G_tmp
    }
  }



  return(subgraphs)
}




##### TODO: hier weitermachen!!
### FC Berechnung und dann als Attribut direkt in den Graphen!


#' Function to callapse peptide nodes (also the peptide ratios if present)
#'
#' @param subgraphs Submatrix list
#' @param sparse TRUE, if sparse matrix type is used
#' @param fc Are peptide ratios included?
#' @param fast TRUE, if combining of column names should be skipped (makes calculation faster if there are large subgraphs)
#'
#' @return Graphs with collapsed peptide nodes
#' @export
#'
#'
#' @examples
#' ### TODO
collapse_peptide_nodes <- function(subgraphs, sparse = FALSE, fc = TRUE, fast = FALSE) {


  submatrix <- lapply(subgraphs, igraph::as_incidence_matrix)
  if(sparse) submatrix <- lapply(submatrix, function(x) as(x, "dgCMatrix"))


  if (fc) {
    peptide_ratios <- lapply(subgraphs, function(x) igraph::V(x)$pep_ratio)
    peptide_ratios <- lapply(peptide_ratios, stats::na.omit)

    peptide_ratios2 <- peptide_ratios
  }



  if (fc) {
    for (i in seq_along(submatrix)) {

      tmp <- submatrix[[i]]

      #fc <- submatrix[[i]]$fc
      ind <- duplicated(tmp, MARGIN = 1)
      tmp2 <- tmp[!ind, , drop = FALSE]

      peptide_ratios_tmp <- peptide_ratios[[i]]

      fc_tmp_geom <- rep(NA, nrow(tmp2))
      fc_tmp <- rep(NA, nrow(tmp2))

      for (j in 1:nrow(tmp2)) {
        ind <- apply(tmp, 1, function(x) all(x == tmp2[j,]))
        #print(peptide_ratios_tmp[ind])
        fc_tmp_geom[j] <- 2^mean(log2(peptide_ratios_tmp[ind]))
        fc_tmp[j] <- BBmisc::collapse(peptide_ratios_tmp[ind], sep = ";")

        groupname <- BBmisc::collapse(rownames(tmp)[ind], sep = ";")
        rownames(tmp2)[j] <- groupname
      }
      submatrix[[i]] <- tmp2

      peptide_ratios[[i]] <- fc_tmp_geom
      peptide_ratios2[[i]] <- fc_tmp

    }
  } else {
    for (i in seq_along(submatrix)) {
      #print(i)
      tmp <- submatrix[[i]]

      if(sparse) {
        ind <- duplicated(tmp, MARGIN = 1) #duplicated.dgCMatrix

        if (!any(ind)) next

        tmp2 <- tmp[!ind, , drop = FALSE]
        tmp3 <- tmp[ind, , drop = FALSE]

        if (!fast) {
          for (j in 1:nrow(tmp2)) {
            for (k in 1:nrow(tmp3)) {

              if (all(tmp3[k,] == tmp2[j,]))  {
                groupname <- BBmisc::collapse(c(rownames(tmp3)[k], rownames(tmp2)[j]), sep = ";")
                rownames(tmp2)[j] <- groupname
              }
            }
          }
        }

        submatrix[[i]]<- tmp2

      } else {
        ind <- duplicated(tmp, MARGIN = 1)
        tmp2 <- tmp[!ind, , drop = FALSE]

        fc_tmp <- rep(NA, nrow(tmp2))
        for (j in 1:nrow(tmp2)) {
          ind <- apply(tmp, 1, function(x) all(x == tmp2[j,]))
          groupname <- BBmisc::collapse(rownames(tmp)[ind], sep = ";")
          rownames(tmp2)[j] <- groupname
        }
        submatrix[[i]] <- tmp2
      }
    }
  }

  subgraphs <- lapply(submatrix, igraph::graph_from_incidence_matrix)

  if (fc) {
    for (i in 1:length(subgraphs)) {
      G_tmp <- subgraphs[[i]]
      peptide_ratios_tmp <- peptide_ratios[[i]]
      peptide_ratios2_tmp <- peptide_ratios2[[i]]

      igraph::vertex_attr(G_tmp, "pep_ratio", index = igraph::V(G_tmp)[!igraph::V(G_tmp)$type]) <- peptide_ratios_tmp
      igraph::vertex_attr(G_tmp, "pep_ratio2", index = igraph::V(G_tmp)[!igraph::V(G_tmp)$type]) <- peptide_ratios2_tmp

      subgraphs[[i]] <- G_tmp
    }
  }

  return(subgraphs)

}




