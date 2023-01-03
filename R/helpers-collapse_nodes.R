
### TODO: umschreiben, wenn Studienprojekt WS22/23 fertig ist, für effizientere Lösung,
###       die auch die Namen der Proteine behält.
###       Diese Version müsste dann auch Peptid-Knoten collapsen können, sodass man nur eine
###       Funktion braucht.


#' collapse protein nodes (merges duplicated columns of biadjacency submatrices)
#'
#' @param submatrix Submatrix list
#' @param sparse TRUE, if sparse matrix type is used
#' @param fast TRUE, if combining of column names should be skipped (makes calculation faster if there are large subgraphs)
#' @param matrix FALSE if submatrix is a list of subgraphs
#'
#' @return submatrizes or subgraphs, depending on matrix argument
#' @export
#'
#' @examples
#' ### TODO
#'
collapse_protein_nodes <- function(submatrix, sparse = FALSE, fast = FALSE, matrix = TRUE) {

  if (!matrix) {
    submatrix <- lapply(submatrix, igraph::as_incidence_matrix)
  }


  for (i in seq_along(submatrix)) { # for each submatrix in submatrix
    print(i)
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

  if (!matrix) {
    submatrix <- lapply(submatrix, igraph::graph_from_incidence_matrix)
  }

  return(submatrix)
}




##### TODO: hier weitermachen!!
### FC Berechnung und dann als Attribut direkt in den Graphen!


#' Function to callapse peptide nodes (also the peptide ratios if present)
#'
#' @param submatrix Submatrix list
#' @param sparse TRUE, if sparse matrix type is used
#' @param fc Are peptide ratios included?
#' @param fast TRUE, if combining of column names should be skipped (makes calculation faster if there are large subgraphs)
#' @param matrix FALSE if submatrix is a list of subgraphs
#'
#' @return Graphs with collapsed peptide nodes
#' @export
#'
#'
#' @examples
#' ### TODO
collapse_peptide_nodes <- function(submatrix, sparse = FALSE, fc = TRUE, fast = FALSE, matrix = TRUE) {


  if (!matrix) {
    submatrix <- lapply(submatrix, igraph::as_incidence_matrix)
  }


  if (fc) {
    for (i in seq_along(submatrix)) {

      tmp <- submatrix[[i]]$X
      fc <- submatrix[[i]]$fc
      ind <- duplicated(tmp, MARGIN = 1)
      tmp2 <- tmp[!ind, , drop = FALSE]

      fc_tmp <- rep(NA, nrow(tmp2))
      for (j in 1:nrow(tmp2)) {
        ind <- apply(tmp, 1, function(x) all(x == tmp2[j,]))
        fc_tmp[j] <- 2^mean(log2(fc[ind]))
        groupname <- BBmisc::collapse(rownames(tmp)[ind], sep = ";")
        rownames(tmp2)[j] <- groupname
      }
      submatrix[[i]]$X <- tmp2
      submatrix[[i]]$fc <- fc_tmp
    }
  } else {
    for (i in seq_along(submatrix)) {
      print(i)
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

  if (!matrix) {
    submatrix <- lapply(submatrix, igraph::graph_from_incidence_matrix)
  }

  return(submatrix)

}




