
# TODO: LTS Normalisierung einbauen, wie in stat_workflows!


automatedNormalization <- function(D, D.name = deparse(substitute(D)),
                                   method = "loess", suffix = method, log = TRUE, id = NULL,
                                   output_path = "", save = FALSE,
                                   groupwise = FALSE, group = NULL) {

  if(method == "loess" | method == "quantile" | method == "median") {

    if(log) {
      D <- log2(D)
    }

    #### choose normalization function
    fun <- switch(method,
                  "loess" = limma::normalizeBetweenArrays,
                  "quantile" = limma::normalizeBetweenArrays,
                  "median" = limma::normalizeBetweenArrays)

    ### choose arguments for normalization function
    args <- switch(method,
                   "loess" = list(object = D, method = "cyclicloess"),
                   "quantile" = list(object = D, method = "quantile"),
                   "median" = list(object = D, method = "scale"))

    if (!groupwise) {
      D_norm <- do.call(fun, args)
      D_norm <- as.data.frame(D_norm)
    } else {
      D_split <- split.default(D, group)
      D_split_norm <- lapply(D_split, limma::normalizeBetweenArrays, method = args$method)
      D_norm <- do.call(cbind, D_split_norm)
    }

    # if (length(id_columns) >= 1 ) {
    D_norm_2 <- data.frame(id, D_norm)
    # }

    #  tryCatch(expr = {
    if (save) {
      openxlsx::write.xlsx(x = D_norm_2, file = paste0(output_path, D.name, "_", suffix, ".xlsx"), keepNA = TRUE, overwrite = TRUE)
      message("Normalized data successfully saved!")
    }
    # },
    # error = function(err) {
    # error handler picks up where error was generated
    #  print(paste("MY_ERROR:  ",err))
    #  beepr::beep(sound = 10)
    ##  user_input <- readline(prompt = paste0("+++ Do you want to overwrite ", paste0(DATA.name,"_",method,".xlsx"), "? +++ [yes/no] "))
    #  if(user_input == "yes"){
    #    write.xlsx(x = DATA_norm_2, file = paste0(output_path, DATA.name,"_",suffix,".xlsx"), overwrite = TRUE, keepNA = TRUE)
    #    message("Normalized data successfully saved!")
    #  } else {
    #    message("Overwriting of normalized data failed. Please allow overwriting, remove the data file or choose different normalization method!")
    #   }
    #  })
  }else{ # if method == "nonorm"
    D_norm <- D
    cat("No normalization applied.")

    #  if (ncol(id_columns) >= 1 ) {
    D_norm_2 <- data.frame(id, D_norm)
    # }
    openxlsx::write.xlsx(x = D_norm_2, file = paste0(output_path, D.name, "_", suffix, ".xlsx"), keepNA = TRUE)
  }

  return(D_norm)
}









#### Function for a single MA-Plot:
# x1: Sample 1
# x2: Sample 2
# log: Should data be log-transformed?
#      TRUE, if not already log-transformed, FALSE, if already log-transformed
# alpha: Should points be transparent?
# col: colours of the data points
# ...: further arguments for ma.plot
MAPlot_single <- function(x1, x2, log = TRUE, alpha = FALSE, col = "black", ...) {

  if(log) {
    x1 <- log2(x1)
    x2 <- log2(x2)
  }
  if(alpha) {
    col = alpha(col, 0.5)
  }

  M <- na.omit(x1 - x2)
  A <- na.omit((x1 + x2)/2)

  if (length(col) > 1) {
    na.ind <- attr(M, "na.action")
    col <- col[-na.ind]
  }


  affy::ma.plot(A = A, M = M, pch = 16, cex = 0.7, col = col, show.statistics = FALSE, ...)
}



# function to check if user is sure to plot more than 1000 plots --> if yes it return 1 and the plots will be created
MAPlots_check <- function(X, maxPlots, ...){
  number_states <- max(as.integer(as.factor(colnames(X))))
  number_plots <- choose(number_states,2)
  return_value <- 2

  if(number_plots >= maxPlots){
    #beepr::beep(sound = 10)
    user_input <- readline(prompt = paste("Are you sure you want to generate",number_plots,"MA-plots? [yes/no]"))
    if (user_input == "yes"){
      return_value <- 1
    }else return_value <- 0
  }else return_value <- 1

  return(return_value)
}



### main function for MA-Plots
# X: Data in wide format
# labels: labels of the samples for the title of the MA-Plot
# labels2: second line in title, e.g. group membership
MAPlots <- function(X, log = TRUE, alpha = FALSE, suffix="nonorm",
                    labels = 1:ncol(X), labels2 = colnames(X), maxPlots = 5000,
                    plot_height=15, plot_width=15, output_path = "", ...) {

  require(limma)
  require(affy)
  require(scales)
  require(beepr)

  number_states <- max(as.integer(as.factor(colnames(X))))
  number_plots <- choose(number_states,2)

  if(MAPlots_check(X, maxPlots) == 1){

    num <- 0

    print("Generating MA-Plots ...")

    ### TODO: auf pbapply umsteigen
    pb <- utils::txtProgressBar(min = 0,max = number_plots,char = "#",style = 3)

    pdf(paste0(output_path, "MA_Plots_", suffix, ".pdf"), height = plot_height/2.54, width = plot_width/2.54)

    for(i in 1:(ncol(X)-1)) {
      for (j in (i + 1):ncol(X)) {

        if (is.null(labels2)) {
          main = paste(labels[i], labels[j])
        } else  {
          main = paste(labels[i], labels[j], "\n", labels2[i], labels2[j])
        }

        num <- num + 1
        utils::setTxtProgressBar(pb, num)

        MAPlot_single(X[,i], X[, j], log = log, main = main, ...)
      }
    }
    # sound chosen, "treasure", "facebook" also cool :)
    #beepr::beep("coin")
    close(pb)
    print("MA-Plots finished!")

    dev.off()

  }

}


