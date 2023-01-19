# affy::ma.plot
# function (A, M, subset = sample(1:length(M), min(c(10000, length(M)))),
#           show.statistics = TRUE, span = 2/3, family.loess = "gaussian",
#           cex = 2, plot.method = c("normal", "smoothScatter",
#                                    "add"), add.loess = TRUE, lwd = 1, lty = 1, loess.col = "red",
#           ...)
# {
#   plot.method <- match.arg(plot.method)
#   fn.call <- list(...)
#   sigma <- IQR(M)
#   mean <- median(M)
#   if (!is.element("ylim", names(fn.call))) {
#     yloc <- max(M)
#   }
#   else {
#     yloc <- max(fn.call$ylim)
#   }
#   if (!is.element("xlim", names(fn.call))) {
#     xloc <- max(A)
#   }
#   else {
#     xloc <- max(fn.call$xlim)
#   }
#   if (plot.method == "smoothScatter") {
#     plotmethod <- "smoothScatter"
#   }
#   else if (plot.method == "add") {
#     plotmethod <- "add"
#   }
#   else {
#     plotmethod <- "normal"
#   }
#   aux <- loess(M[subset] ~ A[subset], degree = 1, span = span,
#                family = family.loess)$fitted
#   if (plotmethod == "smoothScatter") {
#     smoothScatter(A, M, ...)
#   }
#   else if (plotmethod == "add") {
#     points(A, M, cex = cex, ...)
#   }
#   else {
#     plot(A, M, cex = cex, ...)
#   }
#   if (add.loess) {
#     o <- order(A[subset])
#     A <- A[subset][o]
#     M <- aux[o]
#     o <- which(!duplicated(A))
#     lines(approx(A[o], M[o]), col = loess.col, lwd = lwd,
#           lty = lty)
#   }
#   abline(0, 0, col = "blue")
#   if (show.statistics) {
#     txt <- format(sigma, digits = 3)
#     txt2 <- format(mean, digits = 3)
#     text(xloc, yloc, paste(paste("Median:", txt2),
#                            paste("IQR:", txt), sep = "\n"), cex = cex,
#          adj = c(1, 1))
#   }
# }


################################################################################
################################################################################

loessFit2 <- function (y, x, weights = NULL, span = 0.3, iterations = 4L,
          min.weight = 1e-05, max.weight = 1e+05, equal.weights.as.null = TRUE,
          method = "weightedLowess") {
  n <- length(y)
  if (length(x) != n)
    stop("y and x have different lengths")
  out <- list(fitted = rep(NA, n), residuals = rep(NA, n))
  obs <- is.finite(y) & is.finite(x)
  xobs <- x[obs]
  yobs <- y[obs]
  nobs <- length(yobs)
  if (nobs == 0)
    return(out)
  if (span < 1/nobs) {
    out$fitted[obs] <- y[obs]
    out$residuals[obs] <- 0
    return(out)
  }
  if (min.weight < 0)
    min.weight <- 0
  if (!is.null(weights)) {
    if (length(weights) != n)
      stop("y and weights have different lengths")
    wobs <- weights[obs]
    wobs[is.na(wobs)] <- 0
    wobs <- pmax(wobs, min.weight)
    wobs <- pmin(wobs, max.weight)
    if (equal.weights.as.null) {
      r <- range(wobs)
      if (r[2] - r[1] < 1e-15)
        weights <- NULL
    }
  }
  if (is.null(weights)) {
    o <- order(xobs)
    #lo <- lowess(x = xobs, y = yobs, f = span, iter = iterations -
    #               1L)
    lo <- loess(yobs ~ xobs, span = span,
                 degree = 1, parametric = FALSE, normalize = FALSE,
                 statistics = "approximate", surface = "direct",  # interpolate
                 cell = 0.01/span, iterations = iterations, trace.hat = "approximate",
                family = "gaussian")


    # loess(M[subset] ~ A[subset], degree = 1, span = span,
    #       #                family = family.loess)


    out$fitted[obs][o] <- lo$y
    out$residuals[obs] <- yobs - out$fitted[obs]
    out$mod <- lo
    out <<- out
    lo <<- lo
    return(out)
  }
  # if (min.weight > 0)
  #   nwobs <- nobs
  # else nwobs <- sum(wobs > 0)
  # if (nwobs < 4 + 1/span) {
  #   if (nwobs == 1L) {
  #     out$fitted[obs] <- yobs[wobs > 0]
  #     out$residuals[obs] <- yobs - out$fitted[obs]
  #   }
  #   else {
  #     fit <- lm.wfit(cbind(1, xobs), yobs, wobs)
  #     out$fitted[obs] <- fit$fitted
  #     out$residuals[obs] <- fit$residuals
  #   }
  #   return(out)
  # }
  # method <- match.arg(method, c("weightedLowess", "locfit",
  #                               "loess"))
  # switch(method, weightedLowess = {
  #   fit <- weightedLowess(x = xobs, y = yobs, weights = wobs,
  #                         span = span, iterations = iterations, npts = 200)
  #   out$fitted[obs] <- fit$fitted
  #   out$residuals[obs] <- fit$residuals
  #   out$mod <- fit
  #   fit <<- fit
  # }, locfit = {
  #   if (!requireNamespace("locfit", quietly = TRUE)) stop("locfit required but is not installed (or can't be loaded)")
  #   biweights <- rep(1, nobs)
  #   for (i in 1:iterations) {
  #     fit <- locfit::locfit(yobs ~ xobs, weights = wobs *
  #                             biweights, alpha = span, deg = 1)
  #     res <- residuals(fit, type = "raw")
  #     s <- median(abs(res))
  #     biweights <- pmax(1 - (res/(6 * s))^2, 0)^2
  #   }
  #   out$fitted[obs] <- fitted(fit)
  #   out$residuals[obs] <- res
  #   out$mod <- fit
  # }, loess = {
  #   oldopt <- options(warn = -1)
  #   on.exit(options(oldopt))
  #   bin <- 0.01
  #   fit <- loess(yobs ~ xobs, weights = wobs, span = span,
  #                degree = 1, parametric = FALSE, normalize = FALSE,
  #                statistics = "approximate", surface = "interpolate",
  #                cell = bin/span, iterations = iterations, trace.hat = "approximate")
  #   out$fitted[obs] <- fit$fitted
  #   out$residuals[obs] <- fit$residuals
  #   out$mod <- fit
  # })
  # print(fit)
  # out
}



################################################################################
################################################################################
##affy:: normalizeCyclicLoess

normalizeCyclicLoess2 <- function (x, weights = NULL, span = 0.7, iterations = 3, method = "fast", subset = NULL) {
  if (is.null(subset)) {
    subset <- 1:nrow(x)
  }

  x <- as.matrix(x)
  method <- match.arg(method, c("fast", "affy",
                                "pairs"))
  n <- ncol(x)
  if (method == "pairs") {
    for (k in 1:iterations) for (i in 1:(n - 1)) for (j in (i +
                                                            1):n) {
      m <- x[, j] - x[, i]
      a <- 0.5 * (x[, j] + x[, i])
      mod <- loessFit2(m[subset], a[subset], weights = weights, span = span, method = "loess") # "weightedLowess"
      f <- stats:::predict.loess(mod$mod, a) # cbind(m, a)
      x[, i] <- x[, i] + f/2
      x[, j] <- x[, j] - f/2
    }
  }
  if (method == "fast") {
    for (k in 1:iterations) {
      a <- rowMeans(x, na.rm = TRUE)
      for (i in 1:n) {
        m <- x[, i] - a
        mod <- loessFit2(m[subset], a[subset], weights = weights, span = span, method = "loess")#$fitted
        mod <<- mod
        f <- stats:::predict.loess(mod$mod, a)
        print(i)
        print(summary(f))
        x[, i] <- x[, i] - f
      }
    }
  }
  # if (method == "affy") {
  #   g <- nrow(x)
  #   for (k in 1:iterations) {
  #     adjustment <- matrix(0, g, n)
  #     for (i in 1:(n - 1)) for (j in (i + 1):n) {
  #       m <- x[, j] - x[, i]
  #       a <- 0.5 * (x[, j] + x[, i])
  #       mod <- loessFit2(m, a, weights = weights, span = span, method = "loess")#$fitted
  #       f <- stats:::predict.loess(mod$mod, cbind(m, a))
  #       adjustment[, j] <- adjustment[, j] + f
  #       adjustment[, i] <- adjustment[, i] - f
  #     }
  #     x <- x - adjustment/n
  #   }
  # }
  x
}

