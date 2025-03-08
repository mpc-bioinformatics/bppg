% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers-small_helper_functions.R
\name{geom_mean}
\alias{geom_mean}
\title{Geometric mean}
\usage{
geom_mean(x, useprod = FALSE)
}
\arguments{
\item{x}{vector with numbers}

\item{useprod}{if TRUE, prod(x)^(1/n) will be calculated, otherwise exp(mean(log(x)))}
}
\value{
geometric mean of the provided data points
}
\description{
Geometric mean
}
\examples{
# TODO
}
