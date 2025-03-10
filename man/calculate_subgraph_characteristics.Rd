% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tables-subgraph_characteristics.R
\name{calculate_subgraph_characteristics}
\alias{calculate_subgraph_characteristics}
\title{Generates a table with characteristics for each subgraph in a list.}
\usage{
calculate_subgraph_characteristics(
  S,
  fastalevel = TRUE,
  prototype = FALSE,
  file = NULL
)
}
\arguments{
\item{S}{list of subgraphs, where peptide and protein nodes are collapsed}

\item{fastalevel}{Are the subgraphs on fasta level?}

\item{prototype}{Are the subgraphs part of a prototype list?}

\item{file}{where to save the table.}
}
\value{
table
}
\description{
Generates a table with characteristics for each subgraph in a list.
}
\examples{
# TODO
}
