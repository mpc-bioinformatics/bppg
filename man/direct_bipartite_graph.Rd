% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers-isomorphisms.R
\name{direct_bipartite_graph}
\alias{direct_bipartite_graph}
\title{Transform a bipartite graph into a directed graph}
\usage{
direct_bipartite_graph(bip_graph, from_type = FALSE)
}
\arguments{
\item{bip_graph}{a bipartite graph}

\item{from_type}{determines if protein or peptide nodes are the "from" nodes}
}
\value{
a bipartite graph that is know directed
}
\description{
Transform a bipartite graph into a directed graph
}
