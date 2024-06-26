% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_graphs_from_FASTA.R
\name{generate_graphs_from_FASTA}
\alias{generate_graphs_from_FASTA}
\title{Generate graphs from a FASTA file}
\usage{
generate_graphs_from_FASTA(
  fasta,
  collapse_protein_nodes = TRUE,
  collapse_peptide_nodes = TRUE,
  result_path = NULL,
  suffix = NULL,
  save_intermediate = FALSE,
  prot_origin = NULL,
  ...
)
}
\arguments{
\item{fasta}{fasta file, already read into R by seqinr::read.fasta}

\item{collapse_protein_nodes}{collapse protein nodes?}

\item{collapse_peptide_nodes}{collapse peptide nodes?}

\item{result_path}{path whereresults are saved. If NULL, results are not saved}

\item{suffix}{suffix for saving results}

\item{save_intermediate}{Save intermediate results?}

\item{...}{additional arguments to bppg::digest_fasta()}
}
\value{
subgraphs (i.e. connected components) from the graph generated from the FASTA file.
}
\description{
Generate graphs from a FASTA file
}
\examples{
library(seqinr)
file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
graphs <- bppg::generate_graphs_from_FASTA(fasta)

}
