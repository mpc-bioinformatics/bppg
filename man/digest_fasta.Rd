% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers-Digest.R
\name{digest_fasta}
\alias{digest_fasta}
\title{In silico tryptic digestion of whole FASTA file}
\usage{
digest_fasta(fasta, missed_cleavages = 2, min_aa = 6, max_aa = 50, ...)
}
\arguments{
\item{fasta}{List of protein sequences (e.g., imported FASTA file by seqinr::read.fasta).}

\item{missed_cleavages}{Maximal number of missed cleavages.}

\item{min_aa}{Minimal number of amino acids (set to 0 for no filtering).}

\item{max_aa}{Maximal number of amino acids (set to Inf for no filtering).}

\item{...}{Additional arguments for Digest2().}
}
\value{
List of vectors of peptide sequences, filtered for minimal and maximal number of
amino acids.
}
\description{
In silico tryptic digestion of whole FASTA file
}
\examples{
library(seqinr)
file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
res <- digest_fasta(fasta)

}
