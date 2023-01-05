#### modified version of OrgMassSpecR::Digest
#### - delete functionality to calculate peptide masses & enzymes other than trypsin
#### - interpret "missed" argument as maximum number of allowed missed cleavages and
####   only warn if nr or missed cleavages is not possible

#' Digestion of a single protein sequence.
#'
#' @param sequence Protein sequence as character.
#' @param enzyme "trypsin" (does not cut before proline) or "trypsin.strict" ().
#' @param missed Maximal number of missed cleavages.
#' @param warn Print out warnings, e.g. if a protein has no cleavage site.
#'
#' @return vector of peptides
#' @export
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "2020_01_31_proteome_S_cerevisae.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#'
#' digested_proteins <- Digest2(fasta[[1]])


Digest2 <- function (sequence, enzyme = "trypsin", missed = 0, warn = TRUE, remove_initial_M = FALSE) {
  seq_vector <- strsplit(sequence, split = "")[[1]]
  end_position <- length(seq_vector)
  if (enzyme == "trypsin") {
    if (seq_vector[end_position] == "K" | seq_vector[end_position] ==
        "R") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    }
    else seq_string <- sequence
    seq_string <- gsub("KP", "!P", seq_string)
    seq_string <- gsub("RP", "!P", seq_string)
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("K|R", seq_vector)
    start <- stop + 1
  }
  if (enzyme == "trypsin.strict") {
    if (seq_vector[end_position] == "K" | seq_vector[end_position] ==
        "R") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    }
    else seq_string <- sequence
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("K|R", seq_vector)
    start <- stop + 1
  }
  if (enzyme != "trypsin" & enzyme != "trypsin.strict")
    stop("undefined enzyme, defined enzymes are trypsin, trypsin.strict")
  if (length(stop) == 0) {
    if (warn) warning("sequence does not contain cleavage sites")
    return(data.frame(sequence = sequence, start = 1, stop = nchar(sequence), mc = 0))
  }

  if (missed > length(stop)) {
    if (warn) warning("number of specified missed cleavages is greater than the maximum possible")
  }

  cleave <- function(sequence, start, stop, misses) {
    peptide <- substring(sequence, start, stop)
    mc <- rep(misses, times = length(peptide))
    result <- data.frame(sequence = peptide, start, stop, mc, stringsAsFactors = FALSE)
    return(result)
  }
  stop_ <- stop
  start <- c(1, start)
  stop <- c(stop, end_position)
  results <- cleave(sequence, start, stop, 0)
  if (missed > 0) {
    for (i in 1:min(missed, length(stop_))) {
      start_tmp <- start[1:(length(start) - i)]
      stop_tmp <- stop[(1 + i):length(stop)]
      peptide <- cleave(sequence, start_tmp, stop_tmp, i)
      results <- rbind(results, peptide)
    }
  }


  if (remove_initial_M) {

    y2 <- results[results$start == 1,] ## there should be at least 1
    y2 <- y2[substr(y2$sequence, 1, 1) == "M"] ## is first amino acid M?

    if (length(y2) > 0) {
      y2$sequence <- substr(y2$sequence, 2, nchar(y2$sequence))
      y2$start <- 2
      results <- rbind(results, y2)
    }

  }

   return(results)
}


# License: BSD_2_clause + file LICENSE
# Copyright (c) 2011-2017, Nathan Dodder
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#   Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
#   Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in
# the documentation and/or other materials provided with the
# distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.








#' In silico tryptic digestion of whole FASTA file
#'
#' @param fasta List of protein sequences (e.g., imported FASTA file by seqinr::read.fasta).
#' @param missed_cleavages Maximal number of missed cleavages.
#' @param min_aa Minimal number of amino acids (set to 0 for no filtering).
#' @param max_aa Maximal number of amino acids (set to Inf for no filtering).
#' @param ... Additional arguments for Digest2().
#'
#' @return List of vectors of peptide sequences, filtered for minimal and maximal number of
#' amino acids.
#' @export
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "2020_01_31_proteome_S_cerevisae.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' res <- digest_fasta(fasta)
#'
digest_fasta <- function(fasta, missed_cleavages = 2, min_aa = 6, max_aa = 50, ...)  {

  digested_proteins <- pbapply::pblapply(fasta, function(x) {
    sequ <- x
    class(sequ) <- NULL
    y <- try({Digest2(sequ, missed = missed_cleavages, warn = FALSE, remove_initial_M = TRUE, ...)})
    ind <- nchar(as.character(y$sequence)) >= min_aa & nchar(as.character(y$sequence)) <= max_aa
    return(as.character(y$sequence[ind]))
  })

  return(digested_proteins)
}






