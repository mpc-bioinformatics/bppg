library(BBmisc)    # for collapse()
library(seqinr)    # for reading in fasta files
library(limma)     # for strsplit2()
library(pbapply)   # for progress bars for apply functions
library(igraph)    # for graph functionality
library(openxlsx)  # for reading and writing xlsx files
library(bppg)      # for generating the graphs
library(ggplot2)   # for plotting
library(data.table)
library(stringr)   # for string handling

.digest2 <- function(sequence,
    enzyme = "trypsin",
    missed = 0,
    warn = TRUE,
    remove_initial_M = FALSE) {
    seq_vector <- strsplit(sequence, split = "")[[1]]
    end_position <- length(seq_vector)
    if (enzyme == "trypsin") {
        if (seq_vector[end_position] == "K" | seq_vector[end_position] == "R") {
            seq_vector[end_position] <- "!"
            seq_string <- paste(seq_vector, collapse = "")
        }
        else {
            seq_string <- sequence
        }
        seq_string <- gsub("KP", "!P", seq_string)
        seq_string <- gsub("RP", "!P", seq_string)
        seq_vector <- strsplit(seq_string, split = "")[[1]]
        stop <- grep("K|R", seq_vector)
        start <- stop + 1
    }
    if (enzyme == "trypsin.strict") {
        if (seq_vector[end_position] == "K" | seq_vector[end_position] == "R") {
            seq_vector[end_position] <- "!"
            seq_string <- paste(seq_vector, collapse = "")
        }
        else {
            seq_string <- sequence
        }
        seq_vector <- strsplit(seq_string, split = "")[[1]]
        stop <- grep("K|R", seq_vector)
        start <- stop + 1
    }
    if (enzyme != "trypsin" & enzyme != "trypsin.strict")
        stop("undefined enzyme, defined enzymes are trypsin, trypsin.strict")
    if (length(stop) == 0) {
        if (warn) warning("sequence does not contain cleavage sites")
        return(data.frame(sequence = sequence, start = 1,
                stop = nchar(sequence), mc = 0))
    }

    if (missed > length(stop)) {
        if (warn) warning("number of specified missed cleavages is greater than
            the maximum possible")
    }

    cleave <- function(sequence, start, stop, misses) {
        peptide <- substring(sequence, start, stop)
        mc <- rep(misses, times = length(peptide))
        data.frame(sequence = peptide, start, stop, mc,
            stringsAsFactors = FALSE)
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
        y2 <- y2[substr(y2$sequence, 1, 1) == "M", ] ## is first amino acid M?

        if (nrow(y2) > 0) {
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








#' In silico tryptic digestion of whole FASTA file.
#'
#' @param fasta              \strong{list of vector of characters} \cr
#'                           A fasta file, already read into R by
#'                           [seqinr::read.fasta()].
#' @param missed_cleavages   \strong{integer} \cr
#'                           The maximal number of missed cleavages.
#' @param min_aa             \strong{integer} \cr
#'                           The minimal number of amino acids
#'                           (set to 0 for no filtering).
#' @param max_aa             \strong{integer} \cr
#'                           The maximal number of amino acids
#'                           (set to Inf for no filtering).
#' @param ...                Additional arguments for [.digest2()].
#'
#' @return List of vectors of peptide sequences, filtered for minimal
#'         and maximal number of amino acids.
#' @export
#' 
#'
#' @seealso [.digest2()]
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' res <- digestFASTA(fasta)
#'
## TODO USE https://bioconductor.org/packages/3.22/bioc/html/cleaver.html
# cleave("LAAGKVEDSD", enzym = "trypsin", missedCleavages = 0:2)
## by Sebastian Gibb ehemals bei Laurent Gatto
digestFASTA <- function(fasta,
    missed_cleavages = 2,
    min_aa = 6,
    max_aa = 50,
    ...)  {

    digested_proteins <- pbapply::pblapply(fasta, function(x) {
        sequ <- x
        class(sequ) <- NULL
        y <- try({
            .digest2(sequ, missed = missed_cleavages, warn = FALSE,
                remove_initial_M = TRUE, ...)})
        ind <- nchar(as.character(y$sequence)) >= min_aa &
            nchar(as.character(y$sequence)) <= max_aa
        as.character(y$sequence[ind])
    })

    return(digested_proteins)
}


#' Generate edgelist from list of in silico digested proteins.
#'
#' @param digested_proteins   \strong{list of vector of characters} \cr
#'                            The output from [digestFASTA()] 
#'                            (List of vectors of peptide sequences)
#' @param prot_origin         \strong{vector of characters} \cr
#'                            origin of the protein (e.g. organism, 
#'                            spike-in/background etc)
#'
#' @return An edgelist.
#' @export
#'
#' @seealso [digestFASTA()]
#'
#' @examples
#' library(seqinr)
#' file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
#' fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
#' digested_proteins <- digestFASTA(fasta)
#' edgelist <- generateEdgelist(digested_proteins)
#'
#'

generateEdgelist <- function(digested_proteins, prot_origin = NULL) {
    ## calculate necessary number of edges by counting the peptides belonging to 
    ## each protein
    mat_length <- sum(lengths(digested_proteins))

    ## generate empty edge matrix of size (#edges)x2
    if (is.null(prot_origin)) {
        edgelist <- matrix(nrow = mat_length, ncol = 2)
    } else {
        edgelist <- matrix(nrow = mat_length, ncol = 3)
    }

    ## add progress bar to loop
    number_of_iterations <- length(digested_proteins)
    pb <- pbapply::startpb(0, length(digested_proteins))
    on.exit(pbapply::closepb(pb))

    ## add an entry to the edge matrix for each peptide-protein relation in the
    ## digested_proteins matrix
    current_row <- 1
    for (i in 1:length(digested_proteins)){ ## TODO VAPPLY
        if (length(digested_proteins[[i]]) != 0) {
            for (j in 1:length(digested_proteins[[i]])){
                edgelist[current_row, 1] <- names(digested_proteins)[[i]]
                edgelist[current_row, 2] <- digested_proteins[[i]][[j]]

                if (!is.null(prot_origin)) {
                    edgelist[current_row, 3] <- prot_origin[[i]]
                }
                current_row <- current_row + 1
            }
            pbapply::setpb(pb, i)
        }
    }

    #progress bar command
    invisible(NULL)

    #find and remove duplicate rows that would lead to duplicate edges
    duplicate_rows <- duplicated(edgelist, margin = 1)
    edgelist <- edgelist[!duplicate_rows, ]

    edgelist <- as.data.frame(edgelist)
    if (is.null(prot_origin)) {
        colnames(edgelist) <- c("protein", "peptide")
    } else {
        colnames(edgelist) <- c("protein", "peptide", "prot_origin")
    }

    return(edgelist)
}





fasta <- seqinr::read.fasta(file = "/home/wenzelju/Documents/Code/test_data/uniprotkb_MYH.fasta",
                             seqtype = "AA", as.string = TRUE)
names(fasta) <- limma::strsplit2(names(fasta), "\\|")[, 2]
digestedEdgelist <- bppg::digestFASTA(fasta)

system.time(replicate(1000, digested <- digestFASTA(fasta)))
#    user  system elapsed   mit cleave
#  80.384   0.721  86.290   33.464   0.700  37.367
system.time(replicate(1000, digested <- bppg::digestFASTA(fasta)))
# 71.065   0.665  76.161

digested <- digestFASTA(fasta)
tmp_edgelist <- generateEdgelist(digested)
system.time(replicate(1000, generateEdgelist(digested)))
# 58.171   0.604  62.594 

generics::setdiff(digestedEdgelist, tmp_edgelist)


# testthat files:
edgelist <- readRDS(testthat::test_path("testfiles/edgelist_test.rds"))

file <- system.file("extdata", "uniprot_test.fasta", package = "bppg")
fasta <- seqinr::read.fasta(file = file, seqtype = "AA", as.string = TRUE)
res <- bppg::digestFASTA(fasta)

generics::setdiff(res$peptide, edgelist$peptide)
