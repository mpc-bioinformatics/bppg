

### preprocess the quantitative data of D2

library(BBmisc)    # for collapse()
library(seqinr)    # for reading in fasta files
library(limma)     # for strsplit2()
library(pbapply)   # for progress bars for apply functions
library(igraph)
library(bppg)


################################################################################
################################################################################
### read in quantitative peptide data ##########################################
DATA <- bppg::read_MQ_peptidetable("C:/Users/schorkka/UNI/Promotion/promotion_project/data/D2_without_isoforms/D2_quant/peptides_D2.txt",
                                   LFQ = FALSE, remove_contaminants = FALSE)

validvalues <- rowSums(!is.na(DATA[,-1]))
barplot(table(validvalues))

### only keep peptides with at least one valid value
### TODO: Einstellung in read_MQ_peptidetable?
DATA <- DATA[validvalues >= 1, ]

### read in pre-calculated edgelist ############################################
edgelist <- read.table("C:/Users/schorkka/UNI/Promotion/promotion_project/data/D2_without_isoforms/D2_fasta/preprocessed/edgelist_collprot_.txt",
                       sep = "\t", header = TRUE)
edgelist_filtered <- edgelist[edgelist[,2] %in% DATA$Sequence, ]


### use edgelist to assign proteins to the peptide data and save it
for (i in 1:nrow(DATA)) {
  pep_tmp <- DATA$Sequence[i]
  prot_tmp <- edgelist_filtered$protein[edgelist_filtered$peptide == pep_tmp]
  prot_tmp <-
}

proteins <-




#prepare for normalization
id <- DATA[,1]
intensities <- DATA[,-1]
intensities[intensities == 0] <- NA

#normalize Intensities
# TODO: auch Median, Quantilsnormalisierung und LFQ-Normalisierung von MaxQuant erlauben?
intensities <- 2^limma::normalizeBetweenArrays(log2(intensities), method = "cyclicloess")




