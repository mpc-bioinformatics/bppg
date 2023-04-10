

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
DATA <- bppg::read_MQ_peptidetable("C:/Users/schorkka/UNI/Promotion/promotion_project/data/D2_without_isoforms/D2_quant/MQ/peptides_D2.txt",
                                   LFQ = FALSE, remove_contaminants = FALSE, zeroToNA = TRUE,
                                   remove_empty_rows = TRUE)

fasta1 <- seqinr::read.fasta(file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/FASTAs/uniprot_proteome_saccharomycescerevisiae_UP000002311_20230220_v202205.fasta",
                             seqtype = "AA", as.string = TRUE)
names(fasta1) <- limma::strsplit2(names(fasta1), "\\|")[,2]
fasta2 <- seqinr::read.fasta(file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/FASTAs/ups1-ups2-sequences.fasta",
                             seqtype = "AA", as.string = TRUE)
names(fasta2) <- limma::strsplit2(names(fasta2), "\\|")[,1]
fasta3 <- seqinr::read.fasta(file = "C:/Users/schorkka/UNI/Promotion/promotion_project/data/FASTAs/MaxQuant_contaminants_20220220.fasta",
                             seqtype = "AA", as.string = TRUE)
names(fasta3) <- paste0("CON_", names(fasta3))
fasta <- c(fasta1, fasta2, fasta3)




