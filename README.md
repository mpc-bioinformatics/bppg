
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bppg

<!-- badges: start -->
<!-- badges: end -->

The bppg package contains functionality to create and characterize 
bipartite graphs that model the relationship between peptides and proteins 
in bottom-up proteomics. With these graphs, protein ratios (fold change 
between two sample groups) are calculated from the respective measured
peptide ratios. The main aim is to make use of quantitative information
contained in shared peptides and making it possible to quantify proteins
without shared peptides.

## Installation

The current release version can be installed from Bioconductor usin the following 
code:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("bppg")
```


You can also install the development version of bppg from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("mpc-bioinformatics/bppg")
```

## Usage

For details on the usage of bppg please see the vignette.


## Publication

The quantification algorithm in bppg is build on the PhD thesis

"Improvement of protein quantification for proteins with shared peptides by 
using bipartite peptide-protein graphs"

http://dx.doi.org/10.17877/DE290R-25361


## Funding

The development of bppg is funded by the German Research Foundation (DFG, 
project number 532401634), the German Network for Bioinformatics Infrastructure 
(de.NBI & ELIXIR-DE, Federal Ministry of Research, Technology and Space, grant 
number W-de.NBI-005) and the Core Unit for Bioinformatics of the Medical Faculty
of the Ruhr University Bochum (CUBiMed.RUB). 





