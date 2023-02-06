
############### install R packages needed for running LFMM on IFB

.libPaths("/shared/home/baotram/R/x86_64-pc-linux-gnu-library/4.0")

## install packages
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("devtools")) install.packages("devtools")
if (!require("optparse")) install.packages("optparse")
if (!require("Biostrings")) BiocManager::install("Biostrings")
if (!require("lfmm")) devtools::install_github("bcm-uga/lfmm",
                                               upgrade = "never",
                                               dependencies = T)
if (!require("vegan")) install.packages("vegan")
if (!require("BEDMatrix")) install.packages("BEDMatrix")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("data.table")) install.packages("data.table")
if (!require("RSpectra")) install.packages("RSpectra")
if (!require("Rcpp")) install.packages("Rcpp")
if (!require("pcadapt")) install.packages("pcadapt")
if (!require("qvalue")) install.packages("qvalue")
if (!require("qvalue")) install.packages("gaston")