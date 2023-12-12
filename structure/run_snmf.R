library(BEDMatrix)
#!/usr/bin/env Rscript --vanilla
library(LEA)
library(vcfR)
library(tidyverse)
library(data.table)
library(Rsamtools)

## SNMF
snmf_dir <- "/shared/projects/most_kmer/afterkiss/snmf"
dir.create(snmf_dir, showWarnings = F)
geno_file <- file.path(snmf_dir, "aftr_kmer.geno")
# write.geno(submat, geno_file)

aftr_snmf_kmer <- snmf(geno_file, K = 1:10, CPU = 48, repetitions = 10, project = "new", entropy = T, seed = 987654321)