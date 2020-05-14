#!/mnt/aspnas/sw/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop('usage: rundada2TAX.R <seqtab.rds> <taxonomy database> <output dir>')
}

seq.tab.rds <- args[1]
tax_db <- args[2]
out <- args[3]

library(dada2)
library(here)

seq.tab <- readRDS(seq.tab.rds)

taxa <- assignTaxonomy(seq.tab, tax_db, tryRC = TRUE, multithread = TRUE)

saveRDS(taxa, here(paste0(out, "taxa.rds")))