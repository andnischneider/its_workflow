#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop('usage: rundada2TAX.R <seqtab.rds> <taxonomy database> <output dir>')
}

seq.tab.rds <- args[1]
tax_db <- args[2]
out <- args[3]
threads <- args[4]

library(dada2)
library(RcppParallel)

if (threads > 1) {
  setThreadOptions(numThreads = as.integer(threads))
  multithreading <- TRUE
}

seq.tab <- readRDS(seq.tab.rds)

taxa <- assignTaxonomy(seq.tab, tax_db, tryRC = TRUE, multithread = multithreading)

saveRDS(taxa, paste0(out, "taxa.rds"))