#!/usr/bin/env Rscript

require(dada2)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  stop('usage: cutITS.R <input R1> <input R2> <output R1> <output R2> <threads>')
}

R1 <- args[1]
R2 <- args[2]
R1_o <- args[3]
R2_o <- args[4]
threads <- args[5]

filterAndTrim(R1, R1_o, R2, R2_o, maxN = 0, multithread = TRUE)
