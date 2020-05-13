#!/usr/bin/env Rscript

require(dada2, quietly = TRUE)
# IO
R1 <- snakemake@input$R1
R2 <- snakemake@input$R2
R1_o <- snakemake@output$R1
R2_o <- snakemake@output$R2
# Parameters
maxN <- snakemake@params$maxN
truncQ <- snakemake@params$truncQ
truncLen <- snakemake@params$truncLen
maxEE <- snakemake@params$maxEE
rm_phix <- snakemake@params$rm_phix
minLen <- snakemake@params$minLen
# Performance
threads <- snakemake@threads

filterAndTrim(fwd = R1, filt = R1_o, rev = R2, filt.rev = R2_o,
              maxN = maxN, truncQ = truncQ, maxEE = maxEE, rm.phix = rm_phix,
              minLen = minLen, multithread = threads)
