#!/usr/bin/env Rscript

require(dada2, quietly=TRUE)

#Parameters
if (is.null(snakemake@params$maxN)) {snakemake@params$maxN <- 0}
if (is.null(snakemake@params$truncQ)) {snakemake@params$truncQ <- 2}
if (is.null(snakemake@params$truncLen)) {snakemake@params$truncLen <- 0}
if (is.null(snakemake@params$maxEE_R1)) {snakemake@params$maxEE_R1 <- Inf}
if (is.null(snakemake@params$maxEE_R2)) {snakemake@params$maxEE_R2 <- Inf}
if (is.null(snakemake@params$rm_phix)) {snakemake@params$rm_phix <- TRUE}
if (is.null(snakemake@params$minLen)) {snakemake@params$minLen <- 0}

maxEE <- c(snakemake@params$maxEE_R1, snakemake@params$maxEE_R2)

# Performance
threads <- snakemake@threads

out <- filterAndTrim(fwd = snakemake@input$R1, filt = snakemake@output$R1,
              rev = snakemake@input$R2, filt.rev = snakemake@output$R2,
              maxN = snakemake@params$maxN, truncQ = snakemake@params$truncQ,
              maxEE = maxEE, rm.phix = snakemake@params$rm_phix,
              minLen = snakemake@params$minLen, multithread = threads,
              verbose=TRUE)

saveRDS(out, snakemake@output$rds)
