#!/usr/bin/env Rscript

print(snakemake@params$maxN)

if ()

require(dada2, quietly = TRUE)

# Parameters
if (is.null(snakemake@params$maxN)) {snakemake@params$maxN <- 0}
if (is.null(snakemake@params$truncQ)) {snakemake@params$truncQ <- 2}
if (is.null(snakemake@params$truncLen)) {snakemake@params$truncLen <- 0}
if (is.null(snakemake@params$maxEE)) {snakemake@params$maxEE <- Inf}
if (is.null(snakemake@params$rm_phix)) {snakemake@params$rm_phix <- TRUE}
if (is.null(snakemake@params$minLen)) {snakemake@params$minLen <- 0}
# Performance
threads <- snakemake@threads

filterAndTrim(fwd = snakemake@input$R1, filt = snakemake@output$R1,
              rev = snakemake@input$R2, filt.rev = snakemake@output$R2,
              maxN = snakemake@params$maxN, truncQ = snakemake@params$truncQ,
              maxEE = snakemake@params$maxEE, rm.phix = snakemake@params$rm_phix,
              minLen = snakemake@params$minLen, multithread = threads)
