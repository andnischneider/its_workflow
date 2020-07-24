#!/usr/bin/env Rscript

library(dada2); packageVersion("dada2")
library(Biostrings)
library(RcppParallel)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  stop('usage: runDada2.R <Forward filtered> <Reverse filtered> <Output dir> <threads>')
}

fw <- args[1]
rv <- args[2]
out <- args[3]
threads <- args[4]

multithreading <- FALSE

if (threads > 1) {
  setThreadOptions(numThreads = as.integer(threads))
  multithreading <- TRUE
}

# File parsing
filtpathF <- fw # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- rv # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, multithread=multithreading) #nbases=1e8,
# Learn reverse error rates
errR <- learnErrors(filtRs, multithread=multithreading) #nbases=1e8,
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

#Additional step to also collect dada info
dadas_f <- vector("list", length(sample.names))
names(dadas_f) <- sample.names
dadas_r <- vector("list", length(sample.names))
names(dadas_r) <- sample.names

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=multithreading, pool = TRUE)
  dadas_f[[sam]] <- ddF
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=multithreading, pool = TRUE)
  dadas_r[[sam]] <- ddR
  merger <- mergePairs(ddF, derepF, ddR, derepR, maxMismatch = 1)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
seqtab.nc <- removeBimeraDenovo(seqtab, method = "consensus",
                                multithread = TRUE, verbose = TRUE)

saveRDS(dadas_f, paste0(out, "/dada_f.rds"))
saveRDS(dadas_r, paste0(out, "/dada_r.rds"))
saveRDS(mergers, paste0(out, "/mergers.rds"))
saveRDS(seqtab, paste0(out, "/seqtab.rds"))
saveRDS(seqtab.nc, paste0(out, "/seqtab_nc.rds"))

#write out sequences as fasta file for ITSx
dna <- DNAStringSet(colnames(seqtab.nc))
names(dna) <- paste0("ASV", seq(ncol(seqtab.nc)))
writeXStringSet(dna, paste0(out, "/seqs.fasta"))

