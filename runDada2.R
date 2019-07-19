#!/mnt/aspnas/sw/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop('usage: runDada2.R <Forward filtered> <Reverse filtered>')
}

fw <- args[1]
rv <- args[2]

library(dada2); packageVersion("dada2")
library(here)
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
errF <- learnErrors(filtFs, multithread=TRUE) #nbases=1e8, 
# Learn reverse error rates
errR <- learnErrors(filtRs, multithread=TRUE) #nbases=1e8, 
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
dir.create(here("data/R_files/dada2"))
saveRDS(seqtab, here("data/R_files/dada2/seqtab.rds")) # CHANGE ME to where you want sequence table saved