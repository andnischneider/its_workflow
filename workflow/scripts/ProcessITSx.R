#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  stop('usage: processITSx.R <itsx_input> <itsx_output> <seqtab.nc> <Output dir> <threads>')
}

itsx_in <- args[1]
itsx_out <- args[2]
sq.nc <- args[3]
out <- args[4]
threads <- args[5]

library(dada2); packageVersion("dada2")
library(here)
library(dplyr)

dna <- readDNAStringSet(itsx_in)
dna_itsx <- as.character(readDNAStringSet(itsx_out))
seqtab.nc <- readRDS(sq.nc)
colnames(seqtab.nc) <- names(dna)
seqtab.nc <- seqtab.nc[,colnames(seqtab.nc)%in%names(dna_itsx)]

colnames(seqtab.nc) <- dna_itsx
seqtab.nc <- t(seqtab.nc)

#Remove sequences below 50bp
if (any(nchar(getSequences(t(seqtab.nc)))<50)) {
  dna_itsx <- dna_itsx[!(nchar(dna_itsx)<50)]
  seqtab.nc <- seqtab.nc[nchar(rownames(seqtab.nc))>50,]
}


#At this point we can summarise all identical sequences 
seqtab.nc2 <- cbind.data.frame(sequence=rownames(seqtab.nc), seqtab.nc)
seqtab.nc3 <- group_by(seqtab.nc2, sequence) %>% 
  summarise_each(funs(sum)) 
seqtab.nc4 <- seqtab.nc3[,-c(1,2)]
rownames(seqtab.nc4) <- seqtab.nc3$sequence
seqtab.nc4 <- data.matrix(t(seqtab.nc4))

saveRDS(seqtab.nc4, here(paste0(out, "/seqtab.nc_itsx_clean.rds")))