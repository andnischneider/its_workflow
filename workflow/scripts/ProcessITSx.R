#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

itsx_in <- snakemake@input$its_in
itsx_out <- snakemake@input$its_out
sq.nc <- snakemake@input$seqtab
out <- snakemake@output$seqtab
out_mock <- snakemake@output$mock
threads <- snakemake@threads
mock <- snakemake@params$mock

mock_samples <- strsplit(mock, ",")[[1]]

library(dada2); packageVersion("dada2")
library(dplyr)
library(Biostrings)

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

if (length(intersect(colnames(seqtab.nc3), mock_samples)) > 0) {
  mock_cols <- append(c("sequence"), intersect(colnames(seqtab.nc3), mock_samples))
  mock <- seqtab.nc3 %>% select(all_of(mock_cols))
  seqtab.nc4 <- seqtab.nc3 %>% select(!any_of(mock_cols))
  saveRDS(mock, out_mock)
} else {
  seqtab.nc4 <- seqtab.nc3[,colnames(seqtab.nc3)!="sequence"]
}

rownames(seqtab.nc4) <- seqtab.nc3$sequence
seqtab.nc4 <- data.matrix(t(seqtab.nc4))

saveRDS(seqtab.nc4, out)
