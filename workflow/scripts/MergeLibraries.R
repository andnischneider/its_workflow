#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#The idea here is that this script does not have a set number of
#arguments, but depends on the number of libraries that have been
#processed and will now be merged.
#Threads argument is not needed (?)
#Usage: MergeLibraries.R <output_dir> <input_dir1> <input_dir2> ... <input_dirn>

out <- args[1]
library(dada2); packageVersion("dada2")
library(Biostrings)

nlibs <- length(args)-1
lib_list <- list()

for (i in 1:nlibs) {
  lib_list[[i]] <- paste0(args[i+1], "/seqtab.nc_itsx_clean.rds")
}

#Conveniently, dada2 has a function to merge any number of sequence tables!

seqtab.nc_all <- mergeSequenceTables(tables = lib_list)
#Save this 
saveRDS(seqtab.nc_all, paste0(out, "/seqtab.nc_all.rds"))

#Export merged sequences for Swarm

dna_clean_sum <- DNAStringSet(colnames(seqtab.nc_all[,colSums(seqtab.nc_all)>0]))
names(dna_clean_sum) <- paste0("ASV", 1:length(dna_clean_sum))
names(dna_clean_sum) <- paste0(names(dna_clean_sum), ";size=", colSums(seqtab.nc_all[,colSums(seqtab.nc_all)>0]))
writeXStringSet(dna_clean_sum, file = paste0(out, "/seqs_sum.fasta"))