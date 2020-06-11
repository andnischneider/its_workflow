#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop('usage: processSwarm.R <swarm_output_dir> <swarm_input_dir> <out_dir>')
}

swarm_res <- args[1]
swarm_in <- args[2]
outdir <- args[3]


library(dada2); packageVersion("dada2")
library(Biostrings)
library(dplyr)
library(readr)
library(stringr)
library(reshape2)

#Import sequences and merged seqtab
seqtab.nc_all <- readRDS(paste0(swarm_in, "/seqtab.nc_all.rds"))
seqs_merged <- readDNAStringSet(paste0(swarm_in, "/seqs_sum.fasta"))

#Function to import swarm clusters
clustfact <- function(file) {
  clus <- readr::read_tsv(file, col_names = FALSE)
  clus <- stringr::str_split_fixed(clus$X1, pattern = " ", n = max(stringr::str_count(clus$X1, "ASV")))
  rownames(clus) <- paste0("cluster", 1:nrow(clus))
  clus <- melt(as.matrix(clus))[,-2]
  clus <- clus[grepl("ASV", clus$value),]
  clus$value <- gsub(";size=\\d+", "", clus$value)
  #clus2 <- as.factor(clus$Var1)
  #names(clus2) <- clus$value
  colnames(clus) <- c("Cluster", "ASV")
  #clus$Cluster <- as.character(clus$Cluster)
  return(clus)
}
clus_swarm <- clustfact(paste0(swarm_res, "/results.txt"))

seqtab_BOTH_sw <- as.data.frame(cbind(Cluster=clus_swarm$Cluster[match(gsub(";size=\\d+", "", names(seqs_merged)), clus_swarm$ASV)], t(seqtab.nc_all)))

seqtab_BOTH_sw <- group_by(seqtab_BOTH_sw, Cluster) %>% 
  summarise_each(funs(sum)) 
seqtab_BOTH_sw$Cluster <- paste0("Cluster", seqtab_BOTH_sw$Cluster)
seqtab_BOTH_sw2 <- seqtab_BOTH_sw
rownames(seqtab_BOTH_sw2) <- seqtab_BOTH_sw$Cluster
seqtab_BOTH_sw2 <- data.matrix(seqtab_BOTH_sw2[,-1])
#This is now ready to proceed and assign taxonomy to all the Swarm clusters
seeds_swarm <- readDNAStringSet(paste0(swarm_res, "/seeds.fasta"))
names(seeds_swarm) <- seqtab_BOTH_sw$Cluster 

seqtab_BOTH_sw3 <- seqtab_BOTH_sw2
rownames(seqtab_BOTH_sw3) <- as.character(seeds_swarm)
seqtab_BOTH_sw3 <- t(seqtab_BOTH_sw3)
#THis can now be used to assign taxonomy on.
saveRDS(seqtab_BOTH_sw3, paste0(outdir, "/seqtab_final.rds"))







