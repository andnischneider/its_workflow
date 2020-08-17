#This script processes the ITS amplicon and RNA
#data for easier use in Figure scripts.

#This has to be run before any of the Figure scripts

library(here)
library(DESeq2)
source(here("../src/UPSCb-common/src/R/featureSelection.R"))

#FIGURE 2
#Import meta files and add filtering factor
meta_rna_r <- read.csv2(here("../doc/RNA/Roots.csv"), stringsAsFactors = FALSE)
meta_rna_r$Group <- paste(meta_rna_r$treatment, meta_rna_r$date, sep = ".")

meta_rna_n <- read.csv2(here("../doc/RNA/Needles.csv"), stringsAsFactors = FALSE)
meta_rna_n$Group <- paste(meta_rna_n$treatment, meta_rna_n$date, sep = ".")

#Import raw counts (RNA WORKFLOW NEEDS TO BE RUN FIRST AND RESULT FILES COPIED INTO data/RNA folder)
counts_raw_rna <- read.delim(here("../data/RNA/Assembly_2012.raw.tsv"), row.names = 1)
#Import taxonomic annotations and select only the ones identified as fungi
annot_tax <- read.delim(here("../data/RNA/gene_taxonomy.tsv"))
annot_tax <- annot_tax[annot_tax$kingdom=="Fungi",]
counts_raw_rna <- counts_raw_rna[rownames(counts_raw_rna) %in% annot_tax$gene,]
#Extract roots and needles
counts_raw_rna_roots <- data.matrix(counts_raw_rna[,meta_rna_r$SciLifeID])
counts_raw_rna_needl <- data.matrix(counts_raw_rna[,meta_rna_n$SciLifeID])

#Create factors for filtering
ccc_r <- as.factor(meta_rna_r$Group)
names(ccc_r) <- meta_rna_r$SciLifeID
ccc_n <- as.factor(meta_rna_n$Group)
names(ccc_n) <- meta_rna_n$SciLifeID

#Filter roots 
dim(counts_raw_rna_roots)
counts_raw_rna_roots_f <- counts_raw_rna_roots[featureSelect(counts_raw_rna_roots, ccc_r, 5, 2),]
dim(counts_raw_rna_roots_f)
counts_raw_rna_roots_f <- counts_raw_rna_roots_f[featureSelectProp(counts_raw_rna_roots_f, ccc_r, 0.00005),]
dim(counts_raw_rna_roots_f)

#Filter needles
dim(counts_raw_rna_needl)
counts_raw_rna_needl_f <- counts_raw_rna_needl[featureSelect(counts_raw_rna_needl, ccc_n, 5, 2),]
dim(counts_raw_rna_needl_f)
counts_raw_rna_needl_f <- counts_raw_rna_needl_f[featureSelectProp(counts_raw_rna_needl_f, ccc_r, 0.00005),]
dim(counts_raw_rna_needl_f)

##Summarise by technical replicates
#Create new grouping variable
meta_rna_n$SampleID <- paste(meta_rna_n$sample, meta_rna_n$treatment, meta_rna_n$date, meta_rna_n$plot, sep = ".")
meta_rna_r$SampleID <- paste(meta_rna_r$sample, meta_rna_r$treatment, meta_rna_r$date, meta_rna_r$plot, sep = ".")

#Summarise the count matrices
rna_r_sum <- t(apply(counts_raw_rna_roots_f, 1, function(f){tapply(f, meta_rna_r$SampleID, mean)}))
rna_n_sum <- t(apply(counts_raw_rna_needl_f, 1, function(f){tapply(f, meta_rna_n$SampleID, mean)}))

#Summarise meta tables and adjust summed matrices¨
meta_r_sum <- unique(meta_rna_r[,-c(1,8)])
rna_r_sum <- rna_r_sum[,meta_r_sum$SampleID]
meta_n_sum <- unique(meta_rna_n[,-c(1,8)])
rna_n_sum <- rna_n_sum[,meta_n_sum$SampleID]

meta_r_sum$treatment <- factor(meta_r_sum$treatment, levels = c("Control",
                                                                "5_year",
                                                                "25_year"))

dir.create(here("RDS"))
saveRDS(meta_r_sum, here("RDS/meta_r_sum.rds"))
saveRDS(meta_n_sum, here("RDS/meta_n_sum.rds"))

#Import tpm counts
counts_tpm_rna <- read.delim(here("../data/RNA/Assembly_2012.tpm.tsv"), row.names = 1)

counts_tpm_rna <- counts_tpm_rna[rownames(counts_tpm_rna) %in% annot_tax$gene,]
#Extract roots and needles
counts_tpm_rna_roots <- data.matrix(counts_tpm_rna[,meta_rna_r$SciLifeID])
counts_tpm_rna_needl <- data.matrix(counts_tpm_rna[,meta_rna_n$SciLifeID])

#Filter roots 
dim(counts_tpm_rna_roots)
counts_tpm_rna_roots_f <- counts_tpm_rna_roots[featureSelect(counts_tpm_rna_roots, ccc_r, 5, 2),]
dim(counts_tpm_rna_roots_f)
counts_tpm_rna_roots_f <- counts_tpm_rna_roots_f[featureSelectProp(counts_tpm_rna_roots_f, ccc_r, 0.00005),]
dim(counts_tpm_rna_roots_f)

#Filter needles
dim(counts_tpm_rna_needl)
counts_tpm_rna_needl_f <- counts_tpm_rna_needl[featureSelect(counts_tpm_rna_needl, ccc_n, 5, 2),]
dim(counts_tpm_rna_needl_f)
counts_tpm_rna_needl_f <- counts_tpm_rna_needl_f[featureSelectProp(counts_tpm_rna_needl_f, ccc_r, 0.00005),]
dim(counts_tpm_rna_needl_f)

##Summarise by technical replicates
#Summarise the count matrices
rna_r_sum_tpm <- t(apply(counts_tpm_rna_roots_f, 1, function(f){tapply(f, meta_rna_r$SampleID, mean)}))
rna_n_sum_tpm <- t(apply(counts_tpm_rna_needl_f, 1, function(f){tapply(f, meta_rna_n$SampleID, mean)}))

#adjust summed matrices by meta tables
rna_r_sum_tpm <- rna_r_sum_tpm[,meta_r_sum$SampleID]
rna_n_sum_tpm <- rna_n_sum_tpm[,meta_n_sum$SampleID]

saveRDS(rna_r_sum_tpm, here("RDS/rna_r_sum_tpm.rds"))
saveRDS(rna_n_sum_tpm, here("RDS/rna_n_sum_tpm.rds"))

##Extract control samples only
control <- grepl("Control", meta_n_sum$treatment)

rna_n_sum_c <- rna_n_sum[,control]
rna_n_sum_c <- rna_n_sum_c[rowSums(rna_n_sum_c)>0,]
rna_r_sum_c <- rna_r_sum[,control]
rna_r_sum_c <- rna_r_sum_c[rowSums(rna_r_sum_c)>0,]

saveRDS(rna_n_sum_c, here("RDS/rna_n_sum_c.rds"))
saveRDS(rna_r_sum_c, here("RDS/rna_r_sum_c.rds"))

################
#Merge and VST transform (Fig 2d)
rna_sum_all <- merge(rna_n_sum_c, rna_r_sum_c, by = "row.names", all = T)
rownames(rna_sum_all) <- rna_sum_all$Row.names
rna_sum_all <- data.matrix(rna_sum_all[,-1])
rna_sum_all[is.na(rna_sum_all)] <- 0

control <- grepl("Control", meta_n_sum$treatment)

meta_sum_all_c <- rbind(meta_n_sum[control,], meta_r_sum[control,-9])

#sanity check 
all(meta_sum_all_c$SampleID == colnames(rna_sum_all))
rownames(meta_sum_all_c) <- meta_sum_all_c$SampleID

#Import to DESeq for vst 
dds_sum_all <- DESeqDataSetFromMatrix(round(rna_sum_all), colData=meta_sum_all_c, design= ~Group)

vsd_sum_all <- assay(varianceStabilizingTransformation(dds_sum_all))

#Save for Mantel/Procrustes
saveRDS(vsd_sum_all, here("RDS/vsd_sum_all.rds"))
saveRDS(meta_sum_all_c, here("RDS/meta_sum_all_c.rds"))

#################################################
###FIG 4
annot_tax_filt <- annot_tax[annot_tax$gene%in%rownames(counts_raw_rna_roots_f),]
saveRDS(annot_tax_filt, here("RDS/annot_tax_filt.rds"))

saveRDS(rna_r_sum, here("RDS/rna_r_sum.rds"))

#Proceed
dds_sum_r <- DESeqDataSetFromMatrix(round(rna_r_sum), colData = meta_r_sum, design = ~Group)

vsd_sum_r <- assay(varianceStabilizingTransformation(dds_sum_r))

##Print for Mantel/Procrustes
saveRDS(vsd_sum_r, here("RDS/vsd_sum_r.rds"))

##Process ITS data
meta_all <- read.csv(here("../doc/Amplicon/meta_NEW.csv"), row.names = 1)

#Clean a bit
meta_all$Treatment <- factor(meta_all$Treatment, levels = c("ctrl",
                                                            "5NO",
                                                            "25NO"))

meta_all$Date <- factor(meta_all$Date, levels = c("5th.June",
                                                  "24th.June",
                                                  "6th.August",
                                                  "9th.October"))

meta_all$Group <- as.factor(meta_all$Group)

#Import count matrix
count_mat_raw <- t(readRDS(here("../results/swarm/seqtab_final.rds")))
###HERE, FIX THE FILTERING!
taxonomy_its_raw <- readRDS(here("../results/taxa.rds"))
taxonomy_its_raw <- gsub("^[[:alpha:]]__", "", taxonomy_its_raw)

ccc <- meta_all$Group
names(ccc) <- colnames(count_mat_raw)

count_mat <- count_mat_raw
rownames(count_mat) <- rownames(taxonomy_its_raw)
count_mat[count_mat<10] <- 0
count_mat <- count_mat[featureSelect(count_mat, ccc, 5, 2),]
count_mat <- count_mat[featureSelectProp(count_mat, ccc, 0.00005),]

taxonomy_its <- taxonomy_its_raw[rownames(count_mat),]

count_mat <- count_mat[,rownames(meta_all)]

write.csv(count_mat, here("../data/Amplicon/count_mat_NEW.csv"), quote = FALSE)
write.csv(taxonomy_its, here("../data/Amplicon/taxa_full_NEW.csv"), quote = FALSE)

#Booleans to extract roots
roots <- grepl("root", meta_all$SampleID)

#Extract and filter root samples
count_mat_r <- count_mat[,roots]
count_mat_r <- count_mat_r[rowSums(count_mat_r)>0,]

saveRDS(count_mat_r, here("RDS/count_mat_r.rds"))
saveRDS(meta_all[roots,], here("RDS/meta_r_its.rds"))

####CLEAN UP ITS TAXONOMY
###Now proceed with ITS data
#Clean and transform to nicer format
taxonomy_its2 <- as.data.frame(cbind(SOTU=paste0("SOTU", 1:nrow(taxonomy_its)), taxonomy_its))

#Time for some clean-up of the tax table (to make the unknown annotations the same as for the RNA-Seq data)
taxonomy_its2$Family <- ifelse(grepl("Incertae", taxonomy_its2$Family), paste0("Uncertain.", taxonomy_its2$Genus), taxonomy_its2$Family)
taxonomy_its2$Family <- gsub("Uncertain.NA", NA, taxonomy_its2$Family)

taxonomy_its2$Order <- ifelse(grepl("Incertae", taxonomy_its2$Order), paste0("Uncertain.", taxonomy_its2$Family), taxonomy_its2$Order)
taxonomy_its2$Order <- gsub("Uncertain.NA", NA, taxonomy_its2$Order)
taxonomy_its2$Order <- gsub("Uncertain.Uncertain", "Uncertain", taxonomy_its2$Order) 

taxonomy_its2[is.na(taxonomy_its2)] <- "unidentified"

taxonomy_its2$Phylum <- ifelse(grepl("unidentified", taxonomy_its2$Phylum), 
                               paste0("Unclassified.", taxonomy_its2$Kingdom),
                               taxonomy_its2$Phylum)

#Let's make a function to automate this
fixNAtax <- function(tax, rank) {
  coln <- which(colnames(tax)==rank)
  namen <- colnames(tax)[coln]
  namen1 <- colnames(tax)[coln-1]
  tax[,namen] <- ifelse(grepl("Unclassified", tax[,namen1]),
                        tax[,namen1],
                        tax[,namen])
  tax[,namen] <- ifelse(grepl("unidentified", tax[,namen]),
                        paste0("Unclassified.", tax[,namen1]),
                        tax[,namen])
  return(tax)
}
taxonomy_its3 <- fixNAtax(taxonomy_its2, "Class")
taxonomy_its3 <- fixNAtax(taxonomy_its3, "Order")
taxonomy_its3 <- fixNAtax(taxonomy_its3, "Family")
taxonomy_its3 <- fixNAtax(taxonomy_its3, "Genus")
taxonomy_its3 <- fixNAtax(taxonomy_its3, "Species")
saveRDS(taxonomy_its3, here("RDS/taxonomy_cleaned_adjusted.rds"))

#########################
# corresponding supplement
saveRDS(rna_n_sum, here("RDS/rna_n_sum.rds"))
meta_n_sum$treatment <- factor(meta_n_sum$treatment, levels = c("Control",
                                                                "5_year",
                                                                "25_year"))
saveRDS(meta_n_sum, here("RDS/meta_n_sum.rds"))
#And taxonomy
annot_tax_filt_n <- annot_tax[annot_tax$gene%in%rownames(rna_n_sum),] 
saveRDS(annot_tax_filt_n, here("RDS/annot_tax_filt_n.rds"))
#Proceed
dds_sum_n <- DESeqDataSetFromMatrix(round(rna_n_sum), colData = meta_n_sum, design = ~Group)
vsd_sum_n <- assay(varianceStabilizingTransformation(dds_sum_n))
##Print for Mantel/Procrustes
saveRDS(vsd_sum_n, here("RDS/vsd_sum_n.rds"))
#ITS
#Booleans to extract needles
needl <- grepl("needl", meta_all$SampleID)

#Extract and filter needle samples
count_mat_n <- count_mat[,needl]
count_mat_n <- count_mat_n[rowSums(count_mat_n)>0,]

saveRDS(count_mat_n, here("RDS/count_mat_n.rds"))
saveRDS(meta_all[needl,], here("RDS/meta_n_its.rds"))
# 

###############FIG 7 - KEGG ortholog data
#Import KO counts
ko_all <- read.delim(here("../data/RNA/kingdom.Fungi.kos.raw.tsv"))

##Separate KO IDs and annotations
ko_annot <- ko_all[,c(1,2)]
ko_annot$KO_name <- gsub("K\\d\\d\\d\\d\\d ", "", ko_annot$KO_name)

ko_all <- ko_all[,-c(1,2)]
rownames(ko_all) <- ko_annot$ko

#Extract roots and needles
kos_raw_rna_roots <- data.matrix(ko_all[,meta_rna_r$SciLifeID])
kos_raw_rna_needl <- data.matrix(ko_all[,meta_rna_n$SciLifeID])

#Create factors for filtering
ccc_r <- as.factor(meta_rna_r$Group)
names(ccc_r) <- meta_rna_r$SciLifeID
ccc_n <- as.factor(meta_rna_n$Group)
names(ccc_n) <- meta_rna_n$SciLifeID

#Filter roots 
dim(kos_raw_rna_roots)
kos_raw_rna_roots_f <- kos_raw_rna_roots[featureSelect(kos_raw_rna_roots, ccc_r, 5, 2),]
dim(kos_raw_rna_roots_f)
kos_raw_rna_roots_f <- kos_raw_rna_roots_f[featureSelectProp(kos_raw_rna_roots_f, ccc_r, 0.00005),]
dim(kos_raw_rna_roots_f)

#Filter needles
dim(kos_raw_rna_needl)
kos_raw_rna_needl_f <- kos_raw_rna_needl[featureSelect(kos_raw_rna_needl, ccc_n, 5, 2),]
dim(kos_raw_rna_needl_f)
kos_raw_rna_needl_f <- kos_raw_rna_needl_f[featureSelectProp(kos_raw_rna_needl_f, ccc_r, 0.00005),]
dim(kos_raw_rna_needl_f)

##Summarise by technical replicates
#Sum the count matrices
kos_r_sum <- t(apply(kos_raw_rna_roots_f, 1, function(f){tapply(f, meta_rna_r$SampleID, mean)}))
dim(kos_r_sum)
kos_r_sum <- kos_r_sum[rowSums(kos_r_sum)>0,]
dim(kos_r_sum)
kos_n_sum <- t(apply(kos_raw_rna_needl_f, 1, function(f){tapply(f, meta_rna_n$SampleID, mean)}))
dim(kos_n_sum)
kos_n_sum <- kos_n_sum[rowSums(kos_n_sum)>0,]
dim(kos_n_sum)

#Sum meta tables and adjust summed matrices
kos_r_sum <- kos_r_sum[,meta_r_sum$SampleID]
kos_n_sum <- kos_n_sum[,meta_n_sum$SampleID]

#Save cleaned data
saveRDS(kos_r_sum, here("RDS/kos_r_sum.rds"))
saveRDS(kos_n_sum, here("RDS/kos_n_sum.rds"))

saveRDS(ko_annot, here("RDS/kos_annot.rds"))

dir.create(here("Figures"))



