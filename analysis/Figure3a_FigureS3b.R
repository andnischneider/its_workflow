library(here)
library(stringr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
source(here("../src/UPSCb-common/src/R/featureSelection.R"))
source(here("../src/ggplot_format.R"))

####IMPORT RNA DATA

#Import meta files and add filtering factor
meta_rna_r <- read.csv2(here("../doc/RNA/Roots.csv"), stringsAsFactors = FALSE)
meta_rna_r$Group <- paste(meta_rna_r$treatment, meta_rna_r$date, sep = ".")

meta_rna_n <- read.csv2(here("../doc/RNA/Needles.csv"), stringsAsFactors = FALSE)
meta_rna_n$Group <- paste(meta_rna_n$treatment, meta_rna_n$date, sep = ".")

#Import raw counts
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

#Sum the count matrices
rna_r_sum <- t(apply(counts_raw_rna_roots_f, 1, function(f){tapply(f, meta_rna_r$SampleID, mean)}))
rna_n_sum <- t(apply(counts_raw_rna_needl_f, 1, function(f){tapply(f, meta_rna_n$SampleID, mean)}))

#Sum meta tables and adjust summed matrices¨
meta_r_sum <- unique(meta_rna_r[,-c(1,8)])
rna_r_sum <- rna_r_sum[,meta_r_sum$SampleID]
meta_n_sum <- unique(meta_rna_n[,-c(1,8)])
rna_n_sum <- rna_n_sum[,meta_n_sum$SampleID]

#Merge roots and needles
rna_sum_all <- merge(rna_n_sum, rna_r_sum, by = "row.names", all = T)
rownames(rna_sum_all) <- rna_sum_all$Row.names
rna_sum_all <- data.matrix(rna_sum_all[,-1])
rna_sum_all[is.na(rna_sum_all)] <- 0
colnames(rna_sum_all) <- gsub("Roots", "root", colnames(rna_sum_all))
colnames(rna_sum_all) <- gsub("Needles", "needle", colnames(rna_sum_all))

#Filter the taxonomy annotations accordingly
annot_tax_filt <- annot_tax[annot_tax$gene %in% rownames(rna_sum_all),]
colnames(annot_tax_filt) <- str_to_title(colnames(annot_tax_filt))

tax_filt <- unique(annot_tax_filt[,-1])
write.csv(tax_filt, here("../data/RNA/RNA_Unique_Taxa.csv"), quote = FALSE)

#Import ITS data
its_mat <- read.csv(here("../data/Amplicon/count_mat_NEW.csv"), row.names = 1)
tax_its <- read.csv(here("../data/Amplicon/taxa_full_NEW.csv"), row.names = 1)
tax_its$Species <- paste(tax_its$Genus, tax_its$Species)
tax_its$Species <- ifelse(grepl("NA", tax_its$Species), NA, tax_its$Species)
tax_its_un <- unique(tax_its)
write.csv(tax_its_un, here("../data/Amplicon/Amplicon_Unique_Taxa.csv"), quote = FALSE)

##Extract sets of taxa that are non-overlapping at each taxonomic level

#How about a function ?
extract_unique <- function (its, rna, level) {
  its2 <- its[!grepl("Unclassified|Uncertain|unidentified", its[,level]),]
  rna2 <- rna[!grepl("Unclassified|Uncertain|unidentified", rna[,level]),]
  its2 <- its2[!is.na(its2[,level]),]
  res_list <- list()
  res_list[["ITS"]] <- setdiff(unique(its2[,level]), unique(rna2[,level]))
  res_list[["RNA"]] <- setdiff(unique(rna2[,level]), unique(its2[,level]))
  return(res_list)
}

its_phylum <- extract_unique(tax_its_un, tax_filt, "Phylum")[[1]]
rna_phylum <- extract_unique(tax_its_un, tax_filt, "Phylum")[[2]]

its_class <- extract_unique(tax_its_un, tax_filt, "Class")[[1]]
rna_class <- extract_unique(tax_its_un, tax_filt, "Class")[[2]]

its_order <- extract_unique(tax_its_un, tax_filt, "Order")[[1]]
rna_order <- extract_unique(tax_its_un, tax_filt, "Order")[[2]]

its_family <- extract_unique(tax_its_un, tax_filt, "Family")[[1]]
rna_family <- extract_unique(tax_its_un, tax_filt, "Family")[[2]]

its_genus <- extract_unique(tax_its_un, tax_filt, "Genus")[[1]]
rna_genus <- extract_unique(tax_its_un, tax_filt, "Genus")[[2]]

its_species <- extract_unique(tax_its_un, tax_filt, "Species")[[1]]
rna_species <- extract_unique(tax_its_un, tax_filt, "Species")[[2]]

#Before calculating and merging the percentages, I have to make the sample names the same
colnames(rna_sum_all) <- gsub("Control", "ctrl", gsub("_year", "NO", gsub("R", "r", colnames(rna_sum_all))))
colnames(rna_sum_all) <- gsub("Early_June", "t1", gsub("Late_June", "t2", gsub("August", "t3", gsub("October", "t4", colnames(rna_sum_all)))))
colnames(rna_sum_all) <- gsub("13B", "13A", gsub("roots", "root", colnames(rna_sum_all)))
colnames(its_mat) <- gsub("its.", "", colnames(its_mat))

rna_sum_all <- rna_sum_all[,colnames(its_mat)]

######Extract proportions of common/unique and unknown SOTUs/Transcripts per tax level
props_its <- data.frame(cbind(Level=c(colnames(tax_its)[2:7]),
                              Total=c(rep(nrow(tax_its), 6)),
                              Unique=c(nrow(tax_its[tax_its$Phylum%in%its_phylum,]),
                                       nrow(tax_its[tax_its$Class%in%its_class,]),
                                       nrow(tax_its[tax_its$Order%in%its_order,]),
                                       nrow(tax_its[tax_its$Family%in%its_family,]),
                                       nrow(tax_its[tax_its$Genus%in%its_genus,]),
                                       nrow(tax_its[tax_its$Species%in%its_species,])),
                              Unknown=c(nrow(tax_its[is.na(tax_its$Phylum),]),
                                        nrow(tax_its[is.na(tax_its$Class),]),
                                        nrow(tax_its[is.na(tax_its$Order),]),
                                        nrow(tax_its[is.na(tax_its$Family),]),
                                        nrow(tax_its[is.na(tax_its$Genus),]),
                                        nrow(tax_its[is.na(tax_its$Species),]))), stringsAsFactors = FALSE)
props_its$Total <- as.numeric(props_its$Total)
props_its$Unique <- as.numeric(props_its$Unique)
props_its$Unknown <- as.numeric(props_its$Unknown)
props_its$Common <- props_its$Total - (props_its$Unique + props_its$Unknown)
props_its2 <- melt(props_its[,-2], value.name = "Proportion", variable.name = "Type")
props_its2$Proportion <- props_its2$Proportion / nrow(tax_its)
props_its2$Level <- factor(props_its2$Level, levels = c("Phylum",
                                                        "Class",
                                                        "Order",
                                                        "Family",
                                                        "Genus",
                                                        "Species"))
props_its2$Type <- factor(props_its2$Type, levels = c("Unknown",
                                                      "Unique",
                                                      "Common"))

ggplot(props_its2, aes(x = Level, y = Proportion, fill = Type))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("grey", brewer.pal(8, "Dark2")[c(2,6)]))+
  ggformat+
  theme(legend.position = "none")+
  ggtitle("Percentage of SOTUs belonging to \n common/unique/unknown taxonomic units")
ggsave(here("Figures/Fig3_Identified_bars_ITS_SOTULevel.pdf"), width = 5, height = 3)

#The same for the RNA data
props_rna <- data.frame(cbind(Level=c(colnames(annot_tax_filt)[4:9]),
                              Total=c(rep(nrow(annot_tax_filt), 6)),
                              Unique=c(nrow(annot_tax_filt[annot_tax_filt$Phylum%in%rna_phylum,]),
                                       nrow(annot_tax_filt[annot_tax_filt$Class%in%rna_class,]),
                                       nrow(annot_tax_filt[annot_tax_filt$Order%in%rna_order,]),
                                       nrow(annot_tax_filt[annot_tax_filt$Family%in%rna_family,]),
                                       nrow(annot_tax_filt[annot_tax_filt$Genus%in%rna_genus,]),
                                       nrow(annot_tax_filt[annot_tax_filt$Species%in%rna_species,])),
                              Unknown=c(nrow(annot_tax_filt[grepl("Unclassified", annot_tax_filt$Phylum),]),
                                        nrow(annot_tax_filt[grepl("Unclassified", annot_tax_filt$Class),]),
                                        nrow(annot_tax_filt[grepl("Unclassified", annot_tax_filt$Order),]),
                                        nrow(annot_tax_filt[grepl("Unclassified", annot_tax_filt$Family),]),
                                        nrow(annot_tax_filt[grepl("Unclassified", annot_tax_filt$Genus),]),
                                        nrow(annot_tax_filt[grepl("Unclassified", annot_tax_filt$Species),]))), stringsAsFactors = FALSE)
props_rna$Total <- as.numeric(props_rna$Total)
props_rna$Unique <- as.numeric(props_rna$Unique)
props_rna$Unknown <- as.numeric(props_rna$Unknown)
props_rna$Common <- props_rna$Total - (props_rna$Unique + props_rna$Unknown)
props_rna$Level <- str_to_title(props_rna$Level)

props_rna2 <- melt(props_rna[,-2], value.name = "Proportion", variable.name = "Type")
props_rna2$Proportion <- props_rna2$Proportion / nrow(annot_tax_filt)

props_rna2$Level <- factor(props_rna2$Level, levels = c("Phylum",
                                                        "Class",
                                                        "Order",
                                                        "Family",
                                                        "Genus",
                                                        "Species"))
props_rna2$Type <- factor(props_rna2$Type, levels = c("Unknown",
                                                      "Unique",
                                                      "Common"))

ggplot(props_rna2, aes(x = Level, y = Proportion, fill = Type))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("grey", brewer.pal(8, "Dark2")[5:6]))+
  ggformat+
  theme(legend.position = "none")+
  ggtitle("Percentage of Transcripts belonging to \n common/unique/unknown taxonomic units")
ggsave(here("Figures/Fig3_Identified_bars_RNA_TranscriptLevel.pdf"), width = 5, height = 3)

##############################
##########And now get the percentages of reads.
tax_its <- cbind(SOTU=paste0("SOTU", 1:nrow(tax_its)), tax_its)
rownames(its_mat) <- tax_its$SOTU
percentages <- cbind(RNA_Phylum=colSums(rna_sum_all[rownames(rna_sum_all)%in%annot_tax_filt$Gene[annot_tax_filt$Phylum%in%rna_phylum],])/colSums(rna_sum_all),
                     ITS_Phylum=colSums(its_mat[rownames(its_mat)%in%tax_its$SOTU[tax_its$Phylum%in%its_phylum],])/colSums(its_mat),
                     RNA_Class=colSums(rna_sum_all[rownames(rna_sum_all)%in%annot_tax_filt$Gene[annot_tax_filt$Class%in%rna_class],])/colSums(rna_sum_all),
                     ITS_Class=colSums(its_mat[rownames(its_mat)%in%tax_its$SOTU[tax_its$Class%in%its_class],])/colSums(its_mat),
                     RNA_Order=colSums(rna_sum_all[rownames(rna_sum_all)%in%annot_tax_filt$Gene[annot_tax_filt$Order%in%rna_order],])/colSums(rna_sum_all),
                     ITS_Order=colSums(its_mat[rownames(its_mat)%in%tax_its$SOTU[tax_its$Order%in%its_order],])/colSums(its_mat),
                     RNA_Family=colSums(rna_sum_all[rownames(rna_sum_all)%in%annot_tax_filt$Gene[annot_tax_filt$Family%in%rna_family],])/colSums(rna_sum_all),
                     ITS_Family=colSums(its_mat[rownames(its_mat)%in%tax_its$SOTU[tax_its$Family%in%its_family],])/colSums(its_mat),
                     RNA_Genus=colSums(rna_sum_all[rownames(rna_sum_all)%in%annot_tax_filt$Gene[annot_tax_filt$Genus%in%rna_genus],])/colSums(rna_sum_all),
                     ITS_Genus=colSums(its_mat[rownames(its_mat)%in%tax_its$SOTU[tax_its$Genus%in%its_genus],])/colSums(its_mat),
                     RNA_Species=colSums(rna_sum_all[rownames(rna_sum_all)%in%annot_tax_filt$Gene[annot_tax_filt$Species%in%rna_species],])/colSums(rna_sum_all),
                     ITS_Species=colSums(its_mat[rownames(its_mat)%in%tax_its$SOTU[tax_its$Species%in%its_species],])/colSums(its_mat))

##Second mat with all unknown percentages
unknowns <- cbind(RNA_Phylum=colSums(rna_sum_all[rownames(rna_sum_all)%in%annot_tax_filt$Gene[grepl("Unclassified", annot_tax_filt$Phylum)],])/colSums(rna_sum_all),
                  ITS_Phylum=colSums(its_mat[rownames(its_mat)%in%tax_its$SOTU[is.na(tax_its$Phylum)],])/colSums(its_mat),
                  RNA_Class=colSums(rna_sum_all[rownames(rna_sum_all)%in%annot_tax_filt$Gene[grepl("Unclassified", annot_tax_filt$Class)],])/colSums(rna_sum_all),
                  ITS_Class=colSums(its_mat[rownames(its_mat)%in%tax_its$SOTU[is.na(tax_its$Class)],])/colSums(its_mat),
                  RNA_Order=colSums(rna_sum_all[rownames(rna_sum_all)%in%annot_tax_filt$Gene[grepl("Unclassified", annot_tax_filt$Order)],])/colSums(rna_sum_all),
                  ITS_Order=colSums(its_mat[rownames(its_mat)%in%tax_its$SOTU[is.na(tax_its$Order)],])/colSums(its_mat),
                  RNA_Family=colSums(rna_sum_all[rownames(rna_sum_all)%in%annot_tax_filt$Gene[grepl("Unclassified", annot_tax_filt$Family)],])/colSums(rna_sum_all),
                  ITS_Family=colSums(its_mat[rownames(its_mat)%in%tax_its$SOTU[is.na(tax_its$Family)],])/colSums(its_mat),
                  RNA_Genus=colSums(rna_sum_all[rownames(rna_sum_all)%in%annot_tax_filt$Gene[grepl("Unclassified", annot_tax_filt$Genus)],])/colSums(rna_sum_all),
                  ITS_Genus=colSums(its_mat[rownames(its_mat)%in%tax_its$SOTU[is.na(tax_its$Genus)],])/colSums(its_mat),
                  RNA_Species=colSums(rna_sum_all[rownames(rna_sum_all)%in%annot_tax_filt$Gene[grepl("Unclassified", annot_tax_filt$Species)],])/colSums(rna_sum_all),
                  ITS_Species=colSums(its_mat[rownames(its_mat)%in%tax_its$SOTU[is.na(tax_its$Species)],])/colSums(its_mat))

unknowns2 <- melt(unknowns)

percentages2 <- melt(percentages)

percentages2 <- cbind(percentages2, str_split_fixed(percentages2$Var1, "\\.", 4))
percentages2 <- cbind(percentages2, str_split_fixed(percentages2$Var2, "_", 2))[,-2]
colnames(percentages2) <- c("SampleID", "Unique", "Type", "Treatment", "Timepoint", "Plot", "SeqType", "TaxLevel")
percentages2$Treatment <- factor(percentages2$Treatment, levels = c("ctrl",
                                                                    "5NO",
                                                                    "25NO"))

percentages2$TaxLevel <- factor(percentages2$TaxLevel, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))

#Unite unknown and known percentages
percentages2$Unknown <- unknowns2$value

percentages3 <- percentages2
percentages3$Type_Seq_Tax <- paste(percentages3$Type, percentages3$SeqType, percentages3$TaxLevel, sep = "_")

percentages4 <- aggregate(percentages3$Unique, list(percentages3$Type_Seq_Tax), mean)
colnames(percentages4) <- c("Group", "Non_Overlap")

uk4 <- aggregate(percentages3$Unknown, list(percentages3$Type_Seq_Tax), mean)

percentages4$Unknown <- uk4$x

percentages4$Overlap <- 1 - (percentages4$Non_Overlap + percentages4$Unknown)
percentages4 <- melt(percentages4)

percentages4_rna <- percentages4[grepl("RNA", percentages4$Group),]
percentages4_rna$Level <- gsub("^[[:alpha:]]+_RNA_", "", percentages4_rna$Group)
percentages4_rna$Level <- factor(percentages4_rna$Level, levels = c("Phylum",
                                                                    "Class",
                                                                    "Order",
                                                                    "Family",
                                                                    "Genus",
                                                                    "Species"))
percentages4_rna$Type <- gsub("^([[:alpha:]]+)_RNA_.*", "\\1", percentages4_rna$Group)
percentages4_rna$variable <- factor(percentages4_rna$variable, levels = c("Unknown",
                                                                          "Non_Overlap",
                                                                          "Overlap"))

ggplot(percentages4_rna, aes(x = Level, y = value))+
  geom_bar(stat = "identity", aes(fill = variable))+
  scale_fill_manual(values = c("grey", brewer.pal(8, "Dark2")[5:6]))+
  facet_wrap(~ Type)+
  ggformat+
  theme(legend.position = "none")#+
  #ggtitle("Percentage of Transcript read counts \n belonging to common/unique/unknown \n taxonomic units")

ggsave(here("Figures/FigS3_Identified_bars_Perc_RNA.pdf"), width = 10, height = 3)

percentages4_its <- percentages4[grepl("ITS", percentages4$Group),]
percentages4_its$Level <- gsub("^[[:alpha:]]+_ITS_", "", percentages4_its$Group)
percentages4_its$Level <- factor(percentages4_its$Level, levels = c("Phylum",
                                                                    "Class",
                                                                    "Order",
                                                                    "Family",
                                                                    "Genus",
                                                                    "Species"))
percentages4_its$Type <- gsub("^([[:alpha:]]+)_ITS_.*", "\\1", percentages4_its$Group)
percentages4_its$variable <- factor(percentages4_its$variable, levels = c("Unknown",
                                                                          "Non_Overlap",
                                                                          "Overlap"))

ggplot(percentages4_its, aes(x = Level, y = value))+
  geom_bar(stat = "identity", aes(fill = variable))+
  scale_fill_manual(values = c("grey", brewer.pal(8, "Dark2")[c(2,6)]))+
  ggformat+
  facet_wrap(~ Type)+
  theme(legend.position = "none")#+
  #ggtitle("Percentage of SOTU read counts \n belonging to common/unique/unknown \n taxonomic units")

ggsave(here("Figures/FigS3_Identified_bars_Perc_ITS.pdf"), width = 10, height = 3)








