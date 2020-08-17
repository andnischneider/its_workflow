library(here)
library(ggplot2)
library(DESeq2)
library(vegan)
source(here("../src/featureSelection.R"))
source(here("../src/ggplot_format.R"))

cols_new <- c("#f77189",
              "#97a431",
              "#36ada4",
              "#a48cf4")

# #Import meta files and add filtering factor
# meta_rna_n <- read.csv2(here("../doc/RNA/Needles.csv"), stringsAsFactors = FALSE)
# meta_rna_n$Group <- paste(meta_rna_n$treatment, meta_rna_n$date, sep = ".")
# 
# #Import raw counts
# counts_raw_rna <- read.delim(here("../data/RNA/Assembly_2012.raw.tsv"), row.names = 1)
# #Import taxonomic annotations and select only the ones identified as fungi
# annot_tax <- read.delim(here("../data/RNA/gene_taxonomy.tsv"))
# annot_tax <- annot_tax[annot_tax$kingdom=="Fungi",]
# counts_raw_rna <- counts_raw_rna[rownames(counts_raw_rna) %in% annot_tax$gene,]
# #Extract needles
# counts_raw_rna_needl <- data.matrix(counts_raw_rna[,meta_rna_n$SciLifeID])
# 
# #Create factor for filtering
# ccc_n <- as.factor(meta_rna_n$Group)
# names(ccc_n) <- meta_rna_n$SciLifeID
# 
# #Filter
# dim(counts_raw_rna_needl)
# counts_raw_rna_needl_f <- counts_raw_rna_needl[featureSelect(counts_raw_rna_needl, ccc_n, 5, 2),]
# dim(counts_raw_rna_needl_f)
# counts_raw_rna_needl_f <- counts_raw_rna_needl_f[featureSelectProp(counts_raw_rna_needl_f, ccc_n, 0.00005),]
# dim(counts_raw_rna_needl_f)
# 
# ##Summarise by technical replicates
# #Create new grouping variable
# meta_rna_n$SampleID <- paste(meta_rna_n$sample, meta_rna_n$treatment, meta_rna_n$date, meta_rna_n$plot, sep = ".")
# 
# #Sum the count matrix
# rna_n_sum <- t(apply(counts_raw_rna_needl_f, 1, function(f){tapply(f, meta_rna_n$SampleID, mean)}))
# 
# #Sum meta table and adjust summed matrix
# meta_n_sum <- unique(meta_rna_n[,-c(1,8)])
# rna_n_sum <- rna_n_sum[,meta_n_sum$SampleID]
# 
# ###PCA
# #Before I run this, I have to VST transform
# meta_n_sum$treatment <- factor(meta_n_sum$treatment, levels = c("Control",
#                                                                 "5_year",
#                                                                 "25_year"))
# 
# #Print summed matrix and meta
# saveRDS(meta_n_sum, here("RDS/meta_n_sum.rds"))
# saveRDS(rna_n_sum, here("RDS/rna_n_sum.rds"))
# #And taxonomy
# annot_tax_filt_n <- annot_tax[annot_tax$gene%in%rownames(rna_n_sum),] 
# saveRDS(annot_tax_filt_n, here("RDS/annot_tax_filt_n.rds"))

meta_n_sum <- readRDS(here("../Prepare_first/RDS/meta_n_sum.rds"))
rna_n_sum <- readRDS(here("../Prepare_first/RDS/rna_n_sum.rds"))

# #Proceed
# dds_sum_n <- DESeqDataSetFromMatrix(round(rna_n_sum), colData = meta_n_sum, design = ~Group)
# 
# vsd_sum_n <- assay(varianceStabilizingTransformation(dds_sum_n))
# 
# ##Print for Mantel/Procrustes
# saveRDS(vsd_sum_n, here("RDS/vsd_sum_n.rds"))

vsd_sum_n <- readRDS(here("../Prepare_first/RDS/vsd_sum_n.rds"))

#Permanova test for significance, on treatment
#And date level
adonis(dist(t(vsd_sum_n)) ~ meta_n_sum$treatment)
adonis(dist(t(vsd_sum_n)) ~ meta_n_sum$date)

#Run PCA and transform into ggplot2-friendly object
pca_sum_n <- prcomp(t(vsd_sum_n))
percent_sum_n <- round(summary(pca_sum_n)$importance[2,]*100, digits = 1)
comps_sum_n <- as.data.frame(pca_sum_n$x[,1:3])
comps_sum_n <- comps_sum_n[meta_n_sum$SampleID,]
comps_sum_n <- cbind(comps_sum_n, meta_n_sum)
comps_sum_n$date <- factor(comps_sum_n$date, levels = c("Early_June",
                                                        "Late_June",
                                                        "August",
                                                        "October"))

#Plot PCA, color by date, since that had higher significance
ggplot(comps_sum_n, aes(x = PC1, y = PC3))+
  geom_point(size = 4, aes(col = date))+
  xlab(paste0("PC1 [", percent_sum_n[1], "%]"))+
  ylab(paste0("PC2 [", percent_sum_n[2], "%]"))+
  scale_color_manual(values = cols_new)+
  ggformat_pca+
  theme(legend.position = "none")
ggsave(here("FigS4b_PCA_needles.pdf"), width = 3.5, height = 3) 


