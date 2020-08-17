library(here)
library(ggplot2)
#library(DESeq2)
library(vegan)
library(scales)
source(here("../src/featureSelection.R"))
source(here("../src/ggplot_format.R"))

cols_new <- c(Control=rgb(0.30196078431372547, 0.6862745098039216, 0.2901960784313726), 
              "5_year"=rgb(0.21568627450980393, 0.49411764705882355, 0.7215686274509804),
              "25_year"=rgb(0.8941176470588236, 0.10196078431372549, 0.10980392156862745))

# #Import meta files and add filtering factor
# meta_rna_r <- read.csv2(here("../doc/RNA/Roots.csv"), stringsAsFactors = FALSE)
# meta_rna_r$Group <- paste(meta_rna_r$treatment, meta_rna_r$date, sep = ".")
# 
# #Import raw counts
# counts_raw_rna <- read.delim(here("../data/RNA/Assembly_2012.raw.tsv"), row.names = 1)
# #Import taxonomic annotations and select only the ones identified as fungi
# annot_tax <- read.delim(here("../data/RNA/gene_taxonomy.tsv"))
# annot_tax <- annot_tax[annot_tax$kingdom=="Fungi",]
# counts_raw_rna <- counts_raw_rna[rownames(counts_raw_rna) %in% annot_tax$gene,]
# #Extract roots
# counts_raw_rna_roots <- data.matrix(counts_raw_rna[,meta_rna_r$SciLifeID])
# 
# #Create factor for filtering
# ccc_r <- as.factor(meta_rna_r$Group)
# names(ccc_r) <- meta_rna_r$SciLifeID
# 
# #Filter 
# dim(counts_raw_rna_roots)
# counts_raw_rna_roots_f <- counts_raw_rna_roots[featureSelect(counts_raw_rna_roots, ccc_r, 5, 2),]
# dim(counts_raw_rna_roots_f)
# counts_raw_rna_roots_f <- counts_raw_rna_roots_f[featureSelectProp(counts_raw_rna_roots_f, ccc_r, 0.00005),]
# dim(counts_raw_rna_roots_f)
# 
# #Adjust tax annotations
# # annot_tax_filt <- annot_tax[annot_tax$gene%in%rownames(counts_raw_rna_roots_f),]
# # saveRDS(annot_tax_filt, here("RDS/annot_tax_filt.rds"))
# 
# ##Summarise by technical replicates
# #Create new grouping variable
# meta_rna_r$SampleID <- paste(meta_rna_r$sample, meta_rna_r$treatment, meta_rna_r$date, meta_rna_r$plot, sep = ".")
# 
# #Sum the count matrix
# rna_r_sum <- t(apply(counts_raw_rna_roots_f, 1, function(f){tapply(f, meta_rna_r$SampleID, mean)}))
# 
# #Sum meta table and adjust summed matrix
# meta_r_sum <- unique(meta_rna_r[,-c(1,8)])
# rna_r_sum <- rna_r_sum[,meta_r_sum$SampleID]
# 
# ###PCA
# #Before I run this, I have to VST transform
# meta_r_sum$treatment <- factor(meta_r_sum$treatment, levels = c("Control",
#                                                                 "5_year",
#                                                                 "25_year"))
# #Print summed matrix and meta
# saveRDS(meta_r_sum, here("RDS/meta_r_sum.rds"))
# saveRDS(rna_r_sum, here("RDS/rna_r_sum.rds"))

#Import summed matrix and meta
meta_r_sum <- readRDS(here("../Prepare_first/RDS/meta_r_sum.rds"))
rna_r_sum <- readRDS(here("../Prepare_first/RDS/rna_r_sum.rds"))

# #Proceed
# dds_sum_r <- DESeqDataSetFromMatrix(round(rna_r_sum), colData = meta_r_sum, design = ~Group)
# 
# vsd_sum_r <- assay(varianceStabilizingTransformation(dds_sum_r))
# 
# ##Print for Mantel/Procrustes
# saveRDS(vsd_sum_r, here("RDS/vsd_sum_r.rds"))

vsd_sum_r <- readRDS(here("../Prepare_first/RDS/vsd_sum_r.rds"))

#Permanova test for significance, on treatment
#And date level
adonis(dist(t(vsd_sum_r)) ~ meta_r_sum$treatment)
adonis(dist(t(vsd_sum_r)) ~ meta_r_sum$date)

#Run PCA and transform into ggplot2-friendly object
pca_sum_r <- prcomp(t(vsd_sum_r))
percent_sum_r <- round(summary(pca_sum_r)$importance[2,]*100, digits = 1)
comps_sum_r <- as.data.frame(pca_sum_r$x[, 1:3])
comps_sum_r <- comps_sum_r[meta_r_sum$SampleID,]
comps_sum_r <- cbind(comps_sum_r, meta_r_sum)

#now PCA plot
ggplot(comps_sum_r, aes(x = PC1, y = PC2))+
  geom_point(size = 4, aes(col = treatment))+
  xlab(paste0("PC1 [", percent_sum_r[1], "%]"))+
  ylab(paste0("PC2 [", percent_sum_r[2], "%]"))+
  scale_color_manual(values = cols_new)+
  ggformat_pca+
  theme(legend.position = "none")

ggsave(here("Fig4b_PCA_roots.pdf"), width = 3.5, height = 3) 

###NEW TREATMENT COLOURS
pdf(here("TreatmentColsNew.pdf"))
show_col(cols_new)
dev.off()




