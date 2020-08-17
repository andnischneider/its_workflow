library(here)
# library(DESeq2)
library(ggplot2)
library(RColorBrewer)
source(here("../src/ggplot_format.R"))

# #Import filtered and cleaned control matrices created 
# #from Figure2b.R
# 
# rna_n_sum_c <- readRDS(here("RDS/rna_n_sum_c.rds"))
# rna_r_sum_c <- readRDS(here("RDS/rna_r_sum_c.rds"))
# 
# #Combining roots and needles, for control PCA
# rna_sum_all <- merge(rna_n_sum_c, rna_r_sum_c, by = "row.names", all = T)
# rownames(rna_sum_all) <- rna_sum_all$Row.names
# rna_sum_all <- data.matrix(rna_sum_all[,-1])
# rna_sum_all[is.na(rna_sum_all)] <- 0
# 
# #Import filtered meta control tables from Figure2a.R
# meta_n_sum <- readRDS(here("RDS/meta_n_sum.rds"))
# meta_r_sum <- readRDS(here("RDS/meta_r_sum.rds"))
# 
# control <- grepl("Control", meta_n_sum$treatment)
# 
# meta_sum_all_c <- rbind(meta_n_sum[control,], meta_r_sum[control,-9])
# 
# #sanity check 
# all(meta_sum_all_c$SampleID == colnames(rna_sum_all))
# rownames(meta_sum_all_c) <- meta_sum_all_c$SampleID
# 
# #Import to DESeq for vst 
# dds_sum_all <- DESeqDataSetFromMatrix(round(rna_sum_all), colData=meta_sum_all_c, design= ~Group)
# 
# vsd_sum_all <- assay(varianceStabilizingTransformation(dds_sum_all))

vsd_sum_all <- readRDS(here("../Prepare_first/RDS/vsd_sum_all.rds"))
meta_sum_all_c <- readRDS(here("../Prepare_first/RDS/meta_sum_all_c.rds"))

##PCA of this
pca_sum_all <- prcomp(t(vsd_sum_all))
percent_sum_all <- round(summary(pca_sum_all)$importance[2,]*100, digits = 1)
comps_sum_all <- as.data.frame(pca_sum_all$x[,1:15])
comps_sum_all <- comps_sum_all[meta_sum_all_c$SampleID,]
comps_sum_all <- cbind(comps_sum_all, meta_sum_all_c)

ggplot(comps_sum_all, aes(x = PC1, y = PC2))+
  geom_point(size = 5, aes(col = sample, shape = plot))+
  xlab(paste0("PC1 [", percent_sum_all[1], "%]"))+
  ylab(paste0("PC2 [", percent_sum_all[2], "%]"))+
  scale_color_manual(values = brewer.pal(8,"Dark2")[1:2])+
  ggformat_pca
ggsave(here("Fig2d_PCA.pdf"), width = 6, height = 3) 
