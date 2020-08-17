library(here)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
source(here("../src/ggplot_format.R"))

vsd_sum_all <- readRDS(here("RDS/vsd_sum_all.rds"))
meta_sum_all_c <- readRDS(here("RDS/meta_sum_all_c.rds"))

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
ggsave(here("Figures/Fig2d_PCA.pdf"), width = 6, height = 3) 
