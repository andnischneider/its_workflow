library(here)
library(ggplot2)
library(vegan)
source(here("../src/featureSelection.R"))
source(here("../src/ggplot_format.R"))

cols_new <- c("#f77189",
              "#97a431",
              "#36ada4",
              "#a48cf4")

meta_n_sum <- readRDS(here("RDS/meta_n_sum.rds"))
rna_n_sum <- readRDS(here("RDS/rna_n_sum.rds"))

vsd_sum_n <- readRDS(here("RDS/vsd_sum_n.rds"))

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
ggsave(here("Figures/FigS4b_PCA_needles.pdf"), width = 3.5, height = 3) 


