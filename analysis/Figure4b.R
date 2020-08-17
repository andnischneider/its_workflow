library(here)
library(ggplot2)
library(vegan)
library(scales)
source(here("../src/UPSCb-common/src/R/featureSelection.R"))
source(here("../src/ggplot_format.R"))

cols_new <- c(Control=rgb(0.30196078431372547, 0.6862745098039216, 0.2901960784313726), 
              "5_year"=rgb(0.21568627450980393, 0.49411764705882355, 0.7215686274509804),
              "25_year"=rgb(0.8941176470588236, 0.10196078431372549, 0.10980392156862745))


#Import summed matrix and meta
meta_r_sum <- readRDS(here("RDS/meta_r_sum.rds"))
rna_r_sum <- readRDS(here("RDS/rna_r_sum.rds"))

vsd_sum_r <- readRDS(here("RDS/vsd_sum_r.rds"))

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

ggsave(here("Figures/Fig4b_PCA_roots.pdf"), width = 3.5, height = 3) 

# ###NEW TREATMENT COLOURS
# pdf(here("TreatmentColsNew.pdf"))
# show_col(cols_new)
# dev.off()




