library(here)
library(vegan)
library(phyloseq)
library(ggplot2)
source(here("../src/ggplot_format.R"))

cols_new <- c(rgb(0.30196078431372547, 0.6862745098039216, 0.2901960784313726), 
              rgb(0.21568627450980393, 0.49411764705882355, 0.7215686274509804),
              rgb(0.8941176470588236, 0.10196078431372549, 0.10980392156862745))



count_mat_r <- readRDS(here("RDS/count_mat_r.rds"))
meta_r_its <- readRDS(here("RDS/meta_r_its.rds"))

#Rarefy for Bray-Curtis PCoA
count_mat_r_rar <- t(rrarefy(t(count_mat_r), min(colSums(count_mat_r))))

#Permanova test for significance, date and treatment
adonis(vegdist(t(count_mat_r_rar)) ~ meta_r_its[,"Treatment"])
adonis(vegdist(t(count_mat_r_rar)) ~ meta_r_its[,"Timepoint"])

#Ordination and Plot
ps_its_r_rar <- phyloseq(otu_table(count_mat_r_rar, taxa_are_rows = TRUE),
                         sample_data(meta_r_its))

r.ord <- ordinate(ps_its_r_rar, "MDS", "bray")
plot_ordination(ps_its_r_rar, r.ord, type = "samples", color = "Treatment")+
  geom_point(size = 4)+
  xlab(paste0("PCoA 1 [", round(r.ord$values$Relative_eig[1]*100, digits = 1), "%]"))+
  ylab(paste0("PCoA 2 [", round(r.ord$values$Relative_eig[2]*100, digits = 1), "%]"))+
  scale_color_manual(values = cols_new)+
  scale_y_reverse()+
  ggformat_pca+
  theme(legend.position = "none")
ggsave("Figures/Fig4a_PCoA_roots.pdf", width = 3.5, height = 3)


########Mantel & Procrustes
#Import vst transformed RNA data
count_mat_rna <- readRDS(here("RDS/vsd_sum_r.rds"))

#adjust RNA names
colnames(count_mat_rna) <- gsub("Roots", "its.root", colnames(count_mat_rna))
colnames(count_mat_rna) <- gsub("_year", "NO", colnames(count_mat_rna))
colnames(count_mat_rna) <- gsub("Control", "ctrl", colnames(count_mat_rna))
colnames(count_mat_rna) <- gsub("Early_June", "t1", gsub("Late_June", "t2", gsub("August", "t3", gsub("October", "t4", colnames(count_mat_rna)))))
colnames(count_mat_rna) <- gsub("13B", "13A", colnames(count_mat_rna))

#Same order
count_mat_rna <- count_mat_rna[,colnames(count_mat_r_rar)]

#Make ordination object
pca_rna <- prcomp(t(count_mat_rna))

#For its data (NMDS)
#nmds_its <- metaMDS(t(count_mat_r_rar))

##PCoA
pcoa_its <- cmdscale(vegdist(t(count_mat_r_rar)))

#Procrustes test
protest(pcoa_its, pca_rna, symmetric = TRUE)

#Mantel test
mantel(vegdist(t(count_mat_r_rar)), dist(t(count_mat_rna)))


