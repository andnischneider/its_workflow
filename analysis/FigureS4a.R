library(here)
library(RColorBrewer)
library(ggplot2)
library(vegan)
library(phyloseq)
source(here("../src/ggplot_format.R"))

# #Import metadata
# meta_all <- read.csv(here("../doc/Amplicon/meta_NEW.csv"), row.names = 1)
# 
# #Clean a bit
# meta_all$Treatment <- factor(meta_all$Treatment, levels = c("ctrl",
#                                                             "5NO",
#                                                             "25NO"))
# 
# meta_all$Date <- factor(meta_all$Date, levels = c("5th.June",
#                                                   "24th.June",
#                                                   "6th.August",
#                                                   "9th.October"))
# 
# #Import count matrix
# count_mat <- read.csv(here("../data/Amplicon/count_mat_NEW.csv"), row.names = 1)
# 
# #Booleans to extract needles
# needl <- grepl("needl", meta_all$SampleID)
# 
# #Extract and filter needle samples
# count_mat_n <- count_mat[,needl]
# count_mat_n <- count_mat_n[rowSums(count_mat_n)>0,]
# 
# saveRDS(count_mat_n, here("RDS/count_mat_n.rds"))
# saveRDS(meta_all[needl,], here("RDS/meta_n_its.rds"))

count_mat_n <- readRDS(here("../Prepare_first/RDS/count_mat_n.rds"))
meta_n_its <- readRDS(here("../Prepare_first/RDS/meta_n_its.rds"))

#Rarefy for Bray-Curtis PCoA
count_mat_n_rar <- t(rrarefy(t(count_mat_n), min(colSums(count_mat_n))))

#Permanova test for significance, date and treatment
adonis(vegdist(t(count_mat_n_rar)) ~ meta_n_its[,"Treatment"])
adonis(vegdist(t(count_mat_n_rar)) ~ meta_n_its[,"Timepoint"])

#Ordination and Plot
ps_its_n_rar <- phyloseq(otu_table(count_mat_n_rar, taxa_are_rows = TRUE),
                         sample_data(meta_n_its))

n.ord <- ordinate(ps_its_n_rar, "MDS", "bray")
plot_ordination(ps_its_n_rar, n.ord, type = "samples", color = "Date")+
  geom_point(size = 4)+
  xlab(paste0("PCoA 1 [", round(n.ord$values$Relative_eig[1]*100, digits = 1), "%]"))+
  ylab(paste0("PCoA 2 [", round(n.ord$values$Relative_eig[2]*100, digits = 1), "%]"))+
  scale_color_manual(values = cols_new_date)+
  #scale_y_reverse()+
  ggformat_pca+
  theme(legend.position = "none")
ggsave("FigS4a_PCoA_needles.pdf", width = 3.5, height = 3)

########Mantel & Procrustes
#Import vst transformed RNA data
count_mat_n_rna <- readRDS(here("../Prepare_first/RDS/vsd_sum_n.rds"))

#adjust RNA names
colnames(count_mat_n_rna) <- gsub("Needles", "its.needle", colnames(count_mat_n_rna))
colnames(count_mat_n_rna) <- gsub("_year", "NO", colnames(count_mat_n_rna))
colnames(count_mat_n_rna) <- gsub("Control", "ctrl", colnames(count_mat_n_rna))
colnames(count_mat_n_rna) <- gsub("Early_June", "t1", gsub("Late_June", "t2", gsub("August", "t3", gsub("October", "t4", colnames(count_mat_n_rna)))))
colnames(count_mat_n_rna) <- gsub("13B", "13A", colnames(count_mat_n_rna))

#Same order
count_mat_n_rna <- count_mat_n_rna[,colnames(count_mat_n_rar)]

#Make ordination object
pca_rna <- prcomp(t(count_mat_n_rna))

#For its data (NMDS)
#nmds_its <- metaMDS(t(count_mat_r_rar))

##PCoA
pcoa_its <- cmdscale(vegdist(t(count_mat_n_rar)))

#Procrustes test
proc1 <- protest(pcoa_its, pca_rna, symmetric = TRUE)
proc1

#Mantel test
mantel(vegdist(t(count_mat_n_rar)), dist(t(count_mat_n_rna)))




