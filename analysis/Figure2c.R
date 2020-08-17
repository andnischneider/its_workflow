library(here)
library(vegan)
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
source(here("../src/ggplot_format.R"))

#Import metadata
meta_all <- read.csv(here("../doc/Amplicon/meta_NEW.csv"), row.names = 1)

#Clean a bit
meta_all$Treatment <- factor(meta_all$Treatment, levels = c("ctrl",
                                                            "5NO",
                                                            "25NO"))

meta_all$Date <- factor(meta_all$Date, levels = c("5th.June",
                                                  "24th.June",
                                                  "6th.August",
                                                  "9th.October"))

#Import count matrix
count_mat <- read.csv(here("../data/Amplicon/count_mat_NEW.csv"), row.names = 1)

##### ONLY CONTROL SAMPLES ROOTS AND NEEDLES COMBINED
control <- grepl("ctrl", meta_all$SampleID)

count_mat_ctrl <- count_mat[,control]
count_mat_ctrl <- count_mat_ctrl[rowSums(count_mat_ctrl)>0,]

count_mat_ctrl_rar <- t(rrarefy(t(count_mat_ctrl), min(colSums(count_mat_ctrl))))

#Import into phyloseq
ps_ctrl <- phyloseq(otu_table(count_mat_ctrl_rar, taxa_are_rows = TRUE),
                    sample_data(meta_all[control,]))
c.ord <- ordinate(ps_ctrl, "MDS", "bray")
plot_ordination(ps_ctrl, c.ord, type = "samples", color = "SampleType", shape = "SoilPlot")+
  geom_point(size = 5)+
  xlab(paste0("PCoA 1 [", round(c.ord$values$Relative_eig[1]*100, digits = 1), "%]"))+
  ylab(paste0("PCoA 2 [", round(c.ord$values$Relative_eig[2]*100, digits = 1), "%]"))+
  scale_color_manual(values = brewer.pal(8,"Dark2")[1:2])+
  ggformat_pca

ggsave(here("Figures/Fig2c_PCoA.pdf"), width = 6, height = 3)

######Mantel & Procrustes
##Import vst transformed RNA control data (Script Figure2d.R needs to be run first)
count_mat_rna <- readRDS(here("RDS/vsd_sum_all.rds"))

#adjust RNA names
colnames(count_mat_rna) <- gsub("Roots", "its.root", colnames(count_mat_rna))
colnames(count_mat_rna) <- gsub("Needles", "its.needle", colnames(count_mat_rna))
colnames(count_mat_rna) <- gsub("_year", "NO", colnames(count_mat_rna))
colnames(count_mat_rna) <- gsub("Control", "ctrl", colnames(count_mat_rna))
colnames(count_mat_rna) <- gsub("Early_June", "t1", gsub("Late_June", "t2", gsub("August", "t3", gsub("October", "t4", colnames(count_mat_rna)))))
colnames(count_mat_rna) <- gsub("13B", "13A", colnames(count_mat_rna))

#Same order
count_mat_rna <- count_mat_rna[,colnames(count_mat_ctrl_rar)]

#Make ordination object
pca_rna <- prcomp(t(count_mat_rna))

#For its data (NMDS)
#nmds_its <- metaMDS(t(count_mat_ctrl_rar))

##PCoA
pcoa_its <- cmdscale(vegdist(t(count_mat_ctrl_rar)))

#Procrustes test
protest(pcoa_its, pca_rna, symmetric = TRUE)


#Mantel test
mantel(vegdist(t(count_mat_ctrl_rar)), dist(t(count_mat_rna)))





