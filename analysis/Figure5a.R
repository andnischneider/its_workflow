library(here)
library(stringr)
library(reshape2)
library(matrixStats)
library(ggplot2)
library(scales)
source(here("../src/ggplot_format.R"))

##RNA first
#Import summed count matrix
rna_r_sum <- readRDS(here("RDS/rna_r_sum.rds"))
#Taxonomy
rna_tax_r <- readRDS(here("RDS/annot_tax_filt.rds"))
#Metadata
rna_meta_r_sum <- readRDS(here("RDS/meta_r_sum.rds"))
rownames(rna_meta_r_sum) <- rna_meta_r_sum$SampleID

#Turn count mat into proportions for visualization
rna_r_sum_p <- prop.table(rna_r_sum, margin = 2)

#Melt count matrix to enable summing by tax levels
melt_filt2 <- function (mat) {
  melt1 <- melt(mat, varnames = c("gene", "SampleID"), value.name = "Count")
  melt2 <- melt1[melt1$Count>0,]
  return(melt2)
}
rna_r_melt <- melt_filt2(rna_r_sum_p)
rna_r_melt2 <- merge.data.frame(rna_r_melt, rna_meta_r_sum, by = "SampleID")
rna_r_melt2 <- merge.data.frame(rna_r_melt2, rna_tax_r, by = "gene")

#Aggregate the percentage by family
rna_r_fam <- aggregate(rna_r_melt2$Count, list(paste(rna_r_melt2$SampleID, rna_r_melt2$family, sep = "-")), sum)
colnames(rna_r_fam) <- c("SampleID_Family", "Count")

#Split annotations
rna_r_fam[,c(3,4)] <- str_split_fixed(rna_r_fam$SampleID_Family, "-", 2)
rna_r_fam <- rna_r_fam[,-1]
  
rna_r_fam_mat <- acast(rna_r_fam, V4~V3, value.var = "Count")
rna_r_fam_mat[is.na(rna_r_fam_mat)] <- 0
#Extract top 12 families for plotting
rna_fam_top12 <- rownames(rna_r_fam_mat)[order(rowMedians(rna_r_fam_mat), decreasing = TRUE)][1:12]

#Summarise family percentages for area plot
colnames(rna_r_fam)[2:3] <- c("SampleID", "family")
rna_r_fam2 <- merge.data.frame(rna_r_fam, rna_meta_r_sum, by = "SampleID")
rna_r_fam2_top12 <- rna_r_fam2[rna_r_fam2$family%in%rna_fam_top12,]
rna_r_fam2_top12$date <- factor(rna_r_fam2_top12$date, levels = c("Early_June", "Late_June", "August", "October"))

#Take mean per sample and family for plotting
rna_r_fam2_top12_mean <- aggregate(rna_r_fam2_top12$Count, list(paste(rna_r_fam2_top12$Group, rna_r_fam2_top12$family)), mean)
#Format for ggplot
rna_r_fam2_top12_mean[,c(3,4)] <- str_split_fixed(rna_r_fam2_top12_mean$Group.1, "\\s", 2)
rna_r_fam2_top12_mean[,c(5,6)] <- str_split_fixed(rna_r_fam2_top12_mean$V3, "\\.", 2)
rna_r_fam2_top12_mean <- rna_r_fam2_top12_mean[,c(2,4,5,6)]
colnames(rna_r_fam2_top12_mean) <- c("Proportion", "Family", "Condition", "Date")
rna_r_fam2_top12_mean$Condition <- factor(rna_r_fam2_top12_mean$Condition, levels = c("Control", "5_year", "25_year"))
rna_r_fam2_top12_mean$Date <- factor(rna_r_fam2_top12_mean$Date, levels = c("Early_June", "Late_June", "August", "October"))

rna_r_fam2_top12_mean$Family <- factor(rna_r_fam2_top12_mean$Family, levels = c("Gloniaceae", "Hyaloscyphaceae", "Mycenaceae", "Unclassified.Ascomycota", "Unclassified.Basidiomycota", "Unclassified.Paracoccidioides", "Atheliaceae", "Cortinariaceae", "Hygrophoraceae", "Russulaceae", "Strophariaceae", "Venturiaceae"))

ggplot(rna_r_fam2_top12_mean, aes(x = Date, y = Proportion))+
  geom_area(aes(group = Family, fill = Family))+
  facet_wrap(~Condition)+
  scale_fill_manual(values = c("#bc74a7",
                               "#5bbf5f",
                               "#b859c4",
                               "#a8b746",
                               "#6e67cb",
                               "#5b8627",
                               "#cc478a",
                               "#4cbfb2",
                               "#d24359",
                               "#4f8d58",
                               "#cc572c",
                               "#6a8fcd",
                               "#d89c44",
                               "#c46f63",
                               "#8e7435"))+
  ggformat+
  theme(axis.text.x = element_text(angle = 90))

ggsave(here("Figures/Fig5a_RNA_Area.pdf"), width = 8, height = 8)

#######################################################
###Now the same procedure with the ITS data
its_r_sum <- readRDS(here("RDS/count_mat_r.rds"))
#Taxonomy
its_tax_r <- readRDS("RDS/taxonomy_cleaned_adjusted.rds")
its_tax_r <- its_tax_r[rownames(its_tax_r)%in%rownames(its_r_sum),]
rownames(its_r_sum) <- its_tax_r$SOTU
#Metadata
its_meta_r_sum <- readRDS(here("RDS/meta_r_its.rds"))
rownames(its_meta_r_sum) <- its_meta_r_sum$SampleID
its_meta_r_sum$Group <- paste(its_meta_r_sum$Timepoint, its_meta_r_sum$Treatment, sep = "_")

#Turn count mat into proportions for visualization
its_r_sum_p <- prop.table(data.matrix(its_r_sum), margin = 2)

#Now we need to create a summed matrix on Genus level, through
#phyloseq and metagenomeseq is too slow with this big dataset.
melt_filt3 <- function (mat) {
  melt1 <- melt(mat, varnames = c("SOTU", "SampleID"), value.name = "Count")
  melt2 <- melt1[melt1$Count>0,]
  return(melt2)
}
its_r_melt <- melt_filt3(its_r_sum_p)
its_r_melt2 <- merge.data.frame(its_r_melt, its_meta_r_sum, by = "SampleID")
its_r_melt2 <- merge.data.frame(its_r_melt2, its_tax_r, by = "SOTU")

#Aggregate the percentage by family
its_r_fam <- aggregate(its_r_melt2$Count, list(paste(its_r_melt2$SampleID, its_r_melt2$Family, sep = "-")), sum)
colnames(its_r_fam) <- c("SampleID_Family", "Count")

#Split annotations
its_r_fam[,c(3,4)] <- str_split_fixed(its_r_fam$SampleID_Family, "-", 2)
its_r_fam <- its_r_fam[,-1]

its_r_fam_mat <- acast(its_r_fam, V4~V3, value.var = "Count")
its_r_fam_mat[is.na(its_r_fam_mat)] <- 0
#Extract top 12 families for plotting
its_fam_top12 <- rownames(its_r_fam_mat)[order(rowMedians(its_r_fam_mat), decreasing = TRUE)][1:12]

#Summarise family percentages for area plot
colnames(its_r_fam)[2:3] <- c("SampleID", "Family")
its_r_fam2 <- merge.data.frame(its_r_fam, its_meta_r_sum, by = "SampleID")
its_r_fam2_top12 <- its_r_fam2[its_r_fam2$Family%in%its_fam_top12,]

#Take mean per sample and family for plotting
its_r_fam2_top12_mean <- aggregate(its_r_fam2_top12$Count, list(paste(its_r_fam2_top12$Group, its_r_fam2_top12$Family)), mean)
#Format for ggplot
its_r_fam2_top12_mean[,c(3,4)] <- str_split_fixed(its_r_fam2_top12_mean$Group.1, "\\s", 2)
its_r_fam2_top12_mean[,c(5,6)] <- str_split_fixed(its_r_fam2_top12_mean$V3, "_", 2)
its_r_fam2_top12_mean <- its_r_fam2_top12_mean[,c(2,4,5,6)]
colnames(its_r_fam2_top12_mean) <- c("Proportion", "Family", "Timepoint", "Treatment")
its_r_fam2_top12_mean$Treatment <- factor(its_r_fam2_top12_mean$Treatment, levels = c("ctrl", "5NO", "25NO"))
levels(its_r_fam2_top12_mean$Treatment) <- c("Control", "5_year", "25_year")
its_r_fam2_top12_mean$Timepoint <- factor(its_r_fam2_top12_mean$Timepoint, levels = c("t1", "t2", "t3", "t4"))
levels(its_r_fam2_top12_mean$Timepoint) <- c("Early_June", "Late_June", "August", "October")

its_r_fam2_top12_mean$Family <- factor(its_r_fam2_top12_mean$Family, levels = c("Archaeorhizomycetaceae", "Helotiaceae", "Myxotrichaceae", "Thelephoraceae", "Unclassified.Helotiales", "Vibrisseaceae", "Atheliaceae", "Cortinariaceae", "Hygrophoraceae", "Russulaceae", "Strophariaceae", "Venturiaceae"))

ggplot(its_r_fam2_top12_mean, aes(x = Timepoint, y = Proportion))+
  geom_area(aes(group = Family, fill = Family))+
  facet_wrap(~Treatment)+
  scale_fill_manual(values = c("#bc74a7",
                               "#5bbf5f",
                               "#b859c4",
                               "#a8b746",
                               "#6e67cb",
                               "#5b8627",
                               "#cc478a",
                               "#4cbfb2",
                               "#d24359",
                               "#4f8d58",
                               "#cc572c",
                               "#6a8fcd",
                               "#d89c44",
                               "#c46f63",
                               "#8e7435"))+
  ggformat+
  theme(axis.text.x = element_text(angle = 90))

ggsave(here("Figures/Fig5a_ITS_Area.pdf"), width = 8, height = 8)

# ##In total we need 18 colors
# colors18 <- c("#618cce",
#               "#65b344",
#               "#a757c8",
#               "#b8b343",
#               "#606dda",
#               "#dd9435",
#               "#7d5a9e",
#               "#5ebf8a",
#               "#d14698",
#               "#3a814f",
#               "#d44458",
#               "#4bbdd1",
#               "#cd5530",
#               "#cd8cd8",
#               "#797e35",
#               "#9d4660",
#               "#b97a4b",
#               "#df839a")
# 
# pdf(here("colors18.pdf"))
# show_col(colors18)
# dev.off()



#########################
# Extract genera to see which ones to put in bold in the heatmap
genera_its <- unique(its_tax_r$Genus)
genera_rna <- as.character(unique(rna_tax_r$genus))
intersect(genera_its, genera_rna)

