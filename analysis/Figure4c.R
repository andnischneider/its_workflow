library(here)
library(reshape2)
library(vegan)
library(nlme)
library(multcomp)
library(ggplot2)
source(here("../src/ggplot_format.R"))

cols_new <- c(rgb(0.30196078431372547, 0.6862745098039216, 0.2901960784313726), 
              rgb(0.21568627450980393, 0.49411764705882355, 0.7215686274509804),
              rgb(0.8941176470588236, 0.10196078431372549, 0.10980392156862745))


#RNA DATA FIRST
#Shannon diversity on Genus level
annot_tax_filt_ps <- readRDS(here("RDS/annot_tax_filt.rds"))
rownames(annot_tax_filt_ps) <- annot_tax_filt_ps$gene
meta_r_sum <- readRDS(here("RDS/meta_r_sum.rds"))
rownames(meta_r_sum) <- meta_r_sum$SampleID
#Now we need to create a summed matrix on Genus level, through
#phyloseq and metagenomeseq is too slow with this big dataset.
melt_filt2 <- function (mat) {
  melt1 <- melt(mat, varnames = c("gene", "SampleID"), value.name = "Count")
  melt2 <- melt1[melt1$Count>0,]
  return(melt2)
}
rna_r_sum <- readRDS(here("RDS/rna_r_sum.rds"))
counts_sum_melted_rna <- melt_filt2(rna_r_sum)
counts_sum_melted_rna2 <- merge.data.frame(counts_sum_melted_rna, meta_r_sum, by = "SampleID")
counts_sum_melted_rna2 <- merge.data.frame(counts_sum_melted_rna2, annot_tax_filt_ps, by = "gene")

counts_sum_genus_rna <- acast(counts_sum_melted_rna2, genus~SampleID, value.var = "Count", fun.aggregate = sum)
shannon_root_genus_rna <- cbind(meta_r_sum, Shannon=diversity(counts_sum_genus_rna, MARGIN = 2))

taxonomy_its3 <- readRDS(here("RDS/taxonomy_cleaned_adjusted.rds"))

##Continue with meta and count data
meta_r_its <- readRDS(here("RDS/meta_r_its.rds"))
count_mat_r <- readRDS(here("RDS/count_mat_r.rds"))

taxonomy_its3_r <- taxonomy_its3[rownames(count_mat_r),]

#Adjust SOTU names in matrix for easier summarizing
count_mat_r2 <- count_mat_r
rownames(count_mat_r2) <- taxonomy_its3_r$SOTU

counts_sum_melted_its <- melt_filt2(data.matrix(count_mat_r2))
colnames(counts_sum_melted_its)[1] <- "SOTU"
counts_sum_melted_its2 <- merge.data.frame(counts_sum_melted_its, meta_r_its, by = "SampleID")
counts_sum_melted_its2 <- merge.data.frame(counts_sum_melted_its2, taxonomy_its3_r, by = "SOTU")

###Now proceed with ITS data
counts_sum_genus_its <- acast(counts_sum_melted_its2, Genus~SampleID, value.var = "Count", fun.aggregate = sum)
shannon_root_genus_its <- cbind(meta_r_its, Shannon=diversity(counts_sum_genus_its, MARGIN = 2))
shannon_root_genus_its$SampleID <- as.character(shannon_root_genus_its$SampleID)

#Adjust RNA annotations before merging
shannon_root_genus_rna$SampleID <- gsub("Roots", "its.root", shannon_root_genus_rna$SampleID)
shannon_root_genus_rna$SampleID <- gsub("_year", "NO", shannon_root_genus_rna$SampleID)
shannon_root_genus_rna$SampleID <- gsub("Control", "ctrl", shannon_root_genus_rna$SampleID)
shannon_root_genus_rna$SampleID <- gsub("Early_June", "t1", gsub("Late_June", "t2", gsub("August", "t3", gsub("October", "t4", shannon_root_genus_rna$SampleID))))
shannon_root_genus_rna$SampleID <- gsub("13B", "13A", shannon_root_genus_rna$SampleID)
rownames(shannon_root_genus_rna) <- shannon_root_genus_rna$SampleID
shannon_root_genus_rna <- shannon_root_genus_rna[shannon_root_genus_its$SampleID,]
#Merge
root_shannon_genus_both <- cbind(shannon_root_genus_rna[,c(4:7)], Shannon_ITS=shannon_root_genus_its$Shannon, Shannon_RNA=shannon_root_genus_rna$Shannon)

#Plot
ggplot(root_shannon_genus_both, aes(x = Shannon_ITS, y = Shannon_RNA, col = treatment))+
  geom_point(size = 3)+
  scale_color_manual(values = cols_new)+
  ggformat+
  coord_cartesian(xlim = c(0, 3),
                  ylim = c(0, 3))+
  theme(axis.text.x = element_text(angle = 0),
        axis.title.x = element_text(),
        legend.position = "none")

ggsave(here("Figures/Fig4c_Shannon_scatter.pdf"), width = 2.2, height = 2.4)
#Pearson corr
cor(root_shannon_genus_both$Shannon_ITS, root_shannon_genus_both$Shannon_RNA)

#Test for significance; RNA treatment
shannon_model_rna <- lme(Shannon_RNA~treatment, random = ~1|plot, data = root_shannon_genus_both, method = "ML")
anova(shannon_model_rna)
#Not really
#shannon_posthoc <- glht(shannon_model_rna, linfct=mcp(treatment="Tukey"))
#summary(shannon_posthoc)

#Test for significance; ITS treatment
shannon_model_its <- lme(Shannon_ITS~treatment, random = ~1|plot, data = root_shannon_genus_both, method = "ML")
anova(shannon_model_its)
#Yes
shannon_posthoc <- glht(shannon_model_its, linfct=mcp(treatment="Tukey"))
summary(shannon_posthoc)


#Test for significance; RNA timepoint
root_shannon_genus_both$date <- factor(root_shannon_genus_both$date, levels = c("Early_June", "Late_June", "August", "October"))
shannon_model_rna_d <- lme(Shannon_RNA~date, random = ~1|plot, data = root_shannon_genus_both, method = "ML")
anova(shannon_model_rna_d)
#Nope
#shannon_posthoc_d <- glht(shannon_model_rna_d, linfct=mcp(date="Tukey"))
#summary(shannon_posthoc_d)

#Test for significance; ITS timepoint
shannon_model_its_d <- lme(Shannon_ITS~date, random = ~1|plot, data = root_shannon_genus_both, method = "ML")
anova(shannon_model_its_d)
#Nope
#shannon_posthoc_its_d <- glht(shannon_model_its_d, linfct=mcp(date="Tukey"))
#summary(shannon_posthoc_its_d)



