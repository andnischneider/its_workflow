library(here)
library(reshape2)
library(vegan)
library(ggplot2)
library(nlme)
library(multcomp)
source(here("../src/ggplot_format.R"))

#RNA DATA FIRST
#Shannon diversity on Genus level
annot_tax_filt_n <- readRDS(here("RDS/annot_tax_filt_n.rds"))
rownames(annot_tax_filt_n) <- annot_tax_filt_n$gene
meta_n_sum <- readRDS(here("RDS/meta_n_sum.rds"))
rownames(meta_n_sum) <- meta_n_sum$SampleID
#Now we need to create a summed matrix on Genus level, through
#phyloseq and metagenomeseq is too slow with this big dataset.
melt_filt2 <- function (mat) {
  melt1 <- melt(mat, varnames = c("gene", "SampleID"), value.name = "Count")
  melt2 <- melt1[melt1$Count>0,]
  return(melt2)
}
rna_n_sum <- readRDS(here("RDS/rna_n_sum.rds"))
counts_sum_melted_rna_n <- melt_filt2(rna_n_sum)
counts_sum_melted_rna_n2 <- merge.data.frame(counts_sum_melted_rna_n, meta_n_sum, by = "SampleID", all.x = TRUE)
counts_sum_melted_rna_n2 <- merge.data.frame(counts_sum_melted_rna_n2, annot_tax_filt_n, by = "gene", all.x = TRUE)

counts_sum_genus_rna_n <- acast(counts_sum_melted_rna_n2, genus~SampleID, value.var = "Count", fun.aggregate = sum)
shannon_needl_genus_rna <- cbind(meta_n_sum, Shannon=diversity(counts_sum_genus_rna_n, MARGIN = 2))

taxonomy_its3 <- readRDS(here("RDS/taxonomy_cleaned_adjusted.rds"))

##Continue with meta and count data
meta_n_its <- readRDS(here("RDS/meta_n_its.rds"))
count_mat_n <- readRDS(here("RDS/count_mat_n.rds"))

taxonomy_its3_n <- taxonomy_its3[rownames(count_mat_n),]

#Adjust SOTU names in matrix for easier summarizing
count_mat_n2 <- count_mat_n
rownames(count_mat_n2) <- taxonomy_its3_n$SOTU

counts_sum_melted_its_n <- melt_filt2(data.matrix(count_mat_n2))
colnames(counts_sum_melted_its_n)[1] <- "SOTU"
counts_sum_melted_its_n2 <- merge.data.frame(counts_sum_melted_its_n, meta_n_its, by = "SampleID")
counts_sum_melted_its_n2 <- merge.data.frame(counts_sum_melted_its_n2, taxonomy_its3_n, by = "SOTU")

###Cast matrix and calculate shannon
counts_sum_genus_its_n <- acast(counts_sum_melted_its_n2, Genus~SampleID, value.var = "Count", fun.aggregate = sum)
shannon_needl_genus_its_n <- cbind(meta_n_its, Shannon=diversity(counts_sum_genus_its_n, MARGIN = 2))
shannon_needl_genus_its_n$SampleID <- as.character(shannon_needl_genus_its_n$SampleID)

#Adjust RNA annotations before merging
shannon_needl_genus_rna$SampleID <- gsub("Needles", "its.needle", shannon_needl_genus_rna$SampleID)
shannon_needl_genus_rna$SampleID <- gsub("_year", "NO", shannon_needl_genus_rna$SampleID)
shannon_needl_genus_rna$SampleID <- gsub("Control", "ctrl", shannon_needl_genus_rna$SampleID)
shannon_needl_genus_rna$SampleID <- gsub("Early_June", "t1", gsub("Late_June", "t2", gsub("August", "t3", gsub("October", "t4", shannon_needl_genus_rna$SampleID))))
shannon_needl_genus_rna$SampleID <- gsub("13B", "13A", shannon_needl_genus_rna$SampleID)
rownames(shannon_needl_genus_rna) <- shannon_needl_genus_rna$SampleID
shannon_needl_genus_rna <- shannon_needl_genus_rna[shannon_needl_genus_its_n$SampleID,]
#Merge
needl_shannon_genus_both <- cbind(shannon_needl_genus_rna[,c(4:7)], Shannon_ITS=shannon_needl_genus_its_n$Shannon, Shannon_RNA=shannon_needl_genus_rna$Shannon)

needl_shannon_genus_both$date <- factor(needl_shannon_genus_both$date, levels = c("Early_June",
                                                                                  "Late_June",
                                                                                  "August",
                                                                                  "October"))

#Plot
ggplot(needl_shannon_genus_both, aes(x = Shannon_ITS, y = Shannon_RNA, col = date))+
  geom_point(size = 3.2)+
  scale_color_manual(values = cols_new_date)+
  ggformat+
  coord_cartesian(xlim = c(0, 3.5),
                  ylim = c(0, 1.9))+
  theme(axis.text.x = element_text(angle = 0),
        axis.title.x = element_text(),
        legend.position = "none")

ggsave(here("Figures/FigS4c_Shannon_scatter_need.pdf"), width = 2.3, height = 2.4)

#Pearson corr
cor(needl_shannon_genus_both$Shannon_ITS, needl_shannon_genus_both$Shannon_RNA)

#Test for significance; RNA treatment
shannon_model_rna <- lme(Shannon_RNA~treatment, random = ~1|plot, data = needl_shannon_genus_both, method = "ML")
anova(shannon_model_rna)
#No sig.
#shannon_posthoc <- glht(shannon_model_rna, linfct=mcp(treatment="Tukey"))
#summary(shannon_posthoc)

#Test for significance; ITS treatment
shannon_model_its <- lme(Shannon_ITS~treatment, random = ~1|plot, data = needl_shannon_genus_both, method = "ML")
anova(shannon_model_its)
#Yes
shannon_posthoc <- glht(shannon_model_its, linfct=mcp(treatment="Tukey"))
summary(shannon_posthoc)

#Test for significance; RNA timepoint
needl_shannon_genus_both$date <- factor(needl_shannon_genus_both$date, levels = c("Early_June", "Late_June", "August", "October"))
shannon_model_rna_d <- lme(Shannon_RNA~date, random = ~1|plot, data = needl_shannon_genus_both, method = "ML")
anova(shannon_model_rna_d)
#Yes
shannon_posthoc_d <- glht(shannon_model_rna_d, linfct=mcp(date="Tukey"))
summary(shannon_posthoc_d)
#October significantly higher than the rest

#Test for significance; ITS timepoint
shannon_model_its_d <- lme(Shannon_ITS~date, random = ~1|plot, data = needl_shannon_genus_both, method = "ML")
anova(shannon_model_its_d)
#Nope
# shannon_posthoc_its_d <- glht(shannon_model_its_d, linfct=mcp(date="Tukey"))
# summary(shannon_posthoc_its_d)
#No significant differences
