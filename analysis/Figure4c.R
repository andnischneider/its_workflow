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
annot_tax_filt_ps <- readRDS(here("../Prepare_first/RDS/annot_tax_filt.rds"))
rownames(annot_tax_filt_ps) <- annot_tax_filt_ps$gene
meta_r_sum <- readRDS(here("../Prepare_first/RDS/meta_r_sum.rds"))
rownames(meta_r_sum) <- meta_r_sum$SampleID
#Now we need to create a summed matrix on Genus level, through
#phyloseq and metagenomeseq is too slow with this big dataset.
melt_filt2 <- function (mat) {
  melt1 <- melt(mat, varnames = c("gene", "SampleID"), value.name = "Count")
  melt2 <- melt1[melt1$Count>0,]
  return(melt2)
}
rna_r_sum <- readRDS(here("../Prepare_first/RDS/rna_r_sum.rds"))
counts_sum_melted_rna <- melt_filt2(rna_r_sum)
counts_sum_melted_rna2 <- merge.data.frame(counts_sum_melted_rna, meta_r_sum, by = "SampleID")
counts_sum_melted_rna2 <- merge.data.frame(counts_sum_melted_rna2, annot_tax_filt_ps, by = "gene")

counts_sum_genus_rna <- acast(counts_sum_melted_rna2, genus~SampleID, value.var = "Count", fun.aggregate = sum)
shannon_root_genus_rna <- cbind(meta_r_sum, Shannon=diversity(counts_sum_genus_rna, MARGIN = 2))

# ###Now proceed with ITS data
# taxonomy_its <- read.csv(here("../data/Amplicon/taxa_full_NEW.csv"), row.names = 1, stringsAsFactors = FALSE)
# #Clean and transform to nicer format
# taxonomy_its2 <- cbind(SOTU=paste0("SOTU", 1:nrow(taxonomy_its)), taxonomy_its)
# 
# #Time for some clean-up of the tax table (to make the unknown annotations the same as for the RNA-Seq data)
# taxonomy_its2$Family <- ifelse(grepl("Incertae", taxonomy_its2$Family), paste0("Uncertain.", taxonomy_its2$Genus), taxonomy_its2$Family)
# taxonomy_its2$Family <- gsub("Uncertain.NA", NA, taxonomy_its2$Family)
# 
# taxonomy_its2$Order <- ifelse(grepl("Incertae", taxonomy_its2$Order), paste0("Uncertain.", taxonomy_its2$Family), taxonomy_its2$Order)
# taxonomy_its2$Order <- gsub("Uncertain.NA", NA, taxonomy_its2$Order)
# taxonomy_its2$Order <- gsub("Uncertain.Uncertain", "Uncertain", taxonomy_its2$Order) 
# 
# taxonomy_its2[is.na(taxonomy_its2)] <- "unidentified"
# 
# taxonomy_its2$Phylum <- ifelse(grepl("unidentified", taxonomy_its2$Phylum), 
#                                paste0("Unclassified.", taxonomy_its2$Kingdom),
#                                taxonomy_its2$Phylum)
# 
# #Let's make a function to automate this
# fixNAtax <- function(tax, rank) {
#   coln <- which(colnames(tax)==rank)
#   namen <- colnames(tax)[coln]
#   namen1 <- colnames(tax)[coln-1]
#   tax[,namen] <- ifelse(grepl("Unclassified", tax[,namen1]),
#                         tax[,namen1],
#                         tax[,namen])
#   tax[,namen] <- ifelse(grepl("unidentified", tax[,namen]),
#                         paste0("Unclassified.", tax[,namen1]),
#                         tax[,namen])
#   return(tax)
# }
# taxonomy_its3 <- fixNAtax(taxonomy_its2, "Class")
# taxonomy_its3 <- fixNAtax(taxonomy_its3, "Order")
# taxonomy_its3 <- fixNAtax(taxonomy_its3, "Family")
# taxonomy_its3 <- fixNAtax(taxonomy_its3, "Genus")
# taxonomy_its3 <- fixNAtax(taxonomy_its3, "Species")
# saveRDS(taxonomy_its3, here("RDS/taxonomy_cleaned_adjusted.rds"))

taxonomy_its3 <- readRDS(here("../Prepare_first/RDS/taxonomy_cleaned_adjusted.rds"))

##Continue with meta and count data
meta_r_its <- readRDS(here("../Prepare_first/RDS/meta_r_its.rds"))
count_mat_r <- readRDS(here("../Prepare_first/RDS/count_mat_r.rds"))

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

#Try to include formula and R square
# lm_eqn = function(m) {
#   
#   l <- list(a = format(coef(m)[1], digits = 2),
#             b = format(abs(coef(m)[2]), digits = 2),
#             r2 = format(summary(m)$r.squared, digits = 3));
#   
#   if (coef(m)[2] >= 0)  {
#     eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
#   } else {
#     eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
#   }
#   
#   as.character(as.expression(eq));                 
# }

#Plot
ggplot(root_shannon_genus_both, aes(x = Shannon_ITS, y = Shannon_RNA, col = treatment))+
  geom_point(size = 3)+
  #stat_summary(fun.data=mean_cl_normal)+
  #stat_smooth(method="lm", se= FALSE, fullrange = TRUE, aes(group = 1))+
  scale_color_manual(values = cols_new)+
  ggformat+
  coord_cartesian(xlim = c(0, 3),
                  ylim = c(0, 3))+
  theme(axis.text.x = element_text(angle = 0),
        axis.title.x = element_text(),
        legend.position = "none")#+
#geom_text(aes(x = 0.7, y = 0.9, label = lm_eqn(lm(root_shannon_genus_both$Shannon_RNA~root_shannon_genus_both$Shannon_ITS))), parse = TRUE)

ggsave(here("Fig4c_Shannon_scatter.pdf"), width = 2.2, height = 2.4)
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



