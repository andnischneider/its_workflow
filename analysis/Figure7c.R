library(here)
library(DESeq2)
library(ggplot2)
library(reshape2)
source(here("../src/ggplot_format.R"))

##Upregulated KOs
fert25_vs_C.ALLtp.UP <- readRDS(here("RDS/fert25_vs_C.ALLtp.UP.rds"))

#Import taxonomic annotations
annot_tax_filt <- readRDS(here("RDS/annot_tax_filt.rds"))
##Import functional annotations
annot_fun <- read.delim(here("../data/RNA/annotation_results.emapper.annotations"), header = FALSE)
annot_fun$V9 <- gsub("ko:", "", annot_fun$V9)
#Select only relevant/upregulated
annot_fun_up <- annot_fun[annot_fun$V9%in%fert25_vs_C.ALLtp.UP$KO,][,c(1,9)]
#Select filtered transcripts
annot_fun_up <- annot_fun_up[annot_fun_up$V1%in%annot_tax_filt$gene,]
annot_fun_up$Family <- annot_tax_filt$family[match(annot_fun_up$V1, annot_tax_filt$gene)]

#Import meta
meta_r_sum <- readRDS(here("RDS/meta_r_sum.rds"))

#Import transcript level counts
rna_r_sum <- readRDS(here("RDS/rna_r_sum.rds"))
#Extract transcript counts associated with upregulated KOs
rna_r_sum_up <- rna_r_sum[as.character(annot_fun_up$V1),]
#Select only 25year samples
rna_r_sum_up <- rna_r_sum_up[,grepl("25_year", colnames(rna_r_sum_up))]

#How many transcripts have both KO and assignment at Family level
annot_fun_rel <- annot_fun[annot_fun$V1%in%annot_tax_filt$gene,]
annot_fun_rel <- annot_fun_rel[grepl("K", annot_fun_rel$V9),]

nrow(annot_fun_rel)/nrow(annot_tax_filt)

#Rename columns
colnames(annot_fun_up)[1:2] <- c("gene", "KO")


#Melt the subset matrix
melt_filt <- function (mat) {
  melt1 <- melt(mat, varnames = c("gene", "SampleID"), value.name = "Count")
  melt2 <- melt1[melt1$Count>0,]
  return(melt2)
}
rna_r_sum_up_m <- melt_filt(rna_r_sum_up)

rna_r_sum_up_m2 <- merge.data.frame(rna_r_sum_up_m, annot_fun_up, by = "gene")
rna_r_sum_up_m3 <- merge.data.frame(rna_r_sum_up_m2, meta_r_sum, by = "SampleID")

#Filter out top12 families
rna_r_sum_up_fam <- acast(rna_r_sum_up_m3, Family~SampleID, value.var = "Count", fun.aggregate = sum)
#Turn into percentages -> percentage of counts per family belonging to DA KOs 
rna_r_sum_up_fam <- prop.table(data.matrix(rna_r_sum_up_fam), margin = 2)*100
#Extract top 12
rna_r_sum_up_fam <- rna_r_sum_up_fam[order(rowMedians(rna_r_sum_up_fam), decreasing = TRUE)[1:12],]
rna_r_sum_up_fam_m <- melt_filt(rna_r_sum_up_fam)
rna_r_sum_up_fam_m2 <- aggregate(rna_r_sum_up_fam_m$Count, list(rna_r_sum_up_fam_m$gene), mean)
colnames(rna_r_sum_up_fam_m2) <- c("Family", "Percentage")
rna_r_sum_up_fam_m2$Spec <- "25years_UP"
rna_r_sum_up_fam_m2$Family <- factor(rna_r_sum_up_fam_m2$Family, levels = c("Unclassified.Basidiomycota",
                                                                            "Unclassified.Ascomycota",
                                                                            "Atheliaceae",
                                                                            "Gloniaceae",
                                                                            "Hyaloscyphaceae",
                                                                            "Hygrophoraceae",
                                                                            "Cortinariaceae",
                                                                            "Unclassified.Mucoromycota",
                                                                            "Mycenaceae",
                                                                            "Thelephoraceae",
                                                                            "Russulaceae",
                                                                            "Pyronemataceae"))


ggplot(rna_r_sum_up_fam_m2, aes(x = Spec, y = Percentage, fill = Family))+
  geom_bar(stat = "identity")+
  coord_cartesian(ylim = c(0,100))+
  scale_fill_manual(values = c("#bf55b9",
                               "#5cb851",
                               "#6b67cd",
                               "#a3a53e",
                               "#a283c7",
                               "#589863",
                               "#d43e6e",
                               "#49a5cf",
                               "#ca5436",
                               "#9d4b6c",
                               "#c08642",
                               "#df8093"))+
  ggformat

ggsave(here("Figures/Fig7c_UP.pdf"), width = 5, height = 6)

############################################
########DOWN
##Downregulated KOs
fert25_vs_C.ALLtp.DOWN <- readRDS(here("RDS/fert25_vs_C.ALLtp.DOWN.rds"))

#Select only relevant/downregulated
annot_fun_down <- annot_fun[annot_fun$V9%in%rownames(fert25_vs_C.ALLtp.DOWN),][,c(1,9)]
#Select filtered transcripts
annot_fun_down <- annot_fun_down[annot_fun_down$V1%in%annot_tax_filt$gene,]
annot_fun_down$Family <- annot_tax_filt$family[match(annot_fun_down$V1, annot_tax_filt$gene)]

#Turn into percentages¨
rna_r_sum_down <- rna_r_sum[as.character(annot_fun_down$V1),]
#Remove 5 year samples
rna_r_sum_down <- rna_r_sum_down[,grepl("Control", colnames(rna_r_sum_down))]

#Rename columns
colnames(annot_fun_down)[1:2] <- c("gene", "KO")

#Melt the subset matrix
rna_r_sum_down_m <- melt_filt(rna_r_sum_down)

rna_r_sum_down_m2 <- merge.data.frame(rna_r_sum_down_m, annot_fun_down, by = "gene")
rna_r_sum_down_m3 <- merge.data.frame(rna_r_sum_down_m2, meta_r_sum, by = "SampleID")

#Filter out top12 families
rna_r_sum_down_fam <- acast(rna_r_sum_down_m3, Family~SampleID, value.var = "Count", fun.aggregate = sum)
#Turn into percentages -> percentage of counts per family belonging to DA KOs 
rna_r_sum_down_fam <- prop.table(data.matrix(rna_r_sum_down_fam), margin = 2)*100
#Extract top 12
rna_r_sum_down_fam <- rna_r_sum_down_fam[order(rowMedians(rna_r_sum_down_fam), decreasing = TRUE)[1:12],]
rna_r_sum_down_fam_m <- melt_filt(rna_r_sum_down_fam)
rna_r_sum_down_fam_m2 <- aggregate(rna_r_sum_down_fam_m$Count, list(rna_r_sum_down_fam_m$gene), mean)
colnames(rna_r_sum_down_fam_m2) <- c("Family", "Percentage")
rna_r_sum_down_fam_m2$Spec <- "Control_DOWN"
rna_r_sum_down_fam_m2$Family <- factor(rna_r_sum_down_fam_m2$Family, levels = c("Cortinariaceae",
                                                                                "Unclassified.Basidiomycota",
                                                                                "Atheliaceae",
                                                                                "Unclassified.Ascomycota",
                                                                                "Hyaloscyphaceae",
                                                                                "Hygrophoraceae",
                                                                                "Gloniaceae",
                                                                                "Unclassified.Mucoromycota",
                                                                                "Mycenaceae",
                                                                                "Unclassified.Leucogyrophana",
                                                                                "Strophariaceae",
                                                                                "Sphaerobolaceae"))

ggplot(rna_r_sum_down_fam_m2, aes(x = Spec, y = Percentage, fill = Family))+
  geom_bar(stat = "identity")+
  coord_cartesian(ylim = c(0,100))+
  scale_fill_manual(values = c("#bf55b9",
                               "#5cb851",
                               "#6b67cd",
                               "#a3a53e",
                               "#a283c7",
                               "#589863",
                               "#d43e6e",
                               "#49a5cf",
                               "#ca5436",
                               "#9d4b6c",
                               "#c08642",
                               "#df8093"))+
  ggformat

ggsave(here("Figures/Fig7c_DOWN.pdf"), width = 5, height = 6)

# pdf(here("15colors.pdf"))
# scales::show_col(c("#a4b349",
#                    "#7561d0",
#                    "#5bb84d",
#                    "#c158bb",
#                    "#4a864a",
#                    "#d44270",
#                    "#52bea3",
#                    "#d14d33",
#                    "#539dd4",
#                    "#d59c41",
#                    "#7d7bc5",
#                    "#7e722c",
#                    "#dc85b0",
#                    "#c27250",
#                    "#9e4b6d"))
# dev.off()





