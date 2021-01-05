library(here)
library(pheatmap)
library(matrixStats)
library(ggplot2)
library(grid)
library(VennDiagram)
source(here("../src/ggplot_format.R"))
source(here("../src/UPSCb-common/src/R/gopher.R"))
source(here("../src/Rtoolbox/src/plotEnrichedTreemap.R"))

#Import counts
mat_r_sum <- readRDS(here("../Prepare_first/RDS/rna_r_sum.rds"))

tax <- readRDS(here("../Prepare_first/RDS/annot_tax_filt.rds"))
#How many Cortinarius transcripts
nrow(tax[tax$genus=="Cortinarius",])
tax_cort <- tax[tax$genus=="Cortinarius",]
rownames(tax_cort) <- tax_cort$gene

#subset counts to Cortinarius
mat_r_sum_cort <- mat_r_sum[rownames(tax_cort),]

meta_r_sum <- readRDS(here("../Prepare_first/RDS/meta_r_sum.rds"))
#Add factor to meta table
rownames(meta_r_sum) <- meta_r_sum$SampleID
meta_r_sum$Group2 <- as.factor(paste(meta_r_sum$treatment, meta_r_sum$date, sep="_"))
meta_r_sum$Group2 <- factor(meta_r_sum$Group2, levels = c("Control_Early_June","Control_Late_June","Control_August","Control_October",
                                                          "5_year_Early_June",  "5_year_Late_June", "5_year_August", "5_year_October",
                                                          "25_year_Early_June","25_year_Late_June", "25_year_August", "25_year_October"))



#Import functional annotations
annot_fun_all <- read.delim(here("../data/RNA/annotation_results.emapper.annotations"), header = FALSE)
annot_fun <- annot_fun_all[annot_fun_all$V1%in%tax$gene,]
annot_fun$V9 <- gsub("ko:", "", annot_fun$V9)



#Import vst counts
mat_sum_vst <- readRDS("../Prepare_first/RDS/vsd_sum_r.rds")
mat_sum_vst_cort <- mat_sum_vst[rownames(tax_cort),]
#Remove 5 year samples
mat_sum_vst_cort2 <- mat_sum_vst_cort[,!grepl("Roots.5_", colnames(mat_sum_vst_cort))]

#Colors for heatmap
sample_cols <- data.frame(Date=meta_r_sum$date, Treatment=meta_r_sum$treatment)
sample_cols$Date <- factor(sample_cols$Date, levels = c("Early_June", "Late_June", "August", "October"))
rownames(sample_cols) <- colnames(mat_sum_vst_cort)
sample_cols <- sample_cols[colnames(mat_sum_vst_cort2),]
ann_colors <- list(
  Treatment = c(Control=cols_new_treat[1], "25_year"=cols_new_treat[3]),
  Date = c(Early_June=cols_new_date[1], Late_June=cols_new_date[2], August=cols_new_date[3], October=cols_new_date[4])
)

#With row scaling
out2 <- pheatmap(mat_sum_vst_cort2[-which(rowSds(mat_sum_vst_cort2)==0),],
                 color = colorRampPalette(c("#FFFFFF", "#FFFAF7", "#FFF4EE", "#FCBBA1", "#F6553C", "#C8171C", "#6C010E"))(50),
                 show_rownames = F, 
                 show_colnames = F,
                 annotation_col = sample_cols,
                 annotation_colors = ann_colors,
                 clustering_method = "ward.D",
                 scale = "row") 
plot(out2$tree_row)
abline(h=2000, col = "red", lty = 2, lwd = 2)
gen_clus2 <- as.data.frame(cutree(out2$tree_row, h=2000))
colnames(gen_clus2) <- "Cluster"
gen_clus2$Cluster <- as.factor(gen_clus2$Cluster)

ann_colors2 <- list(
  Treatment = c(Control=cols_new_treat[1], "25_year"=cols_new_treat[3]),
  Date = c(Early_June=cols_new_date[1], Late_June=cols_new_date[2], August=cols_new_date[3], October=cols_new_date[4]),
  Cluster = c("1"="#214E61", "2"="#077BAD", "3"="#4EB8E6")
)

pheatmap(mat_sum_vst_cort2[-which(rowSds(mat_sum_vst_cort2)==0),],
         color = colorRampPalette(c("#FFFFFF", "#FFFAF7", "#FFF4EE", "#FCBBA1", "#F6553C", "#C8171C", "#6C010E"))(50),
         show_rownames = F, 
         show_colnames = F,
         annotation_row = gen_clus2,
         annotation_col = sample_cols,
         annotation_colors = ann_colors2,
         clustering_method = "ward.D",
         scale = "row") 
#Looks good, print.
pheatmap(mat_sum_vst_cort2[-which(rowSds(mat_sum_vst_cort2)==0),],
         color = colorRampPalette(c("#FFFFFF", "#FFFAF7", "#FFF4EE", "#FCBBA1", "#F6553C", "#C8171C", "#6C010E"))(50),
         show_rownames = F, 
         show_colnames = F,
         annotation_row = gen_clus2,
         annotation_col = sample_cols,
         annotation_colors = ann_colors2,
         clustering_method = "ward.D",
         scale = "row",
         filename = "Fig_8a_hm.pdf")

#Now the corresponding KEGG orthologs
kos_cluster <- function (cluster, gen_clu) {
  cluster_genes <- rownames(gen_clu)[gen_clu[["Cluster"]]==cluster]
  #Select corresponding KOs
  kos <- annot_fun$V9[annot_fun$V1%in%cluster_genes]
  #process to remove empty elements 
  kos <- kos[kos!=""]
  #process to separate the multiple hits
  kos <- unique(strsplit(paste(kos, collapse = ","), ",")[[1]])
  #kos <- strsplit(paste(kos, collapse = ","), ",")[[1]]
  return(kos)
}
cluster1_kos <- kos_cluster("1", gen_clus2)
cluster2_kos <- kos_cluster("2", gen_clus2)
cluster3_kos <- kos_cluster("3", gen_clus2)

cluster1_genes <- rownames(gen_clus2)[gen_clus2$Cluster=="1"]
cluster2_genes <- rownames(gen_clus2)[gen_clus2$Cluster=="2"]
cluster3_genes <- rownames(gen_clus2)[gen_clus2$Cluster=="3"]

pdf("Fig8a_Venn.pdf")
grid.newpage()
grid.draw(venn.diagram(list(cluster1=cluster1_kos, cluster2=cluster2_kos, cluster3=cluster3_kos), NULL))
dev.off()

##GO enrichment per cluster
#All Cortinarius GOs
gos_cort <- annot_fun$V7[annot_fun$V1%in%tax_cort$gene]
gos_cort <- gos_cort[gos_cort!=""]
gos_cort1 <- unique(strsplit(paste(gos_cort, collapse = ","), ",")[[1]])
gos_cort2 <- strsplit(paste(gos_cort, collapse = ","), ",")[[1]]

#All fungal GOs
gos_all <- annot_fun$V7
gos_all <- gos_all[gos_all!=""]
gos_all1 <- unique(strsplit(paste(gos_all, collapse = ","), ",")[[1]])
gos_all2 <- strsplit(paste(gos_all, collapse = ","), ",")[[1]]


enr.c1_go <- gopher(cluster1_genes, background = rownames(mat_r_sum_cort), task = list("go"), alpha=0.001, url = "fungi2012_cortinarius", endpoint = "enrichment")
enr.c2_go <- gopher(cluster2_genes, background = rownames(mat_r_sum_cort), task = list("go"), alpha=0.001, url = "fungi2012_cortinarius", endpoint = "enrichment")
enr.c3_go <- gopher(cluster3_genes, background = rownames(mat_r_sum_cort), task = list("go"), alpha=0.001, url = "fungi2012_cortinarius", endpoint = "enrichment")

pdf(here("Fig8b_TreemapC1_GO.pdf"), width = 4, height = 2.5)
plotEnrichedTreemap(enr.c1_go, enrichment = "go", namespace = "BP", sizeCol = "nt", clusterColor = "#214E61")
dev.off()


pdf(here("Fig8b_TreemapC3_GO.pdf"), width = 4, height = 2.5)
plotEnrichedTreemap(enr.c3_go, enrichment = "go", namespace = "BP", sizeCol = "nt", clusterColor = "#4EB8E6")
dev.off()