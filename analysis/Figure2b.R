library(here)
library(grid)
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
source(here("../src/ggplot_format.R"))

###A

# #Import meta files and add filtering factor
# meta_rna_r <- read.csv2(here("../doc/RNA/Roots.csv"), stringsAsFactors = FALSE)
# meta_rna_r$Group <- paste(meta_rna_r$treatment, meta_rna_r$date, sep = ".")
# 
# meta_rna_n <- read.csv2(here("../doc/RNA/Needles.csv"), stringsAsFactors = FALSE)
# meta_rna_n$Group <- paste(meta_rna_n$treatment, meta_rna_n$date, sep = ".")
# 
# #Import raw counts
# counts_raw_rna <- read.delim(here("../data/RNA/Assembly_2012.raw.tsv"), row.names = 1)
# #Import taxonomic annotations and select only the ones identified as fungi
# annot_tax <- read.delim(here("../data/RNA/gene_taxonomy.tsv"))
# annot_tax <- annot_tax[annot_tax$kingdom=="Fungi",]
# counts_raw_rna <- counts_raw_rna[rownames(counts_raw_rna) %in% annot_tax$gene,]
# #Extract roots and needles
# counts_raw_rna_roots <- data.matrix(counts_raw_rna[,meta_rna_r$SciLifeID])
# counts_raw_rna_needl <- data.matrix(counts_raw_rna[,meta_rna_n$SciLifeID])
# 
# #Create factors for filtering
# ccc_r <- as.factor(meta_rna_r$Group)
# names(ccc_r) <- meta_rna_r$SciLifeID
# ccc_n <- as.factor(meta_rna_n$Group)
# names(ccc_n) <- meta_rna_n$SciLifeID
# 
# #Filter roots 
# dim(counts_raw_rna_roots)
# counts_raw_rna_roots_f <- counts_raw_rna_roots[featureSelect(counts_raw_rna_roots, ccc_r, 5, 2),]
# dim(counts_raw_rna_roots_f)
# counts_raw_rna_roots_f <- counts_raw_rna_roots_f[featureSelectProp(counts_raw_rna_roots_f, ccc_r, 0.00005),]
# dim(counts_raw_rna_roots_f)
# 
# #Filter needles
# dim(counts_raw_rna_needl)
# counts_raw_rna_needl_f <- counts_raw_rna_needl[featureSelect(counts_raw_rna_needl, ccc_n, 5, 2),]
# dim(counts_raw_rna_needl_f)
# counts_raw_rna_needl_f <- counts_raw_rna_needl_f[featureSelectProp(counts_raw_rna_needl_f, ccc_r, 0.00005),]
# dim(counts_raw_rna_needl_f)
# 
# ##Summarise by technical replicates
# #Create new grouping variable
# meta_rna_n$SampleID <- paste(meta_rna_n$sample, meta_rna_n$treatment, meta_rna_n$date, meta_rna_n$plot, sep = ".")
# meta_rna_r$SampleID <- paste(meta_rna_r$sample, meta_rna_r$treatment, meta_rna_r$date, meta_rna_r$plot, sep = ".")
# 
# #Summarise the count matrices
# rna_r_sum <- t(apply(counts_raw_rna_roots_f, 1, function(f){tapply(f, meta_rna_r$SampleID, mean)}))
# rna_n_sum <- t(apply(counts_raw_rna_needl_f, 1, function(f){tapply(f, meta_rna_n$SampleID, mean)}))
# 
# #Summarise meta tables and adjust summed matrices¨
# meta_r_sum <- unique(meta_rna_r[,-c(1,8)])
# rna_r_sum <- rna_r_sum[,meta_r_sum$SampleID]
# meta_n_sum <- unique(meta_rna_n[,-c(1,8)])
# rna_n_sum <- rna_n_sum[,meta_n_sum$SampleID]
# 
# saveRDS(meta_r_sum, here("RDS/meta_r_sum.rds"))
# saveRDS(meta_n_sum, here("RDS/meta_n_sum.rds"))
# 
# ##Extract control samples only
# control <- grepl("Control", meta_n_sum$treatment)
# 
# rna_n_sum_c <- rna_n_sum[,control]
# rna_n_sum_c <- rna_n_sum_c[rowSums(rna_n_sum_c)>0,]
# rna_r_sum_c <- rna_r_sum[,control]
# rna_r_sum_c <- rna_r_sum_c[rowSums(rna_r_sum_c)>0,]
# 
# saveRDS(rna_n_sum_c, here("RDS/rna_n_sum_c.rds"))
# saveRDS(rna_r_sum_c, here("RDS/rna_r_sum_c.rds"))

#Import needed files
meta_r_sum <- readRDS(here("RDS/meta_r_sum.rds"))
meta_n_sum <- readRDS(here("RDS/meta_n_sum.rds"))
rna_n_sum_c <- readRDS(here("RDS/rna_n_sum_c.rds"))
rna_r_sum_c <- readRDS(here("RDS/rna_r_sum_c.rds"))

#Create list with all transcript names for roots and needles
list_rna <- list()
list_rna[[1]] <- rownames(rna_n_sum_c)
list_rna[[2]] <- rownames(rna_r_sum_c)

pdf(here("Figures/Fig2b_VennControl.pdf"))
par(mfrow=c(1,1))
par(mar=c(0.5,0.5,0.5,0.5))
grid.newpage()
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
grid.draw(venn.diagram(list_rna,
                       filename=NULL,
                       category.names = c("Needles","Roots"),
                       fill=brewer.pal(8,"Dark2")[1:2]))
dev.off()
###Density Curve without pseudocount
#extract log transformed counts
rna_n_sum_c_log1 <- cbind(tissue="Needles", logC=c(log10(rna_n_sum_c)))
rna_r_sum_c_log1 <- cbind(tissue="Roots", logC=c(log10(rna_r_sum_c)))

#Transform for plotting
rna_n_sum_c_log1 <- as.data.frame(rna_n_sum_c_log1)
rna_n_sum_c_log1$logC <- as.numeric(as.character(rna_n_sum_c_log1$logC))
rna_r_sum_c_log1 <- as.data.frame(rna_r_sum_c_log1)
rna_r_sum_c_log1$logC <- as.numeric(as.character(rna_r_sum_c_log1$logC))


#plot
ggplot()+
  #geom_violin(aes(col = tissue, alpha = 0))+
  stat_density(geom = "line", data = rna_n_sum_c_log1, aes(x = logC),
               col = brewer.pal(8, "Dark2")[1], size = 1.1, alpha = .7)+
  stat_density(geom = "line", data = rna_r_sum_c_log1, aes(x = logC),
               col = brewer.pal(8, "Dark2")[2], size = 1.1, alpha = .7)+
  ggformat+
  theme(axis.text.x = element_text(angle = 0))+
  xlab("log10(count)")+
  ylab("Density")+
  coord_cartesian(xlim = c(0, 4))

ggsave(here("Figures/Fig2b_density_RNA.pdf"), width = 3, height = 2.5)

###Density Curve with pseudocount
#extract log transformed counts
rna_n_sum_c_log2 <- cbind(tissue="Needles", logC=c(log10(rna_n_sum_c+1)))
rna_r_sum_c_log2 <- cbind(tissue="Roots", logC=c(log10(rna_r_sum_c+1)))

#Transform for plotting
rna_n_sum_c_log2 <- as.data.frame(rna_n_sum_c_log2)
rna_n_sum_c_log2$logC <- as.numeric(as.character(rna_n_sum_c_log2$logC))
rna_r_sum_c_log2 <- as.data.frame(rna_r_sum_c_log2)
rna_r_sum_c_log2$logC <- as.numeric(as.character(rna_r_sum_c_log2$logC))

#plot
ggplot()+
  #geom_violin(aes(col = tissue, alpha = 0))+
  stat_density(geom = "line", data = rna_n_sum_c_log2, aes(x = logC),
               col = brewer.pal(8, "Dark2")[1], size = 1.1, alpha = .7)+
  stat_density(geom = "line", data = rna_r_sum_c_log2, aes(x = logC),
               col = brewer.pal(8, "Dark2")[2], size = 1.1, alpha = .7)+
  ggformat+
  theme(axis.text.x = element_text(angle = 0))+
  xlab("log10(count+1)")+
  ylab("Density")+
  coord_cartesian(xlim = c(0, 4))

ggsave(here("Figures/FigS2b_density_RNA.pdf"), width = 3, height = 2.5)




