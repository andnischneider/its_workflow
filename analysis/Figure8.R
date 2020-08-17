library(here)
library(DESeq2)
library(ggplot2)
source(here("../src/Rtoolbox/old/utilsDE.r"))
source(here("../src/Rtoolbox/old/plotUpAndDown.R"))
source(here("../src/UPSCb-common/src/R/gopher.R"))
source(here("../src/Rtoolbox/old/plotEnrichedTreemap.R"))
source(here("../src/ggplot_format.R"))


####DE analysis of the genus Cortinarius
#Shows strong response

#Import count matrix and taxonomic annotations (from Fig 7)
mat_r_sum <- readRDS(here("RDS/rna_r_sum.rds"))
meta_r_sum <- readRDS(here("RDS/meta_r_sum.rds"))

tax <- readRDS(here("RDS/annot_tax_filt.rds"))
#How many Cortinarius transcripts
nrow(tax[tax$genus=="Cortinarius",])
tax_cort <- tax[tax$genus=="Cortinarius",]
rownames(tax_cort) <- tax_cort$gene


#subset to Cortinarius
mat_r_sum_cort <- mat_r_sum[rownames(tax_cort),]

#Add factor to meta table
rownames(meta_r_sum) <- meta_r_sum$SampleID
meta_r_sum$Group2 <- as.factor(paste(meta_r_sum$treatment, meta_r_sum$date, sep="_"))
meta_r_sum$Group2 <- factor(meta_r_sum$Group2, levels = c("Control_Early_June","Control_Late_June","Control_August","Control_October",
                                                          "5_year_Early_June",  "5_year_Late_June", "5_year_August", "5_year_October",
                                                         "25_year_Early_June","25_year_Late_June", "25_year_August", "25_year_October"))

#DESeq
dds_sum_r <- DESeqDataSetFromMatrix(round(mat_r_sum_cort), colData = meta_r_sum, design = ~Group2)
dds_sum_r <- DESeq(dds_sum_r)
resultsNames(dds_sum_r)

controlVector <- as.character(unique(meta_r_sum$Group2[grep("Control", meta_r_sum$Group2)]))
fertilised5Vector <- as.character(unique(meta_r_sum$Group2[grep("^5", meta_r_sum$Group2)]))
fertilised25Vector <- as.character(unique(meta_r_sum$Group2[grep("25", meta_r_sum$Group2)]))


fert25_vs_control.res <- mapply(getRes, fertilised25Vector, 
                                controlVector, 
                                MoreArgs = list(localDDS=dds_sum_r, group="Group2"))



names(fert25_vs_control.res) <- c("25_year_Early_June vs Control_Early_June", 
                                  "25_year_Late_June vs Control_Late_June",                        
                                  "25_year_August vs Control_August", 
                                  "25_year_October vs Control_October")

fert25_vs_control.res.filter <- lapply(fert25_vs_control.res, filterDE, p=0.05)


names(fert25_vs_control.res.filter) <- c("25T1 vs CT1", 
                                         "25T2 vs CT2",                        
                                         "25T3 vs CT3", 
                                         "25T4 vs CT4")


plotUpAndDown(fert25_vs_control.res.filter)+
  ylab("Number of orthologs")+
  scale_x_discrete(labels=c("Early June", "Late June", "August", "October"))+
  #coord_cartesian(ylim = c(-750, 750))+
  ggtitle("25 years vs Control")+
  ggformat
ggsave(here("Figures/Fig8a_Cortinarius_25vsC.pdf"), width = 3.6, height = 5)

#Save as separate objects and add transcripts as extra columns, to keep duplicates
fert25_vs_C.T1 <- as.data.frame(fert25_vs_control.res.filter[[1]])
fert25_vs_C.T1 <- cbind(transcript=rownames(fert25_vs_C.T1), fert25_vs_C.T1)
fert25_vs_C.T2 <- as.data.frame(fert25_vs_control.res.filter[[2]])
fert25_vs_C.T2 <- cbind(transcript=rownames(fert25_vs_C.T2), fert25_vs_C.T2)
fert25_vs_C.T3 <- as.data.frame(fert25_vs_control.res.filter[[3]])
fert25_vs_C.T3 <- cbind(transcript=rownames(fert25_vs_C.T3), fert25_vs_C.T3)
fert25_vs_C.T4 <- as.data.frame(fert25_vs_control.res.filter[[4]])
fert25_vs_C.T4 <- cbind(transcript=rownames(fert25_vs_C.T4), fert25_vs_C.T4)

###From this to a treemap
#Merge all 4
fert25_vs_C.ALLtp <- rbind(fert25_vs_C.T1,
                           fert25_vs_C.T2,
                           fert25_vs_C.T3,
                           fert25_vs_C.T4)

##Total unique KOs
length(unique(fert25_vs_C.ALLtp$transcript))
#common between all timepoints
length(intersect(intersect(intersect(fert25_vs_C.T1$transcript, fert25_vs_C.T2$transcript), fert25_vs_C.T3$transcript), fert25_vs_C.T4$transcript))
#Extract common KOs 
genes_common_all <- intersect(intersect(intersect(fert25_vs_C.T1$transcript, fert25_vs_C.T2$transcript), fert25_vs_C.T3$transcript), fert25_vs_C.T4$transcript)

#Separate up and down
fert25_vs_C.ALLtp.UP <- fert25_vs_C.ALLtp[fert25_vs_C.ALLtp$log2FoldChange>0,]
fert25_vs_C.ALLtp.DOWN <- fert25_vs_C.ALLtp[fert25_vs_C.ALLtp$log2FoldChange<0,]

#Extract genes
fert25_vs_C.ALLtp.UP_names <- as.character(fert25_vs_C.ALLtp.UP$transcript)
fert25_vs_C.ALLtp.DOWN_names <- as.character(fert25_vs_C.ALLtp.DOWN$transcript)

#Import functional annotations and add matching KOs to genes (where available)
annot_fun_all <- read.delim(here("../data/RNA/annotation_results.emapper.annotations"), header = FALSE)
annot_fun <- annot_fun_all[annot_fun_all$V1%in%tax$gene,]
annot_fun$V9 <- gsub("ko:", "", annot_fun$V9)

#Turn Cortinarius genes into matching KOs
kos_up <- annot_fun$V9[annot_fun$V1%in%fert25_vs_C.ALLtp.UP_names]
kos_down <- annot_fun$V9[annot_fun$V1%in%fert25_vs_C.ALLtp.DOWN_names]
#process to remove empty elements 
kos_up <- kos_up[kos_up!=""]
kos_down <- kos_down[kos_down!=""]
#process to separate the multiple hits
kos_up <- unique(strsplit(paste(kos_up, collapse = ","), ",")[[1]])
kos_down <- unique(strsplit(paste(kos_down, collapse = ","), ",")[[1]])

kos_cort <- annot_fun$V9[annot_fun$V1%in%tax_cort$gene]
kos_cort <- kos_cort[kos_cort!=""]
kos_cort <- unique(strsplit(paste(kos_cort, collapse = ","), ",")[[1]])


enr.UP <- gopher(kos_up, task = list("ko_pathway"), alpha=0.05, url = "ko", endpoint = "enrichment")
enr.DOWN <- gopher(kos_down, task = list("ko_pathway"), alpha=0.001, url = "ko", endpoint = "enrichment")

pdf(here("Figures/Fig8b_TreemapUp.pdf"), width = 10, height = 5)
plotEnrichedTreemap(enr.UP, enrichment = "ko_pathway", de = "up", sizeCol = "nt")
dev.off()

pdf(here("Figures/Fig8b_TreemapDown.pdf"), width = 10, height = 5)
plotEnrichedTreemap(enr.DOWN, enrichment = "ko_pathway", de = "down", sizeCol = "nt")
dev.off()




