library(here)
library(DESeq2)
library(ggplot2)
source(here("../src/ggplot_format.R"))
source(here("../src/Rtoolbox/old/utilsDE.r"))
source(here("../src/Rtoolbox/old/plotUpAndDown.R"))
source(here("../src/UPSCb-common/src/R/gopher.R"))
source(here("../src/Rtoolbox/old/plotEnrichedTreemap.R"))

#Import KEGG ortholog data
kos_r_sum <- readRDS(here("RDS/kos_r_sum.rds"))

meta_r_sum <- readRDS(here("RDS/meta_r_sum.rds"))
rownames(meta_r_sum) <- meta_r_sum$SampleID
meta_r_sum$Group2 <- as.factor(paste(meta_r_sum$treatment, meta_r_sum$date, sep="_"))
meta_r_sum$Group2 <- factor(meta_r_sum$Group2, levels = c("Control_Early_June","Control_Late_June","Control_August","Control_October",
                                                          "5_year_Early_June",  "5_year_Late_June", "5_year_August", "5_year_October",
                                                          "25_year_Early_June","25_year_Late_June", "25_year_August", "25_year_October"))


dds_sum_r <- DESeqDataSetFromMatrix(round(kos_r_sum), colData = meta_r_sum, design = ~Group2)
dds_sum_r <- DESeq(dds_sum_r)
resultsNames(dds_sum_r)

controlVector <- as.character(unique(meta_r_sum$Group2[grep("Control", meta_r_sum$Group2)]))
fertilised5Vector <- as.character(unique(meta_r_sum$Group2[grep("^5", meta_r_sum$Group2)]))
fertilised25Vector <- as.character(unique(meta_r_sum$Group2[grep("25", meta_r_sum$Group2)]))



fert5_vs_control.res <- mapply(getRes, fertilised5Vector, 
                               controlVector, 
                               MoreArgs = list(localDDS=dds_sum_r, group="Group2"))

fert25_vs_control.res <- mapply(getRes, fertilised25Vector, 
                                controlVector, 
                                MoreArgs = list(localDDS=dds_sum_r, group="Group2"))

fert25_vs_fert5.res <- mapply(getRes, fertilised25Vector, 
                              fertilised5Vector, 
                              MoreArgs = list(localDDS=dds_sum_r, group="Group2"))


names(fert5_vs_control.res) <- c("5_year_Early_June vs Control_Early_June", 
                                 "5_year_Late_June vs Control_Late_June",                        
                                 "5_year_August vs Control_August", 
                                 "5_year_October vs Control_October")

names(fert25_vs_control.res) <- c("25_year_Early_June vs Control_Early_June", 
                                  "25_year_Late_June vs Control_Late_June",                        
                                  "25_year_August vs Control_August", 
                                  "25_year_October vs Control_October")

names(fert25_vs_fert5.res) <- c("25_year_Early_June vs 5_year_Early_June", 
                                "25_year_Late_June vs 5_year_Late_June",                        
                                "25_year_August vs 5_year_August", 
                                "25_year_October vs 5_year_October")


fert5_vs_control.res.filter <- lapply(fert5_vs_control.res, filterDE, p=0.05)
fert25_vs_control.res.filter <- lapply(fert25_vs_control.res, filterDE, p=0.05)
fert25_vs_fert5.res.filter <- lapply(fert25_vs_fert5.res, filterDE, p=0.05)

names(fert5_vs_control.res.filter) <- c("5T1 vs CT1", 
                                        "5T2 vs CT2",                        
                                        "5T3 vs CT3", 
                                        "5T4 vs CT4")

names(fert25_vs_control.res.filter) <- c("25T1 vs CT1", 
                                         "25T2 vs CT2",                        
                                         "25T3 vs CT3", 
                                         "25T4 vs CT4")

names(fert25_vs_fert5.res.filter) <- c("25T1 vs 5T1", 
                                       "25T2 vs 5T2",                        
                                       "25T3 vs 5T3", 
                                       "25T4 vs 5T4")



#Plot and save the 25 year vs control comparison
plotUpAndDown(fert25_vs_control.res.filter)+
  ylab("Number of orthologs")+
  scale_x_discrete(labels=c("Early June", "Late June", "August", "October"))+
  #coord_cartesian(ylim = c(-750, 750))+
  ggtitle("25 years vs Control")+
  ggformat
ggsave(here("Figures/Fig7a_25vsC.pdf"), width = 3.6, height = 5)

##Investigate how many KOs are common between timepoints
fert25_vs_C.T1 <- as.data.frame(fert25_vs_control.res.filter[[1]])
fert25_vs_C.T1 <- cbind(KO=rownames(fert25_vs_C.T1), fert25_vs_C.T1)
fert25_vs_C.T2 <- as.data.frame(fert25_vs_control.res.filter[[2]])
fert25_vs_C.T2 <- cbind(KO=rownames(fert25_vs_C.T2), fert25_vs_C.T2)
fert25_vs_C.T3 <- as.data.frame(fert25_vs_control.res.filter[[3]])
fert25_vs_C.T3 <- cbind(KO=rownames(fert25_vs_C.T3), fert25_vs_C.T3)
fert25_vs_C.T4 <- as.data.frame(fert25_vs_control.res.filter[[4]])
fert25_vs_C.T4 <- cbind(KO=rownames(fert25_vs_C.T4), fert25_vs_C.T4)

###From this to a treemap
fert25_vs_C.ALLtp <- rbind(fert25_vs_C.T1,
                           fert25_vs_C.T2,
                           fert25_vs_C.T3,
                           fert25_vs_C.T4)

##Total unique KOs
length(unique(fert25_vs_C.ALLtp$KO))
#common between all timepoints
length(intersect(intersect(intersect(fert25_vs_C.T1$KO, fert25_vs_C.T2$KO), fert25_vs_C.T3$KO), fert25_vs_C.T4$KO))
#Extract common KOs 
kos_common_all <- intersect(intersect(intersect(fert25_vs_C.T1$KO, fert25_vs_C.T2$KO), fert25_vs_C.T3$KO), fert25_vs_C.T4$KO)

#Separate up and down
fert25_vs_C.ALLtp.UP <- fert25_vs_C.ALLtp[fert25_vs_C.ALLtp$log2FoldChange>0,]
fert25_vs_C.ALLtp.DOWN <- fert25_vs_C.ALLtp[fert25_vs_C.ALLtp$log2FoldChange<0,]

#Extract KOs to be tested for enrichment
fert25_vs_C.ALLtp.UP_names <- unique(fert25_vs_C.ALLtp.UP$KO)
fert25_vs_C.ALLtp.DOWN_names <- unique(fert25_vs_C.ALLtp.DOWN$KO)

#Save RDS objects
ifelse(!dir.exists(here("RDS")), dir.create(here("RDS")), FALSE)
saveRDS(fert25_vs_C.ALLtp.UP, here("RDS/fert25_vs_C.ALLtp.UP.rds"))
saveRDS(fert25_vs_C.ALLtp.DOWN, here("RDS/fert25_vs_C.ALLtp.DOWN.rds"))

#Run gofer to enrich
enr.UP <- gopher(fert25_vs_C.ALLtp.UP_names, task = list("ko_pathway"), alpha=0.05, url = "ko", endpoint = "enrichment")
enr.DOWN <- gopher(fert25_vs_C.ALLtp.DOWN_names, task = list("ko_pathway"), alpha=0.05, url = "ko", endpoint = "enrichment")

pdf(here("Figures/Fig7b_TreemapUp.pdf"), width = 10, height = 7)
plotEnrichedTreemap(enr.UP, enrichment = "ko_pathway", de = "up", sizeCol = "nt")
dev.off()

pdf(here("Figures/Fig7b_TreemapDown.pdf"), width = 10, height = 7)
plotEnrichedTreemap(enr.DOWN, enrichment = "ko_pathway", de = "down", sizeCol = "nt")
dev.off()



