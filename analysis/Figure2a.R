library(here)
library(RColorBrewer)
library(VennDiagram)
library(ggplot2)
source(here("../src/ggplot_format.R"))

#Import metadata
meta_all <- read.csv(here("../doc/Amplicon/meta_NEW.csv"), row.names = 1)

#Clean a bit
meta_all$Treatment <- factor(meta_all$Treatment, levels = c("ctrl",
                                                    "5NO",
                                                    "25NO"))

meta_all$Date <- factor(meta_all$Date, levels = c("5th.June",
                                          "24th.June",
                                          "6th.A ugust",
                                          "9th.October"))

#Import count matrix
count_mat <- read.csv(here("../data/Amplicon/count_mat_NEW.csv"), row.names = 1)

##### ONLY CONTROL SAMPLES ROOTS AND NEEDLES COMBINED
control <- grepl("ctrl", meta_all$SampleID)

count_mat_ctrl <- count_mat[,control]
count_mat_ctrl <- count_mat_ctrl[rowSums(count_mat_ctrl)>0,]

#Booleans to extract roots and needles
roots <- grepl("root", meta_all$SampleID)
needl <- grepl("needl", meta_all$SampleID)

#Extract and filter root samples
count_mat_r <- count_mat[,roots]
count_mat_r <- count_mat_r[rowSums(count_mat_r)>0,]
#Extract and filter needle samples
count_mat_n <- count_mat[,needl]
count_mat_n <- count_mat_n[rowSums(count_mat_n)>0,]

##Venn
##Control only
count_mat_r_c <- count_mat_r[,grepl("ctrl", colnames(count_mat_r))]
count_mat_r_c <- count_mat_r_c[rowSums(count_mat_r_c)>0,]
count_mat_n_c <- count_mat_n[,grepl("ctrl", colnames(count_mat_n))]
count_mat_n_c <- count_mat_n_c[rowSums(count_mat_n_c)>0,]

list_its <- list()
list_its[[1]] <- rownames(count_mat_n_c)
list_its[[2]] <- rownames(count_mat_r_c)

pal=brewer.pal(8,"Dark2")

pdf(here("Figures/Fig2a_VennControl.pdf"))
par(mfrow=c(1,1))
par(mar=c(0.5,0.5,0.5,0.5))
grid.newpage()
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
grid.draw(venn.diagram(list_its,
                       filename=NULL,
                       category.names = c("Needles","Roots"),
                       fill=pal[1:2]))
dev.off()

##Density curve without pseudocount
#Extract log transformed counts
count_mat_r_c_log1 <- cbind(tissue="Roots", logC=c(log10(data.matrix(count_mat_r_c))))
count_mat_n_c_log1 <- cbind(tissue="Needles", logC=c(log10(data.matrix(count_mat_n_c))))

count_mat_r_c_log1 <- as.data.frame(count_mat_r_c_log1)
count_mat_r_c_log1$logC <- as.numeric(as.character(count_mat_r_c_log1$logC))

count_mat_n_c_log1 <- as.data.frame(count_mat_n_c_log1)
count_mat_n_c_log1$logC <- as.numeric(as.character(count_mat_n_c_log1$logC))

ggplot()+
  stat_density(geom = "line", data = count_mat_n_c_log1, aes(x = logC),
               col = brewer.pal(8, "Dark2")[1], size = 1.1, alpha = .7)+
  stat_density(geom = "line", data = count_mat_r_c_log1, aes(x = logC),
               col = brewer.pal(8, "Dark2")[2], size = 1.1, alpha = .7)+
  ggformat+
  theme(axis.text.x = element_text(angle = 0),
        axis.title.x = element_text())+
  xlab("log10(count)")+
  ylab("Density")+
  coord_cartesian(xlim = c(0, 4))

ggsave(here("Figures/Fig2a_density_ITS.pdf"), width = 3, height = 2.5)

##Density curve with pseudocount
#Extract log transformed counts
count_mat_r_c_log <- cbind(tissue="Roots", logC=c(log10(data.matrix(count_mat_r_c+1))))
count_mat_n_c_log <- cbind(tissue="Needles", logC=c(log10(data.matrix(count_mat_n_c+1))))

count_mat_r_c_log <- as.data.frame(count_mat_r_c_log)
count_mat_r_c_log$logC <- as.numeric(as.character(count_mat_r_c_log$logC))

count_mat_n_c_log <- as.data.frame(count_mat_n_c_log)
count_mat_n_c_log$logC <- as.numeric(as.character(count_mat_n_c_log$logC))

ggplot()+
  stat_density(geom = "line", data = count_mat_n_c_log, aes(x = logC),
               col = brewer.pal(8, "Dark2")[1], size = 1.1, alpha = .7)+
  stat_density(geom = "line", data = count_mat_r_c_log, aes(x = logC),
               col = brewer.pal(8, "Dark2")[2], size = 1.1, alpha = .7)+
  ggformat+
  theme(axis.text.x = element_text(angle = 0),
        axis.title.x = element_text())+
  xlab("log10(count+1)")+
  ylab("Density")+
  coord_cartesian(xlim = c(0, 4))

ggsave(here("Figures/FigS2a_density_ITS.pdf"), width = 3, height = 2.5)




