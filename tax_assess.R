taxa_its_all_OLD_test <- taxa_its_all_OLD[,-1]
rownames(taxa_its_all_OLD_test) <- colnames(seq.tab_its_all2)

ps_test <- phyloseq(otu_table(readRDS("data/ITS/Both/dada2/seq.tab_all.rds"), taxa_are_rows = FALSE),
                    sample_data(sample_data(read.csv("doc/ITS_all.csv", row.names = 1))),
                    tax_table(as.matrix(taxa_its_all_OLD_test)))

test <- phyloseq_to_metagenomeSeq(ps_its_all_c)
library(metagenomeSeq)
test_gen <- aggregateByTaxonomy(test, lvl = "Genus")
test_mat <- MRcounts(test_gen)
test_mat <- test_mat[!grepl("no_match|unidentified", rownames(test_mat)),]
test_mat2 <- prop.table(test_mat, margin = 2)*100
test_mat2 <- test_mat2[,grepl("root", colnames(test_mat2))]
test_mat2 <- test_mat2[rowMeans(test_mat2)>0,] 
sort(rowMeans(test_mat2), decreasing = TRUE)


old_stuff <- read.csv("~/TP_v2_dump/Haas2018_counts.csv")
old_stuff_tax <- old_stuff[,c(109:115)]
#old_stuff_tax$OTU <- rownames(old_stuff_tax)
old_stuff_tab <- old_stuff[,c(1:108)]
old_stuff_tab_ctrlrt <- old_stuff_tab[,grepl("root.ctrl", colnames(old_stuff_tab))]
colnames(old_stuff_tab_ctrlrt) <- gsub("24th.June", "t2", colnames(old_stuff_tab_ctrlrt))
colnames(old_stuff_tab_ctrlrt) <- gsub("5th.June", "t1", colnames(old_stuff_tab_ctrlrt))
colnames(old_stuff_tab_ctrlrt) <- gsub("6th.August", "t3", colnames(old_stuff_tab_ctrlrt))
colnames(old_stuff_tab_ctrlrt) <- gsub("9th.October", "t4", colnames(old_stuff_tab_ctrlrt))
colnames(old_stuff_tab_ctrlrt) <- paste0("its.", colnames(old_stuff_tab_ctrlrt))

old_meta <- read.csv("doc/ITS_all.csv", row.names = 1)[c(37:40,45:48,53:56),]

ps_old <- phyloseq(otu_table(data.matrix(old_stuff_tab_ctrlrt), taxa_are_rows = TRUE),
                   tax_table(as.matrix(old_stuff_tax)),
                   sample_data(old_meta))
  
janu <- phyloseq_to_metagenomeSeq(ps_old)
janu_gen <- aggregateByTaxonomy(janu, lvl = "Genus")
janu_mat <- MRcounts(janu_gen)
janu_mat <- janu_mat[!grepl("unidentified", rownames(janu_mat)),]
janu_mat2 <- prop.table(janu_mat, margin = 2)*100
janu_mat2 <- janu_mat2[rowMeans(janu_mat2)>0,]
sort(rowMeans(janu_mat2), decreasing = T)

janu_mat3 <- janu_mat2[,colnames(test_mat2)]
