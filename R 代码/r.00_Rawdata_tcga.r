rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("./00_Rawdata_tcga")) {
  dir.create("./00_Rawdata_tcga")
}
setwd("./00_Rawdata_tcga")

library(tidyverse)
library(dplyr)
library(rtracklayer)

tcga_count <-
  read_tsv(file = './TARGET-OS.htseq_counts.tsv.gz') 
tcga_fpkm <-
  read_tsv(file = "./TARGET-OS.htseq_fpkm.tsv.gz") 
tcga_survival <-
  read_tsv(file = "./TARGET-OS.survival.tsv") 
tcga_phenotype <-
  read_tsv(file = "./TARGET-OS.clinical.tsv.gz") 

tcga_count <- as.data.frame(tcga_count)  

tcga_count <-
  column_to_rownames(tcga_count, var = 'Ensembl_ID') 
tcga_count <- 2 ^ tcga_count - 1 

gene_annotation <- import('./gencode.v22.annotation.gtf.gz')

gene_annotation <- as.data.frame(gene_annotation)
gene_annotation <-
  gene_annotation[gene_annotation$type == 'gene', ] 
gene_annotation <-
  dplyr::select(gene_annotation,
                gene_id,
                gene_name,
                seqnames,
                start,
                end,
                width,
                strand,
                gene_type)

gene_symbol <-
  gene_annotation[, (1:2)] 
colnames(gene_symbol) <-
  c("ID", "symbol") 

data_tcga <- tcga_count
data_tcga$ID <- rownames(data_tcga) 
data_tcga$ID <-
  as.character(data_tcga$ID) 
gene_symbol$ID <-
  as.character(gene_symbol$ID) 

data_tcga <-
  inner_join(data_tcga, gene_symbol, by = "ID") 
data_tcga <- dplyr::select(data_tcga, -ID) 
data_tcga <-
  dplyr::select(data_tcga, symbol, everything()) 
data_tcga <-
  mutate(data_tcga, rowMean = rowMeans(data_tcga[grep('TCGA', names(data_tcga))])) 
data_tcga <- arrange(data_tcga, desc(rowMean)) 
data_tcga <-
  distinct(data_tcga, symbol, .keep_all = T) 
data_tcga <- dplyr::select(data_tcga, -rowMean) 
data_tcga <-
  tibble::column_to_rownames(data_tcga, var = "symbol")   
dim(data_tcga)  

mete <- data.frame(colnames(data_tcga))
for (i in 1:length(mete[, 1])) {
  num <-
    as.numeric(as.character(substring(mete[i, 1], 18, 19))) 
  if (num %in% seq(1, 9)) {
    mete[i, 2] = "T"
  } 
  if (num %in% seq(10, 29)) {
    mete[i, 2] = "N"
  } 
}
colnames(mete) <- c("id", "group")
table(mete$group)
mete$group <- as.factor(mete$group) 
mete <-
  subset(mete, mete$group == "T") 

exp_tumor <-
  data_tcga[, which(colnames(data_tcga) %in% mete$id)] 
exp_tumor <- as.data.frame(exp_tumor)

exp_control <-
  data_tcga[, which(!colnames(data_tcga) %in% mete$id)] 
exp_control <- as.data.frame(exp_control)

data_final <- cbind(exp_control, exp_tumor) 
data_final <- na.omit(data_final) 

clinical_data <- tcga_phenotype[tcga_phenotype$sample_id %in% colnames(data_final), ]
Group <- data.frame(sample = colnames(data_final), group = "Osteosarcoma")

table(Group$group)

dim(data_final)  

tcga_survival <- tcga_survival[, -3]
tcga_survival <- tcga_survival[tcga_survival$sample %in% Group$sample, ]

data_final <- data_final[, colnames(data_final) %in% tcga_survival$sample]
Group <- Group[Group$sample %in% tcga_survival$sample, ]
rownames(Group) <- seq_len(nrow(Group))

dim(data_final)

tumor_sample <- Group[Group$group == "Osteosarcoma", ]
rownames(tumor_sample) <- seq_len(nrow(tumor_sample))

write.csv(tcga_survival, './Survival_Training.csv')
write.csv(clinical_data, './Clinical_Training.csv')
write.csv(data_final, file = './Dat_Training_count.csv')
write.csv(Group, file = './Group_Training.csv')
write.csv(tumor_sample, file = './Tumor_Sample.csv')

tcga_fpkm <- as.data.frame(tcga_fpkm)
tcga_fpkm <- tcga_fpkm[!duplicated(tcga_fpkm$Ensembl_ID),]
tcga_fpkm <- column_to_rownames(tcga_fpkm, var = 'Ensembl_ID')
tcga_fpkm <- 2 ^ tcga_fpkm - 1

dat_fpkm <- tcga_fpkm
dat_fpkm$ID <- rownames(dat_fpkm)
dat_fpkm$ID <- as.character(dat_fpkm$ID)

gene_symbol$ID <- as.character(gene_symbol$ID)

dat_fpkm <- dat_fpkm %>%
  inner_join(gene_symbol, by = 'ID') %>%
  select(-ID) %>%     
  select(symbol, everything()) %>%     
  mutate(rowMean = rowMeans(.[grep('TCGA', names(.))])) %>%    
  arrange(desc(rowMean)) %>%       
  distinct(symbol, .keep_all = T) %>%      
  select(-rowMean) %>%     
  tibble::column_to_rownames(colnames(.)[1])   
dim(dat_fpkm)  

dat_fpkm <- dat_fpkm[rownames(data_final), ]
dat_fpkm <- dat_fpkm[,colnames(data_final)]
dat_fpkm <- na.omit(dat_fpkm)
dat_fpkm <- log2(dat_fpkm + 1)

dim(dat_fpkm)  

range(dat_fpkm)  

if (max(abs(dat_fpkm)) > 50) {
  if (min(dat_fpkm) < 0) {
    dat_fpkm <- log2(dat_fpkm + abs(min(dat_fpkm)) + 1)
    print("log2 transform finished")
  } else if (min(dat_fpkm) == 0) {
    dat_fpkm <- log2(dat_fpkm + 1)
    print("log2 transform finished")
  } else {
    dat_fpkm <- log2(dat_fpkm)
    print("log2 transform finished")
  }
} else{
  print("log2 transform not needed")
}

range(dat_fpkm) 

write.csv(dat_fpkm, file = './Dat_Training_fpkm.csv')

save.image("../ws_00_Rawdata_tcga.RData")
