rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("./00_Rawdata1")) {
  dir.create("./00_Rawdata1")
}
setwd("./00_Rawdata1")

library(GEOquery)
library(tidyverse)

GEO_data <- 'GSE99671'

options(timeout = 30000000)  

exprSet <- read.table("GSE99671_Osteosarcoma_filecounts.txt.gz",
                      comment.char = "!",       
                      stringsAsFactors = F,
                      header = T)

expr <- select(exprSet, -Chr)
rownames(expr) <- NULL

write.csv(expr, paste0(GEO_data, '_raw_count.csv'))

annot <- read.csv('./gene_annoation.csv',
                  check.names = F,
                  row.names = 1)

table(annot$GeneType) 
annot <- annot[annot$GeneType == 'protein-coding', ] 

expr <- expr[expr$symbol %in% annot$Symbol, ]

dat <- expr %>%
  dplyr::select(symbol, everything()) %>%     
  mutate(rowMean = rowMeans(.[grep('OS', names(.))])) %>%    
  arrange(desc(rowMean)) %>%       
  distinct(symbol, .keep_all = T) %>%      
  dplyr::select(-rowMean) %>%     
  tibble::column_to_rownames(colnames(.)[1])   

group <- read.csv("./group_GSE99671.csv", check.names = F)
group <- group[group$sample %in% colnames(dat), ] 

group <- group[order(group$group),]  
dat <- dat[, group$sample]  

write.csv(dat, file = paste0('dat.', GEO_data, '_count.csv'))
write.csv(group, file = paste0('group.', GEO_data, '.csv'))
file.symlink(from = paste0('dat.', GEO_data, '_count.csv'),
             to = 'Dat_Training.csv')  
file.symlink(from = paste0('group.', GEO_data, '.csv'), to = 'Group_Training.csv')

colnames(annot)
gene_length <- dplyr::select(annot, Symbol, Length)
gene_length <- gene_length[gene_length$Symbol %in% rownames(dat),]

countToFpkm <- function(counts, effLen) {
  N <- sum(counts)
  exp(log(counts) + log(1e9) - log(effLen) - log(N))
}

dat <- dat[gene_length$Symbol,]
expr_fpkm <- apply(dat, 2, countToFpkm, effLen = gene_length$Length)

range(expr_fpkm)  

if (max(abs(expr_fpkm)) > 50) {
  if (min(expr_fpkm) < 0) {
    expr_fpkm <- log2(expr_fpkm + abs(min(expr_fpkm)) + 1)
    print("log2 transform finished")
  } else if (min(expr_fpkm) == 0) {
    expr_fpkm <- log2(expr_fpkm + 1)
    print("log2 transform finished")
  } else {
    expr_fpkm <- log2(expr_fpkm)
    print("log2 transform finished")
  }
} else{
  print("log2 transform not needed")
}

range(expr_fpkm) 

expr_fpkm <- data.frame(expr_fpkm)
write.csv(expr_fpkm, file = paste0('dat.', GEO_data, '_fpkm.csv'))
file.symlink(from = paste0('dat.', GEO_data, '_fpkm.csv'),
             to = 'Dat_Training_fpkm.csv')  

save.image("../ws_00_Rawdata1.RData")
