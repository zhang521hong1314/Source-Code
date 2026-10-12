rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("./00_Rawdata2")) {
  dir.create("./00_Rawdata2")
}
setwd("./00_Rawdata2")

library(GEOquery)
library(tidyverse)
library(lance)
library(readr)
library(readxl)  
library(stringr)
library(AnnoProbe)

options(timeout = 300000000)  

gset <- AnnoProbe::geoChina('GSE21257')  

gset <- gset[[1]]  

expr <-
  as.data.frame(exprs(gset))  

gpl <- getGEO(filename = "./GPL10295_family.soft.gz") %>% Table()

probe2symbol <- gpl %>%
  dplyr::select('ID', contains("ymbol"))  

colnames(probe2symbol) <- c('ID', 'Gene_name')  

probe2symbol <- probe2symbol %>%
  filter('Gene_name' != '') %>%  
  separate('Gene_name', c('symbol', 'drop'), sep = '///') %>%  
  dplyr::select(-drop)  

probe2symbol = probe2symbol[probe2symbol$symbol != '',]  
probe2symbol$ID <- as.character(probe2symbol$ID)  

range(expr)  

if (max(abs(expr)) > 50) {
  if (min(expr) < 0) {
    expr <- log2(expr + abs(min(expr)) + 1)
    print("log2 transform finished")
  } else if (min(expr) == 0) {
    expr <- log2(expr + 1)
    print("log2 transform finished")
  } else {
    expr <- log2(expr)
    print("log2 transform finished")
  }
} else{
  print("log2 transform not needed")
}

range(expr) 

dat <- expr
dat$ID <- rownames(dat)  
dat$ID <- as.character(dat$ID)  

dat <- dat %>%
  inner_join(probe2symbol, by = 'ID') %>%  
  dplyr::select(-ID) %>%  
  
  dplyr::select(symbol, everything()) %>%  
  mutate(rowMean = rowMeans(.[grep('GSM', names(.))])) %>%
  
  
  
  
  
  arrange(desc(rowMean)) %>%  
  distinct(symbol, .keep_all = T) %>%  
  dplyr::select(-rowMean) %>%  
  tibble::column_to_rownames(colnames(.)[1])  

protein <-
  read.table(
    "./protein_coding.txt",
    header = T,
    sep = "\t"
  )

dat <-
  dat[rownames(dat) %in% protein$gene_name, ]  

pd <- pData(gset)  

colnames(pd)  

group <- data.frame(sample = pd$geo_accession,
                    group = pd$source_name_ch1)  
table(group$group)  

group$group <- ifelse(
  grepl("hhhhhhhh", group$group),
  'Control',
  ifelse(
    grepl("Osteosarcoma", group$group),
    'Osteosarcoma',
    'Other_Disease'
  )
)  

table(group$group)

group <- group[group$group != "Other_Disease", ]  
group <- group[order(group$group), ]  
dat <- dat[, group$sample]  

survival <- subset(pd, select = c('characteristics_ch1.6', 'characteristics_ch1.6'))

survival <-
  dplyr::rename(survival, "OS.time" = "characteristics_ch1.6", "OS" = "characteristics_ch1.6.1")

survival <- mutate(survival, OS.time = as.numeric(gsub(".*(at|after)\\s*([0-9]+)\\s*months.*", "\\2", OS.time)))  

survival$OS <- ifelse(
  grepl("Deceased", survival$OS),
  '1',
  ifelse(
    grepl("Alive", survival$OS),
    '0',
    'Other'
  )
)  

survival[, c("OS", "OS.time")] <- lapply(survival[, c("OS", "OS.time")], as.numeric)  
survival$OS.time <- round(survival$OS.time * 30)  
survival <- na.omit(survival)

survival <- subset(survival, survival$OS.time > 0)

str(survival)  

tumor_sample <- group[group$group == "Osteosarcoma", ]

tumor_sample <- tumor_sample[tumor_sample$sample %in% rownames(survival), ]
rownames(tumor_sample) <- seq_len(nrow(tumor_sample))

survival <- survival[tumor_sample$sample, ]
survival <- cbind(sample = rownames(survival), survival) 
rownames(survival) <- seq_len(nrow(survival))

group <- tumor_sample  

dat <- dat[, group$sample]  

write.csv(dat, file = 'Dat_Training.csv')  
write.csv(group, file = 'Group_Training.csv')
write.csv(survival, file = 'Survival_Training.csv')
write.csv(tumor_sample, file = 'Tumor_Sample.csv')

save.image("../ws_00_Rawdata2.RData")

