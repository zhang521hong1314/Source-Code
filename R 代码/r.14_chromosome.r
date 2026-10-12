rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("./14_chromosome")) {
  dir.create("./14_chromosome")
}
setwd("./14_chromosome")

library(OmicCircos)
library(S4Vectors)
library(tidyverse)
library(rtracklayer)
library(RCircos)

DEG <- read.csv('../01_DEG/DEGs_all.csv', row.names = 1, check.names = F)
gene <- read.csv('../06_model_validation/lasso_genes.csv', row.names = 1, check.names = F)
colnames(gene)[1] <- 'genes'

dat_rep <- DEG[rownames(DEG) %in% gene$genes, ]
gene <-
  data.frame(Symbol = rownames(dat_rep), change = dat_rep$change)
gene$colour <- ifelse(gene$change == 'UP', 'red', 'blue')

hg19.gencode.gtf <-
  import.gff(con = './gencode.v19.annotation.gtf.gz')
chr <- data.frame(hg19.gencode.gtf)
chr <- subset(chr, source == 'HAVANA' & type == 'gene')  
chr <- chr[chr$gene_name %in% gene$Symbol,]
chr <- chr[chr$gene_status == 'KNOWN',]

RCircos.Gene.Label <- data.frame(chr$seqnames,
                                 chr$start,
                                 chr$end,
                                 chr$gene_name)
data(RCircos.Gene.Label.Data)
colnames(RCircos.Gene.Label) <- colnames(RCircos.Gene.Label.Data)

RCircos.Gene.Label$pos <- NA
cnt <- 1
while (cnt < nrow(chr) + 1) {
  RCircos.Gene.Label$pos[cnt] <-
    mean(c(
      RCircos.Gene.Label$chromStart[cnt],
      RCircos.Gene.Label$chromEnd[cnt]
    ))
  cnt = cnt + 1
}
gene_fus <- data.frame(
  chr = RCircos.Gene.Label$Chromosome,
  pos = RCircos.Gene.Label$pos,
  Symbol = RCircos.Gene.Label$Gene,
  gene = RCircos.Gene.Label$Gene
)

gene_chr <- merge(gene_fus, gene, by = 'Symbol')
gene_chr <- gene_chr[, -1]
colnames(gene_chr)[3] <- 'Symbol'
gene_chr <- gene_chr[order(gene_chr$pos, decreasing = F), ]
gene_chr <- gene_chr[order(gene_chr$chr, decreasing = F), ]
chrom <- unique(gene_fus$chr)
gene_chr$chr <- factor(gene_chr$chr, levels = chrom)

expr <- read.csv('../00_Rawdata1/Dat_Training_fpkm.csv', row.names = 1, check.names = F)
group <-
  read.csv(file = '../00_Rawdata1/Group_Training.csv', row.names = 1, check.names = F)

x <- expr
mat <- t(scale(t(x)))

mat[mat < (-2)] <- (-2)
mat[mat > 2] <- 2
expr <- mat

skcm <- subset(group, group == 'AP')
exp_gene_skcm <-
  subset(expr, select = skcm$sample) %>% as.data.frame() 
exp_gene_skcm$Symbol <- rownames(exp_gene_skcm)
chr_heatmap1 <- merge(gene_fus, exp_gene_skcm, by = 'Symbol')
chr_heatmap1 <- chr_heatmap1[, -1]

control <- subset(group, group == 'Control')
exp_gene_control <-
  subset(expr, select = control$sample) %>% as.data.frame() 
exp_gene_control$Symbol <- rownames(exp_gene_control)
chr_heatmap2 <- merge(gene_fus, exp_gene_control, by = 'Symbol')
chr_heatmap2 <- chr_heatmap2[, -1]

logfc <- DEG
logfc <- data.frame(Symbol = rownames(DEG), logfc = DEG$log2FoldChange)
logfc <- merge(gene_chr, logfc, by = 'Symbol')
logfc <- logfc[, c('chr', 'pos', 'Symbol', 'logfc', 'colour')]
logfc <- logfc[order(logfc$pos, decreasing = F), ]
logfc <- logfc[order(logfc$chr, decreasing = F), ]
chrom <- unique(gene_fus$chr)
logfc$chr <- factor(logfc$chr, levels = chrom)

pdf(
  'gene_OmicCircos.pdf',
  w = 9,
  h = 9,
  family = 'Times'
)
par(
  mfrow = c(1, 1),
  pin = c(5, 5),
  mar = c(6, 6, 6, 5),
  family = 'Times',
  font = 2
)
plot(
  c(1, 800),
  c(1, 800),
  type = 'n',
  axes = F,
  xlab = '',
  ylab = '',
  main = ''
)
circos(
  R = 330,
  cir = 'hg18',
  W = 40,
  mapping = gene_chr,
  type = 'label',
  side = 'out',
  col = gene_chr$colour,
  cex = 0.8,
  lwd = 1
)
circos(
  R = 300,
  cir = 'hg18',
  W = 10,
  type = 'chr',
  print.chr.lab = T,
  scale = T,
  cex = 10
)  

circos(
  R = 140,
  cir = 'hg18',
  W = 80,
  mapping = chr_heatmap2,
  col.v = 4,
  type = 'heatmap2',
  cluster = F,
  col.bar = F,
  lwd = -1,
  col = 'black'
)
circos(
  R = 90,
  cir = "hg18",
  W = 50,
  mapping = logfc,
  col.v = 4,
  type = "s",
  B = TRUE,
  col = logfc$colour
)

dev.off()

png(
  'gene_OmicCircos.png',
  w = 9,
  h = 9,
  units = 'in',
  res = 600,
  family = 'Times'
)
par(
  mfrow = c(1, 1),
  pin = c(5, 5),
  mar = c(6, 6, 6, 5),
  family = 'Times',
  font = 2
)
plot(
  c(1, 800),
  c(1, 800),
  type = 'n',
  axes = F,
  xlab = '',
  ylab = '',
  main = ''
)
circos(
  R = 330,
  cir = 'hg18',
  W = 40,
  mapping = gene_chr,
  type = 'label',
  side = 'out',
  col = gene_chr$colour,
  cex = 0.8,
  lwd = 1
)
circos(
  R = 300,
  cir = 'hg18',
  W = 10,
  type = 'chr',
  print.chr.lab = T,
  scale = T,
  cex = 10
)  

circos(
  R = 140,
  cir = 'hg18',
  W = 80,
  mapping = chr_heatmap2,
  col.v = 4,
  type = 'heatmap2',
  cluster = F,
  col.bar = F,
  lwd = -1,
  col = 'black'
)
circos(
  R = 90,
  cir = "hg18",
  W = 50,
  mapping = logfc,
  col.v = 4,
  type = "s",
  B = TRUE,
  col = logfc$colour
)

dev.off()

save.image("../ws_14_chromosome.RData")
