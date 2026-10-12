rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("./01_hypotype")) {
  dir.create("./01_hypotype")
}
setwd("./01_hypotype")

library(lance)
library(readxl)
library(readr)
library(tidyverse)
library(ConsensusClusterPlus)
library(ggplot2)
library(Rtsne)

data1 <-
  read.csv("../00_Rawdata_tcga/Dat_Training_fpkm.csv",
           row.names = 1,
           check.names = F)

group <- read.csv('../00_Rawdata_tcga/Group_Training.csv')
disease <- subset(group, group == 'Osteosarcoma')

gene <- read.csv('../data/PRGs.csv')

data2 <- data1[, disease$sample]
data <- data2[rownames(data2) %in% gene$genes, ]

exprSet = sweep(data, 1, apply(data, 1, median, na.rm = T))

dim(exprSet)

d <- as.matrix(exprSet)

maxK <- 6

clusterAlg <- "pam"


distance <- "canberra"

ccres <- ConsensusClusterPlus(
  d,
  maxK = maxK,
  reps = 1000,
  pItem = 0.9,
  pFeature = 1,
  tmyPal = c("white", "blue"),
  title = '01.cluster/',
  clusterAlg = clusterAlg,
  distance = distance,
  seed = 12345,
  plot = "pdf"
)
ccres_png <- ConsensusClusterPlus(
  d,
  maxK = maxK,
  reps = 1000,
  pItem = 0.9,
  pFeature = 1,
  tmyPal = c("white", "blue"),
  title = '01.cluster/',
  clusterAlg = clusterAlg,
  distance = distance,
  seed = 12345,
  plot = "png"
)

icl = calcICL(ccres, title = '02.PAM_cluster/', plot = "png")
icl_pdf = calcICL(ccres, title = '02.PAM_cluster/', plot = "pdf")

optK <- 2

sample_subtypes <- ccres[[optK]][["consensusClass"]]

table(sample_subtypes)
sample_subtypes <- as.data.frame(sample_subtypes)
write.csv(sample_subtypes, file = "03.sample.subtype.csv")

annCol <- data.frame(results = paste0("Cluster",
                                      ccres[[optK]][['consensusClass']]),
                     row.names = colnames(data2))

set.seed(123)
tSNE_res <-
  Rtsne(
    t(exprSet),
    dims = 2,
    perplexity = 20,
    verbose = F,
    max_iter = 500,
    check_duplicates = F
  )

tsne <-
  data.frame(tSNE1 = tSNE_res[["Y"]][, 1],
             tSNE2 = tSNE_res[["Y"]][, 2],
             cluster = annCol$results)

custom_colors <- c(
  "Cluster1" = "#1d86b5",
  "Cluster2" = "#c76162",
  "Cluster3" = "#476D87",
  "Cluster4" = "#53A85F",
  "Cluster5" = "#F1BB72",
  "Cluster6" = "#F3B1A0",
  "Cluster7" = "#D6E7A3",
  "Cluster8" = "#57C3F3",
  "Cluster9" = "#E5D2DD"
)

p1 <- ggplot(tsne, aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size = 4.5, alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  stat_ellipse(level = 0.85, show.legend = F) +
  theme(legend.position = "top")
ggsave("t-SNE.pdf", plot = p1, width = 6, height = 4)
ggsave("t-SNE.png", plot = p1, width = 6, height = 4, dpi = 300)

survival <- read.csv("../00_Rawdata_tcga/Survival_Training.csv",
                     check.names = F,
                     row.names = 1
)

dat_os1 <- merge(t(data), annCol, by.x = "row.names", by.y = "row.names") %>%
  tibble::column_to_rownames(var = 'Row.names')
dat_os <- merge(dat_os1, survival, by.x = "row.names", by.y = "sample") %>%
  tibble::column_to_rownames(var = 'Row.names')

library(survminer)
library(survival)

kmfit <- survfit(Surv(OS.time, OS) ~ results, data =  dat_os)

train_survival  <- ggsurvplot(
  kmfit,
  pval = TRUE,
  conf.int = F,

  legend.labs = c("Cluster1","Cluster2"),

  legend.title = "Cluster",

  title = "Kaplan-Meier Curve for Survival",

  font.main = c(15, "bold"),
  risk.table = TRUE,

  risk.table.col = "strata",

  linetype = "strata",

  ggtheme = theme_bw(),
  palette = c("#1d86b5", "#c76162", '#476D87', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#E5D2DD')
)

train_survival

pdf(
  file = paste0("KM.pdf"),
  width = 5,
  height = 5,
  family = 'Times'
)
a <- dev.cur()
png(
  file = paste0("KM.png"),
  width = 5,
  height = 5,
  units = "in",
  res = 600,
  family = 'Times'
)
dev.control("enable")
par(mar = c(2, 2, 2, 2),
    cex = 1.5,
    family = "Times")
print(train_survival)
dev.copy(which = a)
dev.off()
dev.off()

save.image("../ws_01_hypotype_1.RData")

rm(list = ls())
setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3/01_hypotype")

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(rstatix)
library(magrittr)
library(stringr)
library(lance)
library(DESeq2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)

dat = read.csv('../00_Rawdata_tcga/Dat_Training_count.csv',
               check.names = F,
               row.names = 1)

group <- read.csv("03.sample.subtype.csv",  check.names = T)
colnames(group) <- c("sample", "group")

group <- group[order(group$group),]
dat <- dat[, group$sample]

group$group <-
  factor(group$group, levels = c("1", "2"))
table(group$group)
dat <- dat[, group$sample]

dat <- mutate(dat, across(where(is.numeric), as.integer))

dds <-
  DESeqDataSetFromMatrix(countData = dat,
                         colData = group,
                         design = ~ group)

dds <-
  dds[rownames(counts(dds)) > 1,]
dds <-
  estimateSizeFactors(dds)

dds <- DESeq(dds)

res <-
  results(dds, contrast = c("group", "2", "1"))

DEG <- as.data.frame(res)
DEG <- na.omit(DEG)

DEG$change <- as.factor(ifelse(
  DEG$padj < 0.05 & abs(DEG$log2FoldChange) > 1,
  ifelse(DEG$log2FoldChange > 1, 'UP', 'DOWN'),
  'NOT'
))

sig_diff <-
  subset(DEG, DEG$padj < 0.05 & abs(DEG$log2FoldChange) > 1)

write.csv(DEG, file = "DEGs_all.csv")
write.csv(sig_diff, file = "DEGs_sig.csv")

dat_rep <-
  DEG[rownames(DEG) %in%
        rownames(rbind(head(sig_diff[order(sig_diff$log2FoldChange, decreasing = T),], 10),
                       head(sig_diff[order(sig_diff$log2FoldChange, decreasing = F),], 10))),]

volcano_plot <- ggplot(data = DEG,
                       aes(
                         x = log2FoldChange,
                         y = -log10(padj),
                         color = change
                       )) +
  scale_color_manual(values = c(

    "UP" = "red",
    "NOT" = "darkgray",
    "DOWN" = "blue"
  )) +
  geom_point(size = 1.2,

             alpha = 0.4,

             na.rm = T) +
  theme_bw(base_size = 12, base_family = "Times") +

  geom_vline(

    xintercept = c(-1, 1),

    lty = 4,

    col = "darkgray",

    lwd = 0.6
  ) +
  geom_hline(

    yintercept = -log10(0.05),
    lty = 4,
    col = "darkgray",
    lwd = 0.6
  ) +
  theme(

    legend.position = "right",

    panel.grid = element_blank(),

    legend.title = element_blank(),

    legend.text = element_text(

      face = "bold",

      color = "black",

      family = "Times",
      size = 13
    ),
    plot.title = element_text(hjust = 0.5),

    axis.text.x = element_text(

      face = "bold",
      color = "black",
      size = 15
    ),
    axis.text.y = element_text(

      face = "bold",
      color = "black",
      size = 15
    ),
    axis.title.x = element_text(

      face = "bold",
      color = "black",
      size = 15
    ),
    axis.title.y = element_text(

      face = "bold",
      color = "black",
      size = 15
    )
  ) +
  geom_label_repel(

    data = dat_rep,

    aes(label = rownames(dat_rep)),

    max.overlaps = 20,

    size = 4,

    box.padding = unit(0.5, "lines"),

    min.segment.length = 0,

    point.padding = unit(0.8, "lines"),

    segment.color = "black",

    show.legend = FALSE
  ) +
  labs(x = "log2(Fold Change)",
       y = "-log10 (P.adjust)")

volcano_plot
ggsave('volcano.png',
       volcano_plot,
       width = 8,
       height = 6)
ggsave('volcano.pdf',
       volcano_plot,
       width = 8,
       height = 6)

library(ComplexHeatmap)
library(circlize)
group <- group %>% as.data.frame()
gene.expr <- dat[rownames(dat_rep),]

mat <- t(scale(t(gene.expr)))

mat[mat < (-2)] <- -2
mat[mat > 2] <- 2

annotation_col <-
  data.frame(row.names = group$sample, Group = group$group)
annotation_row <-
  data.frame(row.names = rownames(dat_rep), sig = dat_rep$log2FoldChange)
annotation_row <-
  annotation_row[order(annotation_row$sig, decreasing = F), , drop = F]

mat <-
  mat[rownames(annotation_row),]
pdf('heatmap.pdf',

    width = 6,
    height = 6,
    family = 'Times')
col_fun <-
  colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

densityHeatmap(mat ,

               title = "Distribution as heatmap",
               ylab = " ",

               height = unit(6, "cm")) %v%
  HeatmapAnnotation(Group = annotation_col$Group, col = list(Group = c(

    "1" = "blue",
    "2" = "red"
  ))) %v%


  Heatmap(

    mat,

    row_names_gp = gpar(fontsize = 9),

    show_column_names = F,

    name = "expression",

    height = unit(6, "cm"),

    cluster_columns = FALSE,

    cluster_rows = FALSE,

    col = col_fun
  )
dev.off()

png(

  'heatmap.png',
  w = 6,
  h = 6,
  units = 'in',

  res = 600,

  family = 'Times'
)
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
densityHeatmap(mat ,
               title = "Distribution as heatmap",
               ylab = " ",
               height = unit(6, "cm")) %v%
  HeatmapAnnotation(Group = annotation_col$Group, col = list(Group = c(
    "1" = "blue",
    "2" = "red"
  ))) %v%
  Heatmap(
    mat,
    row_names_gp = gpar(fontsize = 9),
    show_column_names = F,
    name = "expression",
    height = unit(6, "cm"),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    col = col_fun
  )
dev.off()

save.image("../ws_01_hypotype_2.RData")
