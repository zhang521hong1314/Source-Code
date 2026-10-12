rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("./01_DEG")) {
  dir.create("./01_DEG")
}
setwd("./01_DEG")

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

dat = read.csv('../00_Rawdata1/Dat_Training.csv',
               check.names = F,
               row.names = 1)

group = read.csv('../00_Rawdata1/Group_Training.csv', row.names = 1,
                 check.names = F)

group$group <-
  factor(group$group, levels = c("Control", "Osteosarcoma"))
table(group$group)
dat <- dat[, group$sample]  

dat <- mutate(dat, across(where(is.numeric), as.integer))  

dds <-
  DESeqDataSetFromMatrix(countData = round(dat),
                         colData = group,
                         design = ~ group)

dds <-
  dds[rownames(counts(dds)) > 1,] 
dds <-
  estimateSizeFactors(dds) 

dds <- DESeq(dds)

res <-
  results(dds, contrast = c("group", "Osteosarcoma", "Control"))  
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
    
    "Control" = "blue",
    "Osteosarcoma" = "red"
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
    "Control" = "blue",
    "Osteosarcoma" = "red"
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

save.image("../ws_01_DEG.RData")
