rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("./10_GSEA")) {
  dir.create("./10_GSEA")
}
setwd("./10_GSEA")

library(clusterProfiler)
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

tarin_dat <- read.csv("../00_Rawdata_tcga/Dat_Training_count.csv",
                      check.names = F,
                      row.names = 1)
kegggmt <-
  read.gmt("./c2.cp.kegg.v2023.1.Hs.symbols.gmt")

risk <- read.csv('../06_model_validation/risk.csv',
                 check.names = F,
                 row.names = 1)
risk$id <- rownames(risk)
risk <- risk[, c("id", "risk")]
dat_final <- tarin_dat[, risk$id]

set.seed(1)
risk$risk <- factor(risk$risk, levels = c('high risk', 'low risk'))

colData <- data.frame(risk)
dds <-
  DESeqDataSetFromMatrix(
    countData = round(dat_final),
    colData = colData,
    design = ~ risk
  )
dds = dds[rownames(counts(dds)) > 1, ]
dds <- DESeq(dds)
dds

res <- results(dds, contrast = c("risk", 'high risk', 'low risk'))

resOrdered <- res[order(res$padj), ]
DEG <- as.data.frame(resOrdered) %>% na.omit()
saveRDS(DEG, file = 'DEG.rds')
DEG <- readRDS("DEG.rds")

risk_gene_symbol <- rownames(DEG)
risk_gene <- DEG

gene <- bitr(
  risk_gene_symbol,
  fromType = "SYMBOL",
  toType = c("ENTREZID"),
  OrgDb = "org.Hs.eg.db"
)
gene <- dplyr::distinct(gene, SYMBOL, .keep_all = TRUE)
colnames(gene) <- c("GeneSymbol", "ENTREZID")
risk_gene$GeneSymbol <- rownames(risk_gene)
data_all <- risk_gene %>% dplyr::inner_join(., gene, by = "GeneSymbol")
data_all_sort <- data_all %>% arrange(desc(log2FoldChange))
genelist <- data_all_sort$log2FoldChange 
names(genelist) <-
  data_all_sort$GeneSymbol 
write.csv(genelist, file = "genelist.csv")

KEGG <- GSEA(genelist, TERM2GENE = kegggmt, pvalueCutoff = 0.05) 
sortKEGG <- data.frame(KEGG)
sortKEGG <-
  sortKEGG[order(sortKEGG$p.adjust, decreasing = F), ]
kegg_result <- KEGG[KEGG$p.adjust < 0.05 & abs(KEGG$NES) > 1] 
write.csv(kegg_result, file = "kegg_result.csv")
library(ggplot2)
if (T) {
  gsea.plot = function(res.kegg, top.hall) {
    gsdata <-
      do.call(rbind,
              lapply(top.hall, enrichplot:::gsInfo, object = res.kegg))
    gsdata$Description = factor(gsdata$Description, levels = top.hall)
    p1 = ggplot(gsdata, aes_(x = ~ x)) + xlab(NULL) +
      
      theme_classic(14) +
      theme(
        panel.grid.major = element_line(colour = "grey92"),
        panel.grid.minor = element_line(colour = "grey92"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
      ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_color_brewer(palette = "Paired") +
      ggtitle(paste0("Gene Set Enrichment Analysis")) +
      
      geom_line(aes_(y = ~ runningScore, color = ~ Description), size = 1) +
      theme(
        legend.position = "top",
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")
      ) +
      guides(color = guide_legend(
        ncol = 1,
        byrow = TRUE,
        reverse = T
      )) +
      ylab("Running Enrichment Score") +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        text = element_text(face = "bold", family = "Times")
      ) +
      theme(
        axis.title.y = element_text(
          size = 15,
          family = "Times",
          face = "bold"
        ),
        axis.text.y = element_text(
          size = 15,
          family = "Times",
          face = "bold"
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    i = 0
    for (term in unique(gsdata$Description)) {
      idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
      gsdata[idx, "ymin"] <- i
      gsdata[idx, "ymax"] <- i + 1
      i <- i + 1
    }
    p2 = ggplot(gsdata, aes_(x = ~ x)) +
      geom_linerange(aes_(
        ymin = ~ ymin,
        ymax = ~ ymax,
        color = ~ Description
      )) +
      xlab(NULL) + ylab(NULL) +
      
      theme_classic(14) +
      theme(
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        text = element_text(face = "bold", family = "Times")
      ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_fermenter() +
      
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        legend.text = element_text(size = 14)
      )
    p3 <- gseaplot2(res.kegg, top.hall, subplots = 3) +

      theme(
        axis.title.x = element_text(
          size = 15,
          family = "Times",
          face = "bold"
        ),
        axis.text.x = element_text(
          angle = 0,
          size = 15,
          hjust = 0,
          family = "Times",
          face = "bold"
        ),
        axis.title.y = element_text(
          size = 15,
          family = "Times",
          face = "bold"
        ),
        axis.text.y = element_text(
          size = 15,
          family = "Times",
          face = "bold"
        )
      ) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 15) ,
        legend.text = element_text(size = 14)
      )
    p <- p1 / p2 / p3 + plot_layout(ncol = 1, height = c(4, 1, 2))
    return(p)
  }
}
library(patchwork)
library(ggplot2)
library(enrichplot)
paths <- rownames(sortKEGG[c(1:5), ])
pal <- gsea.plot(KEGG, paths)
pal
png(
  "KEGG.png",
  w = 10,
  h = 8,
  units = "in",
  res = 300
)
print(pal)
dev.off()
pdf("KEGG.pdf", width = 10, height = 8)
print(pal)
dev.off()

save.image("../ws_10_GSEA.RData")
