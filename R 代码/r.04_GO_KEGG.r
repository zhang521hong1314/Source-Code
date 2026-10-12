rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists('./04_GO_KEGG')) {
  dir.create('./04_GO_KEGG')
}
setwd('./04_GO_KEGG')

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggplot2)

gene1 <- read.csv('../03_candidate/CGs.csv',
                  row.names = 1,
                  check.names = F)

  gene2 <- gene1$CGs 

GO_database <- 'org.Hs.eg.db'  

gene <-
  bitr(gene2,
       fromType = 'SYMBOL',
       toType = 'ENTREZID',
       OrgDb = GO_database)

GO <- enrichGO(
  gene = gene$ENTREZID,
  
  OrgDb = GO_database,
  
  keyType = "ENTREZID",
  
  ont = "ALL",
  
  pvalueCutoff = 0.5,
  
  qvalueCutoff = 0.5,
  
  readable = T
)
go_res1 <- data.frame(GO) 

enrichment_fold <- apply(go_res1,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(GeneRatio/BgRatio,2)
})
enrichment_fold <- as.data.frame(enrichment_fold)

go_res2 <- cbind(go_res1, enrichment_fold)

go_res2 <- go_res2[order(go_res2$pvalue, decreasing = FALSE), ]

go_res <- go_res2[go_res2$pvalue < 0.05,]  

write.csv(go_res, file = "./GO_res.csv")

top10_go_res <- head(go_res, 10)

write.csv(top10_go_res, file = "./top10_go_res.csv")

top10_go_res1 <- top10_go_res %>%
  mutate(Term = factor(Description, levels = rev(Description)))

circle_plot <- ggplot(top10_go_res1, aes(x = Term, y = enrichment_fold)) +  
  geom_point(aes(size = Count, color = pvalue)) +  
  scale_color_gradient(low = "#8B0000", high = "#1E3A5F") + 
  coord_flip() +  
  theme_minimal() +
  labs(
    title = "The Most Enriched GO Terms",
    x = "GO Terms",
    y = "Fold Enrichment",
    size = "Count",  
    color = "pvalue"  
  ) +
  theme(
    axis.text = element_text(size = 15, family = "Times"),
    axis.title = element_text(size = 25, family = "Times"),
    plot.title = element_text(hjust = 0.5, size = 26, family = "Times"),
    plot.background = element_rect(fill = "white", color = "white")  
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))  

ggsave("GO_circle_plot.png", circle_plot, width = 10, height = 8, dpi = 600)
ggsave("GO_circle_plot.pdf", circle_plot, width = 10, height = 8)

bar_plot <- ggplot(top10_go_res1, aes(x = Term, y = Count, fill = pvalue)) +  
  geom_bar(stat = "identity") +  
  scale_fill_gradient(low = "#8B0000", high = "#1E3A5F") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "The Most Enriched GO Terms",
    x = "GO Terms",
    y = "Gene number",
    fill = "pvalue"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +  
  
  theme(
    axis.text = element_text(size = 12, family = "Times"),
    axis.title = element_text(size = 14, family = "Times"),
    plot.title = element_text(hjust = 0.5, size = 16, family = "Times"),
    plot.background = element_rect(fill = "white", color = "white")  
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))  

ggsave("GO_bar_plot.png", bar_plot, width = 10, height = 8,dpi = 600)
ggsave("GO_bar_plot.pdf", bar_plot, width = 10, height = 8)

save.image("../ws_04_GO_KEGG_1.RData")

rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3/04_GO_KEGG")

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

gene1 <- read.csv('../03_candidate/CGs.csv',
                  row.names = 1,
                  check.names = F)

gene2 <- gene1$ 

GO_database <- 'org.Hs.eg.db'  
KEGG_database <- 'hsa' 

gene <-
  bitr(gene2,
       fromType = 'SYMBOL',
       toType = 'ENTREZID',
       OrgDb = GO_database)

KEGG <- enrichKEGG(
  gene = gene$ENTREZID,
  keyType = "kegg",
  organism = KEGG_database,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5
)

kegg_res1 <-
  setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_res2 <- data.frame(kegg_res1)

enrichment_fold <- apply(kegg_res2,1,function(x){
  GeneRatio = eval(parse(text = x["GeneRatio"]))
  BgRatio = eval(parse(text = x["BgRatio"]))
  enrichment_fold = round(GeneRatio/BgRatio,2)
})
enrichment_fold <- as.data.frame(enrichment_fold)

kegg_res3 <- cbind(kegg_res2, enrichment_fold)

kegg_res3 <- kegg_res3[order(kegg_res3$pvalue, decreasing = FALSE), ]

kegg_res <- kegg_res3[kegg_res3$pvalue < 0.05, ]  

write.csv(kegg_res , file = "./KEGG_res.csv")

top5_kegg_res <- head(kegg_res, 5)

write.csv(top5_kegg_res, file = "./top5_kegg_res.csv")

top5_kegg_res1 <- top5_kegg_res %>%
  mutate(Term = factor(Description, levels = rev(Description)))

circle_plot <- ggplot(top5_kegg_res1, aes(x = Term, y = enrichment_fold)) +
  geom_point(aes(size = Count, color = pvalue)) +
  scale_color_gradient(low = "#8B0000", high = "#1E3A5F") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "KEGG Pathways Enrichment",
    x = "Pathways",
    y = "Fold Enrichment",
    size = "Count",
    color = "pvalue"
  ) +
  theme(
    axis.text = element_text(size = 15, family = "Times"),
    axis.title = element_text(size = 25, family = "Times"),
    plot.title = element_text(hjust = 0.5, size = 26, family = "Times"),
    plot.background = element_rect(fill = "white", color = "white")  
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))  

ggsave("KEGG_circle_plot.png", circle_plot, width = 10, height = 8, dpi = 600)
ggsave("KEGG_circle_plot.pdf", circle_plot, width = 10, height = 8)

bar_plot <- ggplot(top5_kegg_res1, aes(x = Term, y = Count, fill = pvalue)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#8B0000", high = "#1E3A5F") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "KEGG Pathways Enrichment",
    x = "Pathways",
    y = "Gene number",
    fill = "pvalue"
  ) +
  theme(
    axis.text = element_text(size = 12, family = "Times"),
    axis.title = element_text(size = 14, family = "Times"),
    plot.title = element_text(hjust = 0.5, size = 16, family = "Times"),
    plot.background = element_rect(fill = "white", color = "white")  
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))  

ggsave("KEGG_bar_plot.png", bar_plot, width = 10, height = 8,dpi = 600)
ggsave("KEGG_bar_plot.pdf", bar_plot, width = 10, height = 8)

save.image("../ws_04_GO_KEGG_2.RData")
