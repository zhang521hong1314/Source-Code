rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("./03_candidate")) {
  dir.create("./03_candidate")
}
setwd("./03_candidate")

library(dplyr)
library(ggvenn)

DEGs1 <- read.csv('../01_DEG/DEGs_sig.csv',
                 check.names = F)
colnames(DEGs1)[1] <- 'genes'

DEGs2 <- read.csv('../01_hypotype/DEGs_sig.csv',
                     check.names = F)
colnames(DEGs2)[1] <- 'genes'

CGs <-
  Reduce(intersect,
         list(DEGs1$genes, DEGs2$genes))
CGs <- data.frame(CGs)
write.csv(CGs, file = 'CGs.csv')

mydata <- list(
  'DEGs1' = DEGs1$genes,
  'DEGs2' = DEGs2$genes
)

png(
  'venn_CGs.png',
  w = 8,
  h = 8,
  units = 'in',
  res = 600,
  family = 'Times'
)
ggvenn(
  mydata,
  c('DEGs1', 'DEGs2'),
  fill_color = c("red", "turquoise"),
  show_percentage = T,
  stroke_alpha = 0.5,
  stroke_size = 1,
  text_size = 5,
  stroke_linetype = "dotdash",
  set_name_color = c("red", "turquoise"),
  text_color = 'black'
)
dev.off()

pdf('venn_CGs.pdf',
    width = 8,
    height = 8,
    family = 'Times')
ggvenn(
  mydata,
  c('DEGs1', 'DEGs2'),
  fill_color = c("red", "turquoise"),
  show_percentage = T,
  stroke_alpha = 0.5,
  stroke_size = 1,
  text_size = 5,
  stroke_linetype = "dotdash",
  set_name_color = c("red", "turquoise"),
  text_color = 'black'
)
dev.off()

save.image("../ws_03_candidate.RData")
