rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("./15_subcell")) {
  dir.create("./15_subcell")
}
setwd("./15_subcell")

library(tidyverse)
library(IOBR)
loc <-
  read.csv('./tableExport.csv') %>% dplyr::select(-c('X', 'ID', 'mRNALocater'))
colnames(loc)
data <- gather(loc, key = mRNALocater, value = score, -c('FastaID'))
colnames(data)
ggplot(data, aes(x = FastaID, y = score, fill = mRNALocater)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  geom_text(aes(label = score),
            position = position_stack(vjust = 0.5),
            size = 0) +
  scale_fill_manual(
    values = palettes(category = "random", 22, show_col = F),
    labels = c(
      "Cytoplasm",
      "Endoplasmic reticulum",
      "Extracellular region",
      "Mitochondria",
      "Nucleus"
    )
  ) +
  labs(y = "Percentage", title = "") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.key.size = unit(0.5, 'cm'),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 20),
    title = element_text(size = 20),
    legend.text = element_text(size = 15)
  )
ggsave("02.mRNAlocator.png", width = 8, height = 6)
ggsave("02.mRNAlocator.pdf", width = 8, height = 6)

save.image("../ws_15_subcell.RData")
