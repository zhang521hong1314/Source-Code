.libPaths("/data/nas2/Software/miniconda3/envs/public_R/lib/R/library")

rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("./12_ic50")) {
  dir.create("./12_ic50")
}
setwd("./12_ic50")

library(tidyverse)
library(lance)
dat <-
  read.csv(
    "../00_Rawdata_tcga/Dat_Training_fpkm.csv",
    row.names = 1,
    check.names = F
  )

risk <-
  read.csv("../06_model_validation/risk.csv",
           row.names = 1,
           check.names = F)

dat <- dat[, rownames(risk)]

library(pRRophetic)
library(ggplot2)
set.seed(12345)
drug <- read.table(file = './drugs.txt', sep = '\t', header = F)
a <- data.frame(row.names = rownames(risk), rownames(risk))
cnt <- 1
while (cnt < 139) {
  predictedPtype <-
    pRRopheticPredict(as.matrix(dat), drug[cnt,], selection = 1)
  Tipifarnib <- data.frame(predictedPtype)
  colnames(Tipifarnib) <- drug[cnt,]
  a <- cbind(a, Tipifarnib)
  cnt = cnt + 1
}
write.csv(a, 'IC50.csv')

a <- read.csv("IC50.csv", header = TRUE, row.names = 1)
colnames(a)[1] <- "sample"

risk <- rownames_to_column(risk, var = "sample")
risk <- risk[, c("sample", "risk")]

drug <- merge(risk, a, by = "sample")

durg <- drug
durg <- rename(durg, 'group' = 'risk')
violin_dat <-
  gather(durg, key = sample, value = ic50, -c("sample", "group"))
objects_to_save <- list(
  durg = durg,
  drug = drug,
  a = a,
  data = dat
)

saveRDS(objects_to_save, file = "IC50.RDS")
library(dplyr)
library(broom)
library(rstatix)
library(ggpubr)

stat.test <- violin_dat %>%
  group_by(sample) %>%
  wilcox_test(ic50 ~ group) %>%
  adjust_pvalue(method = 'fdr')
stat.test
DE_drug <- stat.test[stat.test$p.adj < 0.05,]
DE_drug <- DE_drug %>% arrange(p.adj)

write.csv(file = "01.DE_drug.csv", DE_drug)

de <- DE_drug$sample[1:10]
violin_test <- violin_dat[violin_dat$sample %in% de,]
stat.test <- stat.test[stat.test$sample %in% de,]

stat.test <- stat.test %>% arrange(p.adj)
violin_test$sample <- factor(violin_test$sample, levels = stat.test$sample)

stat.test <- stat.test %>%
  group_by(sample) %>%
  mutate(y.position = max(violin_test$ic50[violin_test$sample == sample], na.rm = TRUE) + 0.1) %>% ungroup()

stat.test <- stat.test %>%
  mutate(p.signif = case_when(
    p.adj <= 0.0001 ~ "****",
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01 ~ "**",
    p.adj <= 0.05 ~ "*",
    TRUE ~ "ns"
  ))

table(violin_test$sample)
paste0(unique(violin_test$sample), collapse = "、")

plot <-
  ggplot(violin_test, aes(x = group, y = ic50, color = group)) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.1,
    color = "black",
    position = position_dodge(0.9)
  ) +
  geom_boxplot(
    width = 0.3,
    color = "black",
    aes(fill = group),
    position = position_dodge(0.9),
    outlier.shape = NA
  ) +




  scale_fill_manual(values = c("#FF6A6A", "#20B2AA"), name = "Group") +
  labs(title = "Drug",
       x = "",
       y = "",
       size = 20) +
  stat_pvalue_manual(
    stat.test,
    label = "p.signif",
    y.position = "y.position",
    size = 3.2,
    family = "Times",
    face = "bold"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      colour = "black",
      face = "bold",
      size = 18
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      colour = "black",
      face = "bold",
      size = 10
    ),
    axis.text.y = element_text(
      hjust = 0.5,
      colour = "black",
      face = "bold",
      size = 12
    ),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.text = element_text(
      face = "bold",
      hjust = 0.5,
      colour = "black",
      size = 12
    ),
    legend.title = element_text(face = "bold", size = 12),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(family = 'Times')
  ) +
  facet_wrap(~ sample,
             scales = "free",
             nrow = 2,
             strip.position = "top")
plot
ggsave(filename = '01.ic50.plot.pdf',
       plot,
       w = 12,
       h = 8)
ggsave(
  filename = '01.ic50.plot.png',
  plot,
  w = 12,
  h = 8
)

save.image("../ws_12_ic50.RData")




