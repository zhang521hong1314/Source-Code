rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("13_checkpoint")) {
  dir.create("13_checkpoint")
}
setwd("13_checkpoint")

dat_1 <- read.csv(
    "../00_Rawdata_tcga/Dat_Training_fpkm.csv",
    row.names = 1,
    check.names = F)

risk <-
  read.csv("../06_model_validation/risk.csv",
           row.names = 1,
           check.names = F)

risk$sample <- rownames(risk)
group <- risk[,c("sample", "risk")]

table(group$risk)
High.sample <- group$sample[which(group$risk == 'high risk')]
Low.sample <- group$sample[which(group$risk == 'low risk')]

check_point <- read.csv("./checkpoint.csv", header = F )

dat2 <- dat_1[check_point$V1, group$sample]

hub_exp2 <- na.omit(dat2)

hub_exp2$Symbol <- rownames(hub_exp2)

library(tidyr)
hub_exp2 <- gather(hub_exp2,
                   key = sample,
                   value = expr,
                   -c("Symbol"))

control.sample <- group$sample[which(group$risk == 'high risk')]
hub_exp2$Group <-
  ifelse(hub_exp2$sample %in% control.sample, 'high risk', 'low risk')

library(rstatix)

stat.test <- hub_exp2 %>%
  group_by(Symbol) %>%
  wilcox_test(expr ~ Group) %>%
  adjust_pvalue(method = 'fdr')

stat.test$p1 <- ifelse(stat.test$p < 0.0001,
                       "****",
                       ifelse(
                         stat.test$p < 0.001,
                         "***",
                         ifelse(stat.test$p < 0.01, "**",
                                ifelse(stat.test$p < 0.05, "*", 'ns'))
                       ))

hub_exp2$Group <-
  factor(hub_exp2$Group, levels = c('high risk', 'low risk'))

library(ggpubr)

hub_exp2$expr <- as.numeric(hub_exp2$expr)
hub_exp2$Group <-
  factor(hub_exp2$Group, levels =  c('high risk', 'low risk'))
x1 = 1:1 
if (sample(x1, 1) == 1) {
  p <- ggboxplot(
    hub_exp2,
    x = "Symbol",
    y = "expr",
    color = "Group",
    palette = c("seagreen4", "darkorange3"),
    outlier.shape = NA,
    bxp.errorbar = T
  ) +
    stat_compare_means(aes(group = Group), label = "p.signif") +
    labs(x = "", y = "Expression", color = "") +
    theme(
      axis.title = element_text(
        size = 20,
        face = "bold",
        family = "Times"
      ),
      axis.text.x = element_text(
        size = 15,
        face = "bold",
        family = "Times",
        angle = 45,
        vjust = 1,
        hjust = 1
      ),
      axis.text.y = element_text(
        size = 15,
        face = "bold",
        family = "Times"
      ),
      legend.text = element_text(size = 15, family = "Times"),
      text = element_text(family = "Times")
    )
} else{
  p <-  ggplot(hub_exp2,
               aes(
                 x = Symbol,
                 y = expr,
                 fill = Group,
                 outlier.shape = NA,
                 bxp.errorbar = T
               )) +
    geom_boxplot(width = 0.8,
                 alpha = 1,
                 position = position_dodge(0.9)) +
    
    
    stat_compare_means(aes(group = Group), label = "p.signif") +
    theme_bw() +
    theme(legend.position = "top") +
    theme(
      axis.title.x = element_text(
        size = 16,
        family = "Times",
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 45,
        size = 16,
        hjust = 1,
        family = "Times",
        face = "bold"
      ),
      axis.title.y = element_text(
        size = 20,
        family = "Times",
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 16,
        family = "Times",
        face = "bold"
      )
    ) +
    theme(legend.title = element_text(size = 15) ,
          legend.text = element_text(size = 14)) +
    labs(x = "", y = "Expression") + ggsci::scale_fill_d3()
}
ggsave(filename = "boxplot.pdf",
       width = 14,
       height = 10,
       p)
ggsave(filename = "boxplot.png",
       width = 14,
       height = 10,
       p)

diff_check_point <- subset(stat.test, p < 0.05)
dat3 <- dat_1[diff_check_point$Symbol, group$sample]
dat3 <- t(dat3) %>% as.data.frame()
riskscore <- risk[, c("riskScore"), drop = FALSE]

library(tibble)
library(magrittr)
library(corrplot)
library(reshape2)
cor_dat <- cbind(riskscore, dat3)
if (T) {
  check_exp <- cor_dat[, -1, drop = F]
  risk_score <- cor_dat[, 1, drop = F]
  library(psych)
  
  d <-
    corr.test(check_exp, risk_score, use = "complete", method = 'spearman')
  r <- data.frame(d$r)
  p <- data.frame(d$p)
  correlation <- data.frame(rownames(p), r$riskScore, p$riskScore)
  colnames(correlation) <- c("Cells", "Correlation", "p.value")
  correlation <- correlation[correlation$p.value < 0.05, ]
  correlation$sig[correlation$p.value > 0.05] <- "ns"
  correlation$sig[correlation$p.value < 0.05 &
                    correlation$p.value > 0.01] <- "*"
  correlation$sig[correlation$p.value < 0.01 &
                    correlation$p.value > 0.001] <- "**"
  correlation$sig[correlation$p.value < 0.001 &
                    correlation$p.value > 0.0001] <- "***"
  correlation$sig[correlation$p.value < 0.0001] <- "****"
  correlation$'correlation' <- abs(correlation$Correlation)
  write.table(
    correlation,
    paste0("01.riskScore_Corr_immuCheckPoint.txt"),
    sep = "\t",
    row.names = F
  )
  library(viridis)
  library(ggpubr)
  
  pdf(
    file = paste0("01.riskScore_Corr_check_point.pdf"),
    width = 10,
    height = 10
  )
  a <- dev.cur()
  png(
    file = paste0("01.riskScore_Corr_check_point.png"),
    width = 10,
    height = 10,
    units = "in",
    res = 300
  )
  dev.control("enable")
  par(mar = c(0, 2, 2, 2),
      cex = 1.5,
      family = "Times")
  
  p0 <- ggdotchart(
    correlation,
    x = "Cells",
    y = "Correlation",
    dot.size = "correlation",
    
    color = "p.value",
    label = 'sig',
    y.label = 0.7,
    sorting = "descending",
    add = "segments",
    rotate = TRUE,
    ggtheme = theme_pubr(margin = F),
    
    ylab = "",
    xlab = '',
    title = paste0("RiskScore Correlation with ImmuneCheckPoint")
  ) +
    scale_colour_gradient(high = "#D1E5F0", low = "#00A087FF") +
    theme(legend.position = "right", panel.background = element_blank()) +
    guides(dot.size = guide_legend(order = 2), color = guide_legend(order = 1)) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", lwd = 0.2) +
    theme(
      axis.title.x = element_text(
        size = 20,
        family = "Times",
        face = "bold"
      ),
      axis.text.x = element_text(
        size = 16,
        family = "Times",
        face = "bold"
      ),
      axis.title.y = element_text(
        size = 20,
        family = "Times",
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 16,
        family = "Times",
        face = "bold"
      ),
      plot.title = element_text(
        size = 20,
        family = "Times",
        face = "bold",
        hjust = 0.5
      ),
      legend.text = element_text(size = 16, family = "Times"),
      legend.title = element_text(
        size = 20,
        family = "Times",
        face = "bold"
      )
    ) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(p0)
  dev.copy(which = a)
  dev.off()
  dev.off()
}
save.image("temp1.RData")

library(psych)
library(corrplot)
library(dplyr)

check <- dat3
prog <- select(risk, -c("OS", "OS.time", "risk", "sample"))

cor_result <- corr.test(check, prog, method = "spearman")  
cor_matrix <- cor_result$r
p_matrix <- cor_result$p

write.csv(cor_matrix, file = 'gene_cor.csv')
write.csv(p_matrix, file = 'gene_pvalue.csv')

library(ggplot2)
library(reshape2)
library(stringr)

cor_data <- melt(cor_matrix, na.rm = TRUE)
colnames(cor_data) <- c("Pathway", "Gene", "Correlation")

p_data <- melt(p_matrix, na.rm = TRUE)
colnames(p_data) <- c("Pathway", "Gene", "pvalue")
p_data$stars <- ifelse(p_data$pvalue < 0.0001, "****",
                       ifelse(p_data$pvalue < 0.001, "***",
                              ifelse(p_data$pvalue < 0.01, "**",
                                     ifelse(p_data$pvalue < 0.05, "*", "ns"))))

cor_data <- merge(cor_data, p_data, by = c("Pathway", "Gene"))

p <- ggplot(cor_data, aes(x = Gene, y = Pathway, fill = Correlation)) +
  geom_tile(color = "white") +  
  scale_fill_gradient2(low = "#85b4d4", mid = "white", high = "#ea5a5b", midpoint = 0, limit = c(-1, 1)) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 10),  
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),  
    axis.title.x = element_text(size = 20, face = "bold"),  
    axis.title.y = element_text(size = 20, face = "bold"),  
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12)  
  ) +
  labs(title = "Significant Correlation", x = "", y = "Checkpoints") +
  guides(
    fill = guide_colorbar(title = "Correlation", barwidth = 1, barheight = 15)
  ) +
  
  
  geom_text(aes(label = round(Correlation, 2)), vjust = 1.5, size = 3) +
  
  geom_text(aes(label = stars), vjust = -0.5, size = 4, color = "black")

ggsave("significant_correlation.pdf", p, width = 8, height = 14)

ggsave("significant_correlation.png", p, width = 8, height = 14, dpi = 300)

save.image("../ws_13_checkpoint_1.RData")

rm(list = ls())
setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("13_checkpoint")) {
  dir.create("13_checkpoint")
}
setwd("13_checkpoint")

library(tidyverse)
library(lance)

dat_1 <- read.csv(
  "../00_Rawdata_tcga/Dat_Training_fpkm.csv",
  row.names = 1,
  check.names = F)

risk <-
  read.csv("../06_model_validation/risk.csv",
           row.names = 1,
           check.names = F)

risk$sample <- rownames(risk)
group <- risk[,c("sample", "risk", "riskScore")]

colnames(group) <- c('sample', 'label', 'riskScore')

Low.sample <- group$sample[which(group$label == 'low risk')]
High.sample <- group$sample[which(group$label == 'high risk')]

dat <- dat_1[, group$sample]

group_estimate <- group$label %>% as.factor()
design <- model.matrix( ~ 0 + group_estimate)
rownames(design) <- group$sample
colnames(design) <- levels(group_estimate)
design <- as.data.frame(design)
Low <- rownames(design)[which(design$`low risk` == 1)]
High <- rownames(design)[which(design$`high risk` == 1)]
length(Low)
length(High)

library(estimate)

expr_train <- dat
write.table(
  expr_train,
  'expr.txt',
  col.names = T,
  row.names = T,
  quote = F,
  sep = "\t"
)

filterCommonGenes(input.f = './expr.txt',
                  output.f = 'expr_train.gct',
                  id = 'GeneSymbol')

estimateScore('expr_train.gct', 'train_purity.gct', platform = "illumina")
es_score <-
  read.table(
    'train_purity.gct',
    skip = 2,
    header = T,
    check.names = F
  )
immu_score <- es_score[, 3:length(es_score)]
rownames(immu_score) <- es_score$NAME
write.table(
  es_score,
  file = "es_score.csv",
  sep = "\t",
  quote = F,
  row.names = F
)

violin_dat <- data.frame(t(immu_score))
rownames(violin_dat) <- group$sample
violin_dat$sample <- rownames(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% Low.sample,
                           "Low risk", "High risk")
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
{
  p1 <- ggplot(violin_dat, aes(x = group, y = StromalScore, fill = group)) +
    
    stat_boxplot(geom = "errorbar",
                 width = 0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width = 0.5,
                 position = position_dodge(0.6),
                 outlier.shape = NA) + 
    
    
    
    scale_fill_manual(values = c("#dc143c", "#6495ed"), name = "Group") +
    stat_compare_means(aes(group = group), method = 'wilcox') +
    theme_bw() + 
    ylim(-2500, 2500) +
    theme(
      panel.background = element_rect(
        colour = 'black',
        size = 1,
        fill = 'white'
      ),
      panel.grid = element_blank()
    ) +
    theme(
      axis.text.x = element_text(
        family = "Times",
        size = 12,
        face = "bold"
      ),
      
      axis.text.y = element_text(
        family = "Times",
        size = 12,
        face = "bold"
      ),
      
      axis.title.y = element_text(
        family = "Times",
        size = 15,
        face = "bold"
      ),
      
      axis.title.x = element_text(
        family = "Times",
        size = 10,
        face = "bold"
      ),
      plot.title = element_text(
        hjust = 0.5,
        family = "Times",
        size = 15,
        face = "bold"
      ),
      legend.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    ) +
    labs(title = "Stromal Score", x = "", y = "Stromal Score") +
    guides(fill = 'none')
  p1
  p2 <- ggplot(violin_dat, aes(x = group, y = ImmuneScore, fill = group)) +
    
    stat_boxplot(geom = "errorbar",
                 width = 0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width = 0.5,
                 position = position_dodge(0.6),
                 outlier.shape = NA) + 
    
    
    
    scale_fill_manual(values = c("#dc143c", "#6495ed"), name = "Group") +
    stat_compare_means(aes(group = group), method = 'wilcox') +
    theme_bw() + 
    ylim(-2000, 4000) +
    theme(
      panel.background = element_rect(
        colour = 'black',
        size = 1,
        fill = 'white'
      ),
      panel.grid = element_blank()
    ) +
    theme(
      axis.text.x = element_text(
        family = "Times",
        size = 12,
        face = "bold"
      ),
      
      axis.text.y = element_text(
        family = "Times",
        size = 12,
        face = "bold"
      ),
      
      axis.title.y = element_text(
        family = "Times",
        size = 15,
        face = "bold"
      ),
      
      axis.title.x = element_text(
        family = "Times",
        size = 10,
        face = "bold"
      ),
      plot.title = element_text(
        hjust = 0.5,
        family = "Times",
        size = 15,
        face = "bold"
      ),
      legend.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    ) +
    labs(title = "Immune Score", x = "", y = "Immune Score") +
    guides(fill = 'none')
  p2
  p3 <-
    ggplot(violin_dat, aes(x = group, y = ESTIMATEScore, fill = group)) +
    
    stat_boxplot(geom = "errorbar",
                 width = 0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width = 0.5,
                 position = position_dodge(0.5),
                 outlier.shape = NA) + 
    
    
    
    scale_fill_manual(values = c("#dc143c", "#6495ed"), name = "Group") +
    stat_compare_means(aes(group = group), method = 'wilcox') +
    theme_bw() + 
    theme(
      panel.background = element_rect(
        colour = 'black',
        size = 1,
        fill = 'white'
      ),
      panel.grid = element_blank()
    ) +
    ylim(-3500, 6000) +
    theme(
      axis.text.x = element_text(
        family = "Times",
        size = 12,
        face = "bold"
      ),
      
      axis.text.y = element_text(
        family = "Times",
        size = 12,
        face = "bold"
      ),
      
      axis.title.y = element_text(
        family = "Times",
        size = 15,
        face = "bold"
      ),
      
      axis.title.x = element_text(
        family = "Times",
        size = 10,
        face = "bold"
      ),
      plot.title = element_text(
        hjust = 0.5,
        family = "Times",
        size = 15,
        face = "bold"
      ),
      legend.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    ) +
    labs(title = "ESTIMATE Score", x = "", y = "ESTIMATE Score") +
    guides(fill = 'none')
  p3
  p5 <- cowplot::plot_grid(p1,
                           p2,
                           p3,
                           nrow = 1,
                           align = 'h',
                           vjust = -0.3)
  p5
}
ggsave(filename = '01.estimate.pdf', p5, w = 8, h = 4)
ggsave(filename = '01.estimate.png', p5, w = 8, h = 4)

library(ggcorrplot)
library(corrplot)
violin_dat1 <- merge(violin_dat, group, by = "sample")
violin_dat2 <- violin_dat1[, c(2:4, 7)]
violin_dat2 <- t(violin_dat2)

cor.dat <- violin_dat2
corr <- round(cor(t(cor.dat), method = 'spearman'), 3)
write.csv(corr, file = 'correlation_r.csv', row.names = T)

corp <- round(cor_pmat(t(cor.dat)), 3)
write.csv(corp, file = 'correlation_p.csv', row.names = T)

col1 <- colorRampPalette(c("#2d649e", 'white', "#dc636c"))

pdf(
  file = 'correlation-1.pdf',
  w = 6,
  h = 6,
  family = 'Times'
)

corr_plot <- corrplot(
  corr,
  method = "circle",
  
  type = "lower",
  
  insig = "blank",
  outline = "white",
  addCoef.col = "black",
  col = col1(200),
  tl.col = 'black',
  tl.offset = 0.4,
  number.font = 2,
  cl.cex = 0.8,
  number.cex = 1.2
)
dev.off()
png(
  file = 'correlation-1.png',
  w = 6,
  h = 6,
  units = 'in',
  res = 600,
  family = 'Times'
)
corr_plot <- corrplot(
  corr,
  method = "circle",
  is.corr = T,
  type = "lower",
  
  insig = "blank",
  outline = "white",
  addCoef.col = "black",
  col = col1(200),
  tl.col = 'black',
  tl.offset = 0.4,
  number.font = 2,
  cl.cex = 0.8,
  number.cex = 1.2
)
dev.off()

save.image("temp2.RData")
save.image("../ws_13_checkpoint_2.RData")

