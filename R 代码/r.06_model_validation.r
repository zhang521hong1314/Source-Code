rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("./06_model_validation")) {
  dir.create("./06_model_validation")
}
setwd("./06_model_validation")

library(survival)
library(glmnet)
library(dplyr)
library(forestplot)

dat_raw <-
  read.csv(
    "../00_Rawdata_tcga/Dat_Training_fpkm.csv",
    check.names = F,
    row.names = 1
  )

survival <-
  read.csv(
    "../00_Rawdata_tcga/Survival_Training.csv",
    check.names = F,
    row.names = 1
  )

tumor_sample <- read.csv("../00_Rawdata_tcga/Tumor_Sample.csv",
                         check.names = F,
                         row.names = 1)

hub_genes <- read.csv("../05_ppi/ppi.csv",
                      check.names = F,
                      row.names = 1)
colnames(hub_genes) <- "genes"

dat_t <- t(dat_raw)
dat_hub <- dat_t[, colnames(dat_t) %in% hub_genes$genes] %>% as.data.frame()

dat <- merge(survival, dat_hub, by.x = "sample", by.y = "row.names") %>%
  tibble::column_to_rownames(var = 'sample')

dat <- dat[tumor_sample$sample, ]

cox_results <- c()
cox_coefficients <- list()
cox_conf.int <- list()

for (gene in colnames(dat[,3:ncol(dat)])) {

  cox_model <- coxph(Surv(OS.time, OS) ~ dat[,gene], data = dat)
  summary_model <- summary(cox_model)


  if (summary_model$coefficients[5] < 0.05 &
      summary_model$coefficients[2] != 1) {
    cox_coefficients[[gene]] <- summary_model$coefficients[1, ]
    cox_conf.int[[gene]] <- summary_model$conf.int[1, ]
    cox_results <- append(cox_results, gene)
  }
}

print(cox_results)
print(cox_coefficients)
print(cox_conf.int)

table1 <- data.frame(
  Gene = cox_results,
  HR = sapply(cox_results, function(x) cox_coefficients[[x]][2]),
  pvalue = sapply(cox_results, function(x) cox_coefficients[[x]][5]),
  HR.95L = sapply(cox_results, function(x) cox_conf.int[[x]][3]),
  HR.95H = sapply(cox_results, function(x) cox_conf.int[[x]][4])
)
rownames(table1) <- seq_len(nrow(table1))

ph_test_results <- data.frame()

for (gene in cox_results) {
  cox_model <- coxph(Surv(OS.time, OS) ~ dat[,gene], data = dat)
  ph_test <- cox.zph(cox_model)
  ph_test_results = rbind(ph_test_results,
                          cbind(Gene = gene,
                                p = ph_test$table[1, "p"]))
}
sigPH = ph_test_results[as.numeric(as.vector(ph_test_results$p)) > 0.05, ]

table1 <- table1[table1$Gene %in% sigPH$Gene, ]
sigPH <- sigPH[table1$Gene %in% sigPH$Gene, ]
rownames(table1) <- seq_len(nrow(table1))
rownames(sigPH) <- seq_len(nrow(sigPH))

write.csv(table1, "Cox_Results.csv")
write.csv(sigPH, "PH.Sig.txt", row.names = F )

variable <- table1$Gene
HR <- sprintf("%.3f", as.numeric(table1$HR))
HRLow <-
  sprintf("%.3f", as.numeric(table1$HR.95L))
HRHigh <-
  sprintf("%.3f", as.numeric(table1$HR.95H))
HR95 <-
  paste0(HR, "(", HRLow, "-", HRHigh, ")")
pValue <-
  ifelse(table1$pvalue < 0.001,
         "<0.001",
         sprintf("%.3f", as.numeric(table1$pvalue)))

pdf(
  file = "forest.pdf",
  height = 8,
  width = 8,
  family = 'Times'
)

n <- nrow(table1)
nRow <- n + 1
ylim <- c(1, nRow)

layout(matrix(c(1, 2), ncol = 2), widths = c(2, 2))

xlim = c(0, 2.5)
par(mar = c(4, 2.5, 2, 1))
plot(
  0,
  xlim = xlim,
  ylim = ylim,
  type = "n",
  axes = F,
  xlab = "",
  ylab = ""
)

text.cex = 1
text(0, n:1, variable, adj = 0, cex = text.cex)
text(1.1, n:1, pValue, adj = 1, cex = text.cex)
text(1.1,
     n + 1,
     'pValue',
     cex = text.cex,
     font = 2,
     adj = 1)
text(2.5, n:1, HR95, adj = 1, cex = text.cex)
text(2.5,
     n + 1,
     'HR(95% CI)',
     cex = text.cex,
     font = 2,
     adj = 1,
)

abline(h = nRow + 0.15, col = "black", lwd = 2)
abline(h = nRow - 0.15, col = "black", lwd = 2)
abline(h = 0.9, col = "black", lwd = 2)

par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
xlim = c(0, max(c(as.numeric(HRLow), as.numeric(HRHigh))) + 1)
plot(
  0,
  xlim = xlim,
  ylim = ylim,
  type = "n",
  axes = F,
  ylab = "",
  xaxs = "i",
  xlab = "Hazard Ratio"
)
arrows(
  as.numeric(HRLow),
  n:1,
  as.numeric(HRHigh),
  n:1,
  angle = 90,
  code = 3,
  length = 0.05,
  col = "darkblue",
  lwd = 2.5
)
abline(
  v = 1,
  col = "black",
  lty = 2,
  lwd = 2
)
boxcolor = ifelse(as.numeric(HR) > 1, 'red', 'green')
points(as.numeric(HR),
       n:1,
       pch = 15,
       col = boxcolor,
       cex = 1.3)
axis(1)
dev.off()

png(
  filename = "forest.png",
  height = 8,
  width = 8,
  units = 'in',
  res = 600,
  family = 'Times'
)

n <- nrow(table1)
nRow <- n + 1
ylim <- c(1, nRow)

layout(matrix(c(1, 2), ncol = 2), widths = c(2, 2))

xlim = c(0, 2.5)
par(mar = c(4, 2.5, 2, 1))
plot(
  0,
  xlim = xlim,
  ylim = ylim,
  type = "n",
  axes = F,
  xlab = "",
  ylab = ""
)
text.cex = 1
text(0, n:1, variable, adj = 0, cex = text.cex)
text(1.1, n:1, pValue, adj = 1, cex = text.cex)
text(1.1,
     n + 1,
     'pValue',
     cex = text.cex,
     font = 2,
     adj = 1)
text(2.5, n:1, HR95, adj = 1, cex = text.cex)
text(2.5,
     n + 1,
     'HR(95% CI)',
     cex = text.cex,
     font = 2,
     adj = 1,
)

abline(h = nRow + 0.15, col = "black", lwd = 2)
abline(h = nRow - 0.15, col = "black", lwd = 2)
abline(h = 0.9, col = "black", lwd = 2)

par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
xlim = c(0, max(c(as.numeric(HRLow), as.numeric(HRHigh))) + 1)
plot(
  0,
  xlim = xlim,
  ylim = ylim,
  type = "n",
  axes = F,
  ylab = "",
  xaxs = "i",
  xlab = "Hazard Ratio"
)
arrows(
  as.numeric(HRLow),
  n:1,
  as.numeric(HRHigh),
  n:1,
  angle = 90,
  code = 3,
  length = 0.05,
  col = "darkblue",
  lwd = 2.5
)
abline(
  v = 1,
  col = "black",
  lty = 2,
  lwd = 2
)
boxcolor = ifelse(as.numeric(HR) > 1, 'red', 'green')
points(as.numeric(HR),
       n:1,
       pch = 15,
       col = boxcolor,
       cex = 1.3)
axis(1)
dev.off()

x_all <- subset(dat, select = -c(OS.time, OS))
x_all <- x_all[, sigPH$Gene]
y_all <- subset(dat, select = c(OS.time, OS))

fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time, y_all$OS),
              family = "cox")

set.seed(18)
cvfit = cv.glmnet(
  as.matrix(x_all),
  Surv(y_all$OS.time, y_all$OS),
  nfold = 10,
  family = "cox"
)

coef.min = coef(cvfit, s = "lambda.min")
coef_min_cal <- as.data.frame(as.matrix(coef.min))
colnames(coef_min_cal) <- 'coef'
coef_min_cal <- t(coef_min_cal) %>% as.data.frame()

cvfit$lambda.min

pdf('lasso.CV.pdf', width = 6, height = 5, family = 'Times')

plot(cvfit)

abline(v = log(cvfit$lambda.min), col = "blue", lwd = 2, lty = 2)

text(
  x = log(cvfit$lambda.min),
  y = max(cvfit$cvm),
  labels = "lambda.min",
  col = "blue",
  cex = 0.8,
  pos = 2
)

dev.off()

png('lasso.CV.png', width = 6, height = 5, units = 'in', res = 300)

plot(cvfit)

abline(v = log(cvfit$lambda.min), col = "blue", lwd = 2, lty = 2)

text(
  x = log(cvfit$lambda.min),
  y = max(cvfit$cvm),
  labels = "lambda.min",
  col = "blue",
  cex = 0.8,
  pos = 2
)

dev.off()

library(RColorBrewer)

gene_names <- rownames(cvfit$glmnet.fit$beta)

num_genes <- length(gene_names)

gene_colors <- colorRampPalette(brewer.pal(min(num_genes, 9), "Set1"))(num_genes)

pdf('lasso.Coef.pdf', width = 6, height = 5, family = 'Times')
plot(cvfit$glmnet.fit, xvar = 'lambda', label = FALSE, col = "grey")

for (i in seq_along(gene_names)) {
  lines(
    log(cvfit$glmnet.fit$lambda),
    cvfit$glmnet.fit$beta[i, ],
    col = gene_colors[i],
    lwd = 1
  )
}

abline(v = log(cvfit$lambda.min), lty = 3, lwd = 1, col = "blue")

text(
  x = log(cvfit$lambda.min),
  y = 0.2,
  labels = "lambda.min",
  col = "blue",
  cex = 0.8,
  pos = 2
)

legend(
  "topright",
  legend = gene_names,
  col = gene_colors,
  lty = 1,
  lwd = 2,
  cex = 0.5
)

dev.off()

png('lasso.Coef.png', width = 6, height = 5, units = 'in', res = 300)
plot(cvfit$glmnet.fit, xvar = 'lambda', label = FALSE, col = "grey")

for (i in seq_along(gene_names)) {
  lines(
    log(cvfit$glmnet.fit$lambda),
    cvfit$glmnet.fit$beta[i, ],
    col = gene_colors[i],
    lwd = 1
  )
}

abline(v = log(cvfit$lambda.min), lty = 3, lwd = 1, col = "blue")

par(xpd = TRUE)

text(
  x = log(cvfit$lambda.min),
  y = 0.2,
  labels = "lambda.min",
  col = "blue",
  cex = 0.8,
  pos = 2
)

legend(
  "topright",
  legend = gene_names,
  col = gene_colors,
  lty = 1,
  lwd = 2,
  cex = 0.5
)

dev.off()

length(coef.min)

lasso_gene <-
  coef.min@Dimnames[[1]][coef.min@i + 1] %>% as.data.frame()
colnames(lasso_gene) <- 'genes'
lasso_gene

coef_min_cal <- coef_min_cal[, lasso_gene$genes]

write.csv(lasso_gene, "lasso_genes.csv")
write.csv(coef_min_cal, file = 'coef_min_cal.csv')

riskScore_cal <- x_all[, lasso_gene$genes]
coef_min_cal <- coef_min_cal[, colnames(riskScore_cal)]
for (i in 1:nrow(riskScore_cal)) {
  riskScore_cal[i, ] <- riskScore_cal[i, ] * coef_min_cal
}
riskScore_cal$riskScore <- rowSums(riskScore_cal)
riskScore_cal <- merge(riskScore_cal, y_all, by.x = "row.names", by.y = "row.names")
riskScore_cal <- rename(riskScore_cal, "id" = "Row.names")

library(survminer)
dat_score <- select(riskScore_cal, 'riskScore', 'OS.time', 'OS')
res.cut <-
  surv_cutpoint(
    dat_score,
    time = 'OS.time',
    event = 'OS',
    minprop = 0.1,

    variables = 'riskScore'
  )
res.cut
cutpoint <- res.cut$cutpoint[1]
cutpoint$cutpoint

riskScore_cal$risk <- ifelse(riskScore_cal$riskScore > cutpoint$cutpoint, 1, 0)

risk <- riskScore_cal
risk <- na.omit(risk)
riskScore <- risk$riskScore

risk2 <- risk
risk2$risk <- ifelse(risk2$risk == 1, 'high risk', 'low risk')
table(risk2$risk)

write.csv(risk2,
          file = 'risk.csv',
          quote = F,
          row.names = F)

library(ggplot2)
library(ggthemes)
library(Ipaper)

risk_dis <- ggplot(risk, aes(
  x = reorder(id, riskScore),
  y = riskScore,
  color = factor(
    risk,
    levels = c(1, 0),
    labels = c("High Risk", "Low Risk")
  )
)) +
  geom_point() +
  scale_color_manual(values = c("#c76162", "#1d86b5")) +
  scale_x_discrete(
    breaks = risk[order(risk$riskScore), ]$id[c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110)],
    labels = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110),
    expand = c(0.02, 0)
  ) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == 0), ]) + 0.5,
             lty = 2) +
  geom_hline(yintercept = cutpoint$cutpoint,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Train risk and survival Distribution") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.justification = c(1, 0),
    legend.position.inside = c(1, 0),

    legend.background = element_rect(color = "black", linewidth = .3),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title = element_text(family = 'Times'),
    text = element_text(family = 'Times')

  )
risk_dis

surv_stat <- ggplot(risk, aes(
  x = reorder(id, riskScore),
  y = OS.time / 365,
  color = factor(OS,
                 levels = c(0, 1),
                 labels = c("Alive", "Dead"))
)) +
  geom_point() +
  scale_color_manual(values = c("#1d86b5", "#c76162")) +
  scale_x_discrete(
    breaks = risk[order(risk$riskScore), ]$id[c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110)],
    labels = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110),
    expand = c(0.02, 0)
  ) +
  ylim(c(0, 20)) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == 0), ]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Time (Years)",
       title = "") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.justification = c(0, 1),
    legend.position.inside = c(0, 1),

    legend.background = element_rect(color = "black", linewidth = .3),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title = element_text(family = 'Times'),
    text = element_text(family = 'Times')
  )
surv_stat
library(gridExtra)
p <- grid.arrange(risk_dis, surv_stat)

ggsave(filename = 'Train_risk_and_survival_Distribution.pdf', p, w =
         8, h = 8)
ggsave(
  filename = 'Train_risk_and_survival_Distribution.png',
  p,
  w = 8,
  h = 8,
  units = 'in',
  dpi = 600
)

library(survminer)

kmfit <- survfit(Surv(OS.time, OS) ~ risk, data =  risk)
train_survival_median <- ggsurvplot(
  kmfit,
  pval = TRUE,
  conf.int = F,

  legend.labs = c("Low risk", "High risk"),

  legend.title = "Risk score",

  title = "Kaplan-Meier Survival Curve in Training Set",

  font.main = c(15, "bold"),
  risk.table = TRUE,

  risk.table.col = "strata",

  linetype = "strata",

  ggtheme = theme_bw(),
  palette = c("#1d86b5", "#c76162")
)
train_survival_median

pdf(
  file = paste0("Train_KM.pdf"),
  width = 5,
  height = 5,
  family = 'Times'
)
a <- dev.cur()
png(
  file = paste0("Train_KM.png"),
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
print(train_survival_median)
dev.copy(which = a)
dev.off()
dev.off()

multi_ROC <- function(time_vector, risk_score_table) {
  library(survivalROC)
  single_ROC <- function(single_time) {
    for_ROC <- survivalROC(
      Stime = risk_score_table$OS.time,
      status = risk_score_table$OS,
      marker = risk_score_table$riskScore,
      predict.time = single_time,
      method = 'KM'
    )
    data.frame(
      'True_positive' = for_ROC$TP,
      'False_positive' = for_ROC$FP,
      'Cut_values' = for_ROC$cut.values,
      'Time' = rep(single_time, length(for_ROC$TP)),
      'AUC' = rep(for_ROC$AUC, length(for_ROC$TP))
    )
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC(time_vector = c(365 * seq(1, 3, 1)),
                           risk_score_table = risk)
table(for_multi_ROC$AUC)
a <-
  unique(data.frame(Time = for_multi_ROC$Time, AUC = for_multi_ROC$AUC))
a

for_multi_ROC$Time <- factor(for_multi_ROC$Time)

library(scales)
library(geomROC)

library(plotROC)

auc_y1 <-
  round(for_multi_ROC[which(for_multi_ROC$Time == 365), 5][1], 3)
auc_y2 <-
  round(for_multi_ROC[which(for_multi_ROC$Time == 730), 5][1], 3)
auc_y3 <-
  round(for_multi_ROC[which(for_multi_ROC$Time == 1095), 5][1], 3)

auc_y4 <-
  round(for_multi_ROC[which(for_multi_ROC$Time == 1460), 5][1], 3)
auc_y5 <-
  round(for_multi_ROC[which(for_multi_ROC$Time == 1825), 5][1], 3)

ROC <- ggplot(
  for_multi_ROC,
  aes(
    x = False_positive,
    y = True_positive,
    label = Cut_values,
    color = Time
  )
) +
  scale_color_manual(
    breaks = c("365", "730", "1095"),
    labels = c("1 year", "2 year", "3 year"),
    values = c("#c76162", "#51ad51", '#944499')
  ) +
  geom_roc(labels = F, stat = 'identity') +
  style_roc() +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = 'gray',
    linetype = 2
  ) +
  theme_bw() +
  labs(title = "ROC Curve in Training Set") +
  theme(
    panel.background = element_rect(
      colour = 'black',
      linewidth = 1,
      fill = 'white'
    ),
    panel.grid = element_blank(),
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold")
  ) +



  annotate(
    "text",
    x = 0.75,
    y = c(0.25, 0.15, 0.05),
    label = c(
      paste('AUC of 1 year =', format(auc_y1, nsmall = 2)),
      paste('AUC of 2 year =', format(auc_y2, nsmall = 2)),
      paste('AUC of 3 year =', format(auc_y3, nsmall = 2))


    )
  )
ROC

pdf('Train_ROC.pdf',
    width = 6,
    height = 5,
    family = 'Times')
ROC
dev.off()
png(
  'Train_ROC.png',
  w = 6,
  h = 5,
  units = 'in',
  res = 600,
  family = 'Times'
)
ROC
dev.off()

save.image("../ws_06_model_validation_1.RData")

rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3/06_model_validation")

library(dplyr)

dat_raw <- read.csv('../00_Rawdata2/Dat_Training.csv',
                     check.names = F,
                     row.names = 1)

survival <- read.csv(
  "../00_Rawdata2/Survival_Training.csv",
  check.names = F,
  row.names = 1
)

tumor_sample <- read.csv("../00_Rawdata2/Tumor_Sample.csv",
                         check.names = F,
                         row.names = 1)

coef <- read.csv("./coef_min_cal.csv",
                         check.names = F,
                         row.names = 1)

lasso_gene <- read.csv("./lasso_genes.csv",
                       check.names = F,
                       row.names = 1)

dat_t <- as.data.frame(t(dat_raw))

dat <- merge(dat_t, survival, by.x = "row.names", by.y = "sample") %>%
  tibble::column_to_rownames(var = 'Row.names')

dat <- dat[tumor_sample$sample, ]

x_all <- subset(dat, select = -c(OS.time, OS))
x_all <- x_all[, lasso_gene$genes]
y_all <- subset(dat, select = c(OS.time, OS))

riskScore_cal <- x_all
coef <- coef[, colnames(riskScore_cal)]
for (i in 1:nrow(riskScore_cal)) {
  riskScore_cal[i, ] <- riskScore_cal[i, ] * coef
}
riskScore_cal$riskScore <- rowSums(riskScore_cal)
riskScore_cal <- merge(riskScore_cal, y_all, by.x = "row.names", by.y = "row.names")
riskScore_cal <- rename(riskScore_cal, "id" = "Row.names")

library(survminer)
dat_score <- select(riskScore_cal, 'riskScore', 'OS.time', 'OS')
res.cut <-
  surv_cutpoint(
    dat_score,
    time = 'OS.time',
    event = 'OS',
    minprop = 0.1,

    variables = 'riskScore'
  )
res.cut
cutpoint <- res.cut$cutpoint[1]
cutpoint$cutpoint

riskScore_cal$risk <- ifelse(riskScore_cal$riskScore > cutpoint$cutpoint, 1, 0)

risk <- riskScore_cal
risk <- na.omit(risk)
riskScore <- risk$riskScore

risk2 <- risk
risk2$risk <- ifelse(risk2$risk == 1, 'high risk', 'low risk')
table(risk2$risk)

write.csv(risk2,
          file = 'risk_test.csv',
          quote = F,
          row.names = F)

library(ggplot2)
library(ggthemes)
library(Ipaper)

risk_dis <- ggplot(risk, aes(
  x = reorder(id, riskScore),
  y = riskScore,
  color = factor(
    risk,
    levels = c(1, 0),
    labels = c("High Risk", "Low Risk")
  )
)) +
  geom_point() +
  scale_color_manual(values = c("#c76162", "#1d86b5")) +
  scale_x_discrete(
    breaks = risk[order(risk$riskScore), ]$id[c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110)],
    labels = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110),
    expand = c(0.02, 0)
  ) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == 0), ]) + 0.5,
             lty = 2) +
  geom_hline(yintercept = cutpoint$cutpoint,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Test risk and survival Distribution") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.justification = c(1, 0),
    legend.position.inside = c(1, 0),

    legend.background = element_rect(color = "black", linewidth = .3),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title = element_text(family = 'Times'),
    text = element_text(family = 'Times')

  )
risk_dis

surv_stat <- ggplot(risk, aes(
  x = reorder(id, riskScore),
  y = OS.time / 365,
  color = factor(OS,
                 levels = c(0, 1),
                 labels = c("Alive", "Dead"))
)) +
  geom_point() +
  scale_color_manual(values = c("#1d86b5", "#c76162")) +
  scale_x_discrete(
    breaks = risk[order(risk$riskScore), ]$id[c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110)],
    labels = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110),
    expand = c(0.02, 0)
  ) +
  ylim(c(0, 20)) +
  geom_vline(xintercept = nrow(risk[which(risk$risk == 0), ]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Time (Years)",
       title = "") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.justification = c(0, 1),
    legend.position.inside = c(0, 1),

    legend.background = element_rect(color = "black", linewidth = .3),
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.title = element_text(family = 'Times'),
    text = element_text(family = 'Times')
  )
surv_stat
library(gridExtra)
p <- grid.arrange(risk_dis, surv_stat)

ggsave(filename = 'Test_risk_and_survival_Distribution.pdf', p, w =
         8, h = 8)
ggsave(
  filename = 'Test_risk_and_survival_Distribution.png',
  p,
  w = 8,
  h = 8,
  units = 'in',
  dpi = 600
)

library(survminer)
library(survival)

kmfit <- survfit(Surv(OS.time, OS) ~ risk, data =  risk)

train_survival_median <- ggsurvplot(
  kmfit,
  pval = TRUE,
  conf.int = F,

  legend.labs = c("Low risk", "High risk"),

  legend.title = "Risk score",

  title = "Kaplan-Meier Survival Curve in Testing Set",

  font.main = c(15, "bold"),
  risk.table = TRUE,

  risk.table.col = "strata",

  linetype = "strata",

  ggtheme = theme_bw(),
  palette = c("#1d86b5", "#c76162")
)
train_survival_median

pdf(
  file = paste0("Test_KM.pdf"),
  width = 5,
  height = 5,
  family = 'Times'
)
a <- dev.cur()
png(
  file = paste0("Test_KM.png"),
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
print(train_survival_median)
dev.copy(which = a)
dev.off()
dev.off()

multi_ROC <- function(time_vector, risk_score_table) {
  library(survivalROC)
  single_ROC <- function(single_time) {
    for_ROC <- survivalROC(
      Stime = risk_score_table$OS.time,
      status = risk_score_table$OS,
      marker = risk_score_table$riskScore,
      predict.time = single_time,
      method = 'KM'
    )
    data.frame(
      'True_positive' = for_ROC$TP,
      'False_positive' = for_ROC$FP,
      'Cut_values' = for_ROC$cut.values,
      'Time' = rep(single_time, length(for_ROC$TP)),
      'AUC' = rep(for_ROC$AUC, length(for_ROC$TP))
    )
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC(time_vector = c(365 * seq(1, 3, 1)),
                           risk_score_table = risk)
table(for_multi_ROC$AUC)
a <-
  unique(data.frame(Time = for_multi_ROC$Time, AUC = for_multi_ROC$AUC))
a

for_multi_ROC$Time <- factor(for_multi_ROC$Time)

library(scales)
library(geomROC)

library(plotROC)

auc_y1 <-
  round(for_multi_ROC[which(for_multi_ROC$Time == 365), 5][1], 3)
auc_y2 <-
  round(for_multi_ROC[which(for_multi_ROC$Time == 730), 5][1], 3)
auc_y3 <-
  round(for_multi_ROC[which(for_multi_ROC$Time == 1095), 5][1], 3)

auc_y4 <-
  round(for_multi_ROC[which(for_multi_ROC$Time == 1460), 5][1], 3)
auc_y5 <-
  round(for_multi_ROC[which(for_multi_ROC$Time == 1825), 5][1], 3)

ROC <- ggplot(
  for_multi_ROC,
  aes(
    x = False_positive,
    y = True_positive,
    label = Cut_values,
    color = Time
  )
) +
  scale_color_manual(
    breaks = c("365", "730", "1095"),
    labels = c("1 year", "2 year", "3 year"),
    values = c("#c76162", "#51ad51", '#944499')
  ) +
  geom_roc(labels = F, stat = 'identity') +
  style_roc() +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = 'gray',
    linetype = 2
  ) +
  theme_bw() +
  labs(title = "ROC Curve in Testing Set") +
  theme(
    panel.background = element_rect(
      colour = 'black',
      linewidth = 1,
      fill = 'white'
    ),
    panel.grid = element_blank(),
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold")
  ) +



  annotate(
    "text",
    x = 0.75,
    y = c(0.25, 0.15, 0.05),
    label = c(
      paste('AUC of 1 year =', format(auc_y1, nsmall = 2)),
      paste('AUC of 2 year =', format(auc_y2, nsmall = 2)),
      paste('AUC of 3 year =', format(auc_y3, nsmall = 2))


    )
  )
ROC

pdf('Test_ROC.pdf',
    width = 6,
    height = 5,
    family = 'Times')
ROC
dev.off()
png(
  'Test_ROC.png',
  w = 6,
  h = 5,
  units = 'in',
  res = 600,
  family = 'Times'
)
ROC
dev.off()

save.image("../ws_06_model_validation_2.RData")
