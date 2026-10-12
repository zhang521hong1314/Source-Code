rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists("./09_nomogram")) {
  dir.create("./09_nomogram")
}
setwd("./09_nomogram")

library(rms)
library(survival)
library(regplot)

model_expr <- read.csv("../06_model_validation/risk.csv",
                       row.names = 1,
                       check.names = F)
model_expr$risk <- ifelse(model_expr$risk == 'high risk',1,0)

hubgene <-
  read.csv(
    "../06_model_validation/lasso_genes.csv",
    row.names = 1,
    check.names = F
  )

multi_cov <- hubgene$genes
cox_data_prog <- as.formula(paste0('Surv(OS.time, OS)~',
                                   paste(
                                     multi_cov,
                                     sep = '',
                                     collapse = '+'
                                   )))
cox_more_prog <- coxph(cox_data_prog,
                       data = as.data.frame(model_expr))
ncox <- coxph(cox_data_prog, data = as.data.frame(model_expr))

regplot(
  ncox,
  observation = model_expr[10, ],
  
  points = TRUE,
  
  title = "Nomogram",
  clickable = F,  
  subticks = T,
  failtime = c(365 * 3,  365 * 2, 365 * 1),
  prfail = T,  
  showP = F,  
  droplines = T, 
  
  interval = "confidence",
  plots = c("density", "boxes")
)

save.image(file = "testimage.Rdata")

ddist <- datadist(model_expr) 
options(datadist = 'ddist')

coxm_1 <- cph(
  cox_data_prog,
  data = model_expr,
  surv = T,
  x = T,
  y = T,
  time.inc = 1 * 365
)
cal_1 <- calibrate(
  coxm_1,
  u = 1 * 365,
  cmethod = 'KM',
  m = 10,  
  B = 1000  
)

coxm_2 <- cph(
  cox_data_prog,
  data = model_expr,
  surv = T,
  x = T,
  y = T,
  time.inc = 2 * 365
)
cal_2 <- calibrate(
  coxm_2,
  u = 2 * 365,
  cmethod = 'KM',
  m = 10,
  B = 1000
)

coxm_3 <- cph(
  cox_data_prog,
  data = model_expr,
  surv = T,
  x = T,
  y = T,
  time.inc = 3 * 365
)
cal_3 <- calibrate(
  coxm_3,
  u = 3 * 365,
  cmethod = 'KM',
  m = 10,
  B = 1000
)

pdf(file = '01.calibrate.pdf', w = 9, h = 8)
par(mar = c(7, 4, 4, 3), cex = 1.5)
plot(
  cal_1,
  subtitles = F,
  lwd = 2,
  lty = 1,
  
  errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
  
  xlab = 'Nomogram-Predicted Probability of 1,2,3 year Survival Probability',
  
  ylab = 'Actual 1,2,3 year Survival Probability',
  
  col = "#42b540",
  
  xlim = c(0, 1),
  ylim = c(0, 1),
  
  family = "serif",
  font = 2
)

plot(
  cal_2,
  add = T,
  subtitles = F,
  lwd = 2,
  lty = 1,
  
  errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
  
  xlab = 'Nomogram-Predicted Probability of 1,2,3 year Survival Probability',
  
  ylab = 'Actual 1,2,3 year Survival Probability',
  
  col = "#00468b",
  
  xlim = c(0, 1),
  ylim = c(0, 1),
  
  family = "serif",
  font = 2
)
plot(
  cal_3,
  add = T,
  subtitles = F,
  lwd = 2,
  lty = 1,
  
  errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
  
  xlab = 'Nomogram-Predicted Probability of 1,2,3 year Survival Probability',
  
  ylab = 'Actual 1,2,3 year Survival Probability',
  
  col = "#ed0000",
  
  xlim = c(0, 1),
  ylim = c(0, 1),
  
  family = "serif",
  font = 2
)

legend(
  "bottomright",
  legend = c("1-year", "2-year", "3-year"),
  col = c("#42b540", "#00468b", "#ed0000"),
  lwd = 2,
  text.font = 1
)

abline(0,
       1,
       lty = 5,
       lwd = 2,
       col = "grey")
dev.off()

png(
  file = '01.calibrate.png',
  w = 9,
  h = 8,
  units = 'in',
  res = 600
)
plot(
  cal_1,
  subtitles = F,
  lwd = 2,
  lty = 1,
  
  errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
  
  xlab = 'Nomogram-Predicted Probability of 1,2,3 year Survival Probability',
  
  ylab = 'Actual 1,2,3 year Survival Probability',
  
  col = "#42b540",
  
  xlim = c(0, 1),
  ylim = c(0, 1),
  
  family = "serif",
  font = 2
)

plot(
  cal_2,
  add = T,
  subtitles = F,
  lwd = 2,
  lty = 1,
  
  errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
  
  xlab = 'Nomogram-Predicted Probability of 1,2,3 year Survival Probability',
  
  ylab = 'Actual 1,2,3 year Survival Probability',
  
  col = "#00468b",
  
  xlim = c(0, 1),
  ylim = c(0, 1),
  
  family = "serif",
  font = 2
)
plot(
  cal_3,
  add = T,
  subtitles = F,
  lwd = 2,
  lty = 1,
  
  errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
  
  xlab = 'Nomogram-Predicted Probability of 1,2,3 year Survival Probability',
  
  ylab = 'Actual 1,2,3 year Survival Probability',
  
  col = "#ed0000",
  
  xlim = c(0, 1),
  ylim = c(0, 1),
  
  family = "serif",
  font = 2
)

legend(
  "bottomright",
  legend = c("1-year", "2-year", "3-year"),
  col = c("#42b540", "#00468b", "#ed0000"),
  lwd = 2,
  text.font = 1
)

abline(0,
       1,
       lty = 5,
       lwd = 2,
       col = "grey")
dev.off()

library(rmda)

dca.dat <- model_expr
colnames(dca.dat)

nomogram <-
  decision_curve(
    OS ~ C2 + DLGAP1 ,  
    data = dca.dat,
    family = binomial(link = 'logit'),
    thresholds = seq(0, 1, by = 0.01),
    confidence.intervals = 0.95,
    study.design = 'case-control',
    population.prevalence = 0.3
  )
C2.dca <- decision_curve(
  OS ~ C2,
  data = dca.dat,
  family = binomial(link = 'logit'),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = 0.95,
  study.design = 'case-control',
  population.prevalence = 0.3
)
DLGAP1.dca <-
  decision_curve(
    OS ~ DLGAP1,
    data = dca.dat,
    family = binomial(link = 'logit'),
    thresholds = seq(0, 1, by = 0.01),
    confidence.intervals = 0.95,
    study.design = 'case-control',
    population.prevalence = 0.3
  )

dca.list <-
  list(C2.dca, DLGAP1.dca, nomogram)  

pdf("03.DCA.pdf", width = 5, height = 5)
plot_decision_curve(
  dca.list,
  curve.names = c('nomogram', 'C2', 'DLGAP1'),   
  cost.benefit.axis = FALSE,
  col = c(
    "#FF4858",
    "#FF9933",
    "seagreen",
    "blue",
    "purple",
    "#0099DD",
    'pink'
  ),
  confidence.intervals = FALSE,
  standardize = FALSE
)
dev.off()
png(
  "03.DCA.png",
  width = 5,
  height = 5,
  res = 600,
  units = 'in'
)
plot_decision_curve(
  dca.list,
  curve.names = c('nomogram', 'C2', 'DLGAP1'),  
  cost.benefit.axis = FALSE,
  col = c(
    "#FF4858",
    "#FF9933",
    "seagreen",
    "blue",
    "purple",
    "#0099DD",
    'pink'
  ),
  confidence.intervals = FALSE,
  standardize = FALSE
)
dev.off()

library(survivalROC)
library(tidyverse)

res.mul <- coxph(cox_data_prog, data = model_expr)
model_expr$riskScore_prog <-
  predict(res.mul, newdata = model_expr, type = "lp")
train_risk_clinical2 <- model_expr

library(survival)
library(survminer)

riskscore <- function(survival_cancer_df,
                      candidate_genes_for_cox,
                      cox_report) {
  library("dplyr")
  risk_score_table <- survival_cancer_df[, candidate_genes_for_cox]
  for (each_sig_gene in colnames(risk_score_table)) {
    risk_score_table$each_sig_gene <- risk_score_table[, each_sig_gene] *
      (summary(cox_report)$coefficients[each_sig_gene, 1])
  }
  risk_score_table <- cbind(risk_score_table,
                            "total_risk_score" = exp(rowSums(risk_score_table))) %>%
    cbind(survival_cancer_df[, c("OS.time", "OS")])
  risk_score_table <- risk_score_table[, c("OS.time",
                                           "OS",
                                           candidate_genes_for_cox,
                                           "total_risk_score")]
}

multi_ROC_out <- function(time_vector, risk_score_table) {
  library(survivalROC)
  single_ROC <- function(single_time) {
    for_ROC <- survivalROC(
      Stime = risk_score_table$OS.time,
      status = risk_score_table$OS,
      marker = risk_score_table$riskScore_prog,
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

for_multi_ROC <- multi_ROC_out(time_vector = c(365 * seq(1, 3, 1)),
                               risk_score_table = train_risk_clinical2)
table(for_multi_ROC$AUC)
a <-
  unique(data.frame(Time = for_multi_ROC$Time, AUC = for_multi_ROC$AUC))
a
for_multi_ROC$Time <- factor(for_multi_ROC$Time)

library(scales)
library(geomROC)
library(plotROC)

auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time == 365),5][1],2)  
auc_y2 <- round(for_multi_ROC[which(for_multi_ROC$Time == 730),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time == 1095),5][1],2)
auc_y4 <- round(for_multi_ROC[which(for_multi_ROC$Time == 1460),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time == 1825),5][1],2)
auc_y7 <- round(for_multi_ROC[which(for_multi_ROC$Time == 2555),5][1],2)

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
    values = c("#4682B4", "#FF4040", "#20B2AA")
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
  labs(title = "ROC for Prognosis Model") +
  theme(
    panel.background = element_rect(
      colour = 'black',
      size = 1,
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
    labels = c("1 years", "2 years", "3 years"),
    values = c("#4682B4", "#FF4040", "#20B2AA")
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
  labs(title = "ROC for Prognosis Model", family = 'serif') +
  theme(
    panel.background = element_rect(
      colour = 'black',
      size = 1,
      fill = 'white'
    ),
    panel.grid = element_blank(),
    plot.title = element_text(
      size = 15,
      hjust = 0.5,
      face = "bold",
      family = 'serif'
    ),
    axis.title.x = element_text(family = 'serif'),
    axis.title.y = element_text(family = 'serif'),
    legend.text = element_text(family = 'serif'),
    legend.title = element_text(family = 'serif')
  ) +
  annotate(
    "text",
    x = 0.75,
    y = c(0.25, 0.15, 0.05),
    label = c(
      paste('AUC of 1 years =', format(auc_y1, nsmall = 2)),
      paste('AUC of 2 years =', format(auc_y2, nsmall = 2)),
      paste('AUC of 3 years =', format(auc_y3, nsmall = 2))
    ),
    family = 'serif'
  )

ROC
ggsave(
  '06.nomogram_ROC.png',
  ROC,
  width = 5,
  height = 4,
  units = 'in',
  dpi = 600
)
ggsave('06.nomogram_ROC.pdf',
       ROC,
       width = 5,
       height = 4)

save.image("../ws_09_nomogram.RData")
