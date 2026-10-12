rm(list = ls())

setwd("/data/nas1/huangyifeng/project/10_TJZK-10322-3")
if (!dir.exists('./11_CIBERSORT')) {
  dir.create('./11_CIBERSORT')
}
setwd('./11_CIBERSORT')

library(GSVA)
library(magrittr)
library(IOBR)
library(ggplot2)
library(tibble)
library(ggpubr)
library(tidyr)

ssgsea = 1
xcell = 0
ciber = 0
datType = F
bar.cell = 10
box.pic = 1
cor.map = 1
singlecor = 1
cell.cor = 1
group.case <- "high risk"
group.ctrl <- "low risk"

if (T) {
  expr <- read.csv("../00_Rawdata_tcga/Dat_Training_fpkm.csv",
                   check.names = F,
                   row.names = 1)

  group <- read.csv("../00_Rawdata_tcga/Group_Training.csv",
                    check.names = F,
                    row.names = 1)
  risk.score <-
    read.csv("../06_model_validation/risk.csv",
             check.names = F,
             row.names = 1)[, "riskScore", drop = F]


  tryCatch({
    used.sample <- read.csv("../00_Rawdata_tcga/Tumor_Sample.csv",
                            check.names = F,
                            row.names = 1)
    expr <- expr[, used.sample$sample]
    risk <-
      read.csv("../06_model_validation/risk.csv",
               check.names = F,
               row.names = 1)[, "risk", drop = F]
    group <- risk[match(used.sample$sample, rownames(risk)), ,drop = F]
    colnames(group) <- 'group'
  }, error = function(e) {
    next
  })

  group <- rownames_to_column(group, var = "sample")

  Treat_sample <-
    group$sample[which(group$group == paste0(group.case))]
  Control_sample <-
    group$sample[which(group$group == paste0(group.ctrl))]

  if (ciber == 1) {
    tryCatch({
      tiics_result <- readRDS('../00_rawdat/00_temp/ciber_score.rds')
    }, error = function(e) {
      cibersort <-
        deconvo_tme(
          eset = expr,
          method = "cibersort",
          arrays = datType,
          perm = 200
        )
      cibersort <- column_to_rownames(cibersort, var = "ID")
      colnames(cibersort) %<>% gsub("_CIBERSORT", "", .)
      colnames(cibersort) %<>% gsub("_", " ", .)
      cibersort_result <- cibersort[1:22] %>% t
      saveRDS(cibersort_result, "../00_rawdat/00_temp/score.rds")
      tiics_result <- readRDS('../00_rawdat/00_temp/score.rds')


    })
  }
  if (ssgsea == 1) {
    gene_set <-
      read.table(
        "./mmc3.txt",
        header = T,
        sep = "\t"
      )
    gene_list <- split(as.matrix(gene_set)[, 1], gene_set[, 2])
    ssgsea_score <-
      gsva(
        as.matrix(expr),
        gene_list,
        method = "ssgsea",
        kcdf = 'Gaussian',
        abs.ranking = TRUE
      )
    saveRDS(ssgsea_score, 'ssgsea_score.rds')
    tiics_result <- ssgsea_score  %>% as.matrix()
  }
  if (xcell == 1) {
    xcell.dat <-
      deconvo_tme(
        eset = expr,
        method = "xcell",
        arrays = datType,
        perm = 200
      )
    xcell.dat <- column_to_rownames(xcell.dat, var = "ID")
    colnames(xcell.dat) %<>% gsub("_xCell", "", .)
    colnames(xcell.dat) %<>% gsub("_", " ", .)
    xcell.dat <- xcell.dat[1:22] %>% t
    saveRDS(xcell.dat, "ciber_score.rds")
    tiics_result <- readRDS('ciber_score.rds')
  }
  pvalue = padj = log2FoldChange <- matrix(0, nrow(tiics_result), 1)
  for (i in 1:nrow(tiics_result)) {
    pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, Treat_sample], tiics_result[i, Control_sample])$p.value
    log2FoldChange[i, 1] = mean(tiics_result[i, Treat_sample]) - mean(tiics_result[i, Control_sample])
  }
  padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
  rTable <-
    data.frame(log2FoldChange, pvalue, padj, row.names = rownames(tiics_result))
  Treat <-
    signif(apply(tiics_result[rownames(rTable), Treat_sample], 1, mean), 4)
  Control <-
    signif(apply(tiics_result[rownames(rTable), Control_sample], 1, mean), 4)
  rTable <-
    data.frame(Treat, Control, rTable[, c("pvalue", "pvalue", "log2FoldChange")])
  rTable$immune_cell <- rownames(rTable)
  rTable$sig <- ifelse(
    rTable$pvalue < 0.05,
    ifelse(
      rTable$pvalue < 0.01,
      ifelse(
        rTable$pvalue < 0.001,
        ifelse(
          rTable$pvalue < 0.0001,
          paste(rTable$immune_cell, "****",  sep = ""),
          paste(rTable$immune_cell, "***", sep = "")
        ),
        paste(rTable$immune_cell, "**", sep = "")
      ),
      paste(rTable$immune_cell, "*",  sep = "")
    ),
    rTable$immune_cell
  )
  write.table(
    rTable,
    file = "immucell_wilcox_test.xls",
    quote = F,
    row.names = F,
    sep = '\t'
  )
  write.table(
    tiics_result,
    file = 'Score_result.xls',
    sep = '\t',
    row.names = T,
    quote = F
  )
  diff_cibersort_Table <- rTable[which(rTable$pvalue < 0.05), ]
  sigcell <- rownames(diff_cibersort_Table)
  saveRDS(sigcell, 'sigCell.rds')
  write.table(
    diff_cibersort_Table,
    file = 'DiffTable.xls',
    quote = F,
    sep = '\t'
  )
}

bar.dat <- as.data.frame(tiics_result)
bar.dat$cell <- rownames(bar.dat)
box.dat <- gather(bar.dat, key = sample, value = 'score', -c("cell"))
box.dat$group <-
  ifelse(box.dat$sample %in% Control_sample,
         paste0(group.ctrl),
         paste0(group.case))

fixed_colors <- c(
  "#8DD3C7",
  "#FFFFB3",
  "#BEBADA",
  "#FB8072",
  "#80B1D3",
  "#FDB462",
  "#B3DE69",
  "#FCCDE5",
  "#D9D9D9",
  "#BC80BD",
  "#CCEBC5",
  "#FFED6F",
  "#A6CEE3",
  "#1F78B4",
  "#B2DF8A",
  "#33A02C",
  "#FB9A99",
  "#E31A1C",
  "#FDBF6F",
  "#FF7F00",
  "#CAB2D6",
  "#6A3D9A",
  "#FF7F7F",
  "#999999",
  "#B15928",
  "#FFCC00",
  "#66C2A5",
  "#FC8D62"
)

library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(7, "Paired"))

pdf('cibersort.box.pdf', width = 10, height = 8)

box.dat$score <- abs(box.dat$score)
ggplot(box.dat, aes(x = sample, y = abs(score), fill = cell)) +
  geom_bar(position = 'stack', stat = 'identity') +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  labs(x = '',
       y = 'scores',
       fill = '') +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'top'
  ) +

  scale_fill_manual(values = fixed_colors) +
  facet_grid( ~ box.dat$group, scales = "free", space = "free")
dev.off()

png(
  'cibersort.box.png',
  w = 10,
  h = 8,
  units = 'in',
  res = 600
)

box.dat$score <- abs(box.dat$score)
ggplot(box.dat, aes(x = sample, y = abs(score), fill = cell)) +
  geom_bar(position = 'stack', stat = 'identity') +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  labs(x = '',
       y = 'scores',
       fill = '') +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'top'
  ) +

  scale_fill_manual(values = fixed_colors) +
  facet_grid( ~ box.dat$group, scales = "free", space = "free")
dev.off()

if (cell.cor == 1) {
  cellcor.dat <- tiics_result %>% t
  sig.cell <- readRDS('sigCell.rds')
  cellcor.dat <- cellcor.dat[, sig.cell]

  library(ggcorrplot)
  library(ggthemes)
  library(psych)

  cor_result <- corr.test(cellcor.dat, method = 'spearman')
  cor_matrix <- cor_result$r
  p_matrix <- cor_result$p

  write.csv(cor_matrix, file = 'cor.index.csv')
  write.csv(p_matrix, file = 'cor.pvalue.csv')

  corr <- as.data.frame(round(cor_matrix, 3))


  p <-
    ggcorrplot(
      corr,
      method = "circle",
      hc.order = TRUE,
      hc.method = "complete",
      outline.col = 'white',
      ggtheme = theme_bw(),
      type = "full",
      lab = TRUE,
      lab_size = 2,
      insig = "blank"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        size = 14,
        hjust = 1,
        family = "Times",
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 14,
        family = "Times",
        face = "bold"
      )
    )
  print(p)
  nump <- 10
  ggsave(filename = "cell.corr.pdf",
         width = nump,
         height = nump,
         p)
  ggsave(filename = "cell.corr.png",
         width = nump,
         height = nump,
         p)
}

if (box.pic == 1) {
  library(ggpubr)
  dat <-
    tiics_result %>% t %>% as.data.frame() %>% tibble::rownames_to_column(var = "sample")
  dat <- merge(group, dat, by = "sample")
  dat2 <-
    tidyr::gather(dat, ImmuneCell, Score, -c("sample", "group"))
  dat2$Score <- as.numeric(dat2$Score)
  dat2$group <-
    factor(dat2$group, levels = c(paste0(group.ctrl), paste0(group.case)))
  x1 = 1:2
  if (sample(x1, 1) == 1) {
    p <- ggboxplot(
      dat2,
      x = "ImmuneCell",
      y = "Score",
      color = "group",
      palette = c("seagreen4", "darkorange3"),
      outlier.shape = NA,
      bxp.errorbar = T
    ) +
      stat_compare_means(aes(group = group), label = "p.signif") +
      labs(x = "", y = "Score", color = "") +
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
    p <-  ggplot(dat2,
                 aes(
                   x = ImmuneCell,
                   y = Score,
                   fill = group,
                   outlier.shape = NA,
                   bxp.errorbar = T
                 )) +
      geom_boxplot(width = 0.8,
                   alpha = 1,
                   position = position_dodge(0.9)) +


      stat_compare_means(aes(group = group), label = "p.signif") +
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
      labs(x = "", y = "Score") + ggsci::scale_fill_d3()
  }
  ggsave(filename = "ImmuneCell_boxplot.pdf",
         width = 12,
         height = 8,
         p)
  ggsave(filename = "ImmuneCell_boxplot.png",
         width = 12,
         height = 8,
         p)
}

if (cor.map == 1) {
  hubgene <-
    read.csv("../06_model_validation/lasso_genes.csv",
             check.names = F,
             row.names = 1)
  hubgene <- hubgene$genes
  library(tibble)
  library(magrittr)
  library(corrplot)
  library(reshape2)

  group <- na.omit(group)
  filte_dat <- expr[, match(group$sample, colnames(expr))]
  filte_dat <-
    filte_dat %>% as.data.frame() %>% rownames_to_column(., var = "Gene.name")
  tem_12 <- filte_dat[match(hubgene, filte_dat$Gene.name), ]
  keyRNA_exp <- tem_12[, -1]
  rownames(keyRNA_exp) <- tem_12$Gene.name
  keyRNA_exp <- keyRNA_exp %>% t %>% as.data.frame()
  keyRNA_exp <- cbind(risk.score, keyRNA_exp)



  sigcell_exp <- tiics_result  %>% t %>% as.data.frame()
  sig.cell <- readRDS('sigCell.rds')
  sigcell_exp <- sigcell_exp[group$sample, sig.cell]
  cor <- stats::cor

  library(psych)

  cor_result2 <- corr.test(sigcell_exp, keyRNA_exp, method = 'spearman')
  cor_matrix2 <- cor_result2$r
  p_matrix2 <- cor_result2$p


  # corr <-
  #   round(cor_matrix2, 3)
  # cor_p <-
  #   round(p_matrix2, 4)

  corr <- cor_matrix2
  cor_p <- p_matrix2

  write.csv(corr,file = "cor_with_cell.csv")
  write.csv(cor_p, file = "corp_with_cell.csv")




  corr <- corr[, colnames(corr) != "riskScore"]
  cor_p <- cor_p[, colnames(cor_p) != "riskScore"]

  library(reshape2)
  temcor <- melt(corr, value.name = "cor")
  temcor <-  temcor[abs(temcor$cor) > 0.3, ]
  temcorp <- melt(cor_p, value.name = "pvalue")
  temcorp <- temcorp[temcorp$pvalue < 0.05, ]
  colnames(temcorp)[1:2] <- c("Var1", "Var2")
  tem <- merge(temcor, temcorp, all = T) %>% na.omit()
  tem <- tem[order(tem$Var2, decreasing = F), ]
  tem$pvalue[tem$pvalue == 0] <- "<0.001"

  library(RColorBrewer)
  library(pheatmap)
  col1 <-
    colorRampPalette(c("#008B45B2" , "#D1E5F0", "#EE0000B2"), alpha = TRUE)(100)
  cor <- corr
  pvalue <- cor_p
  myfun <- function(pval) {
    stars = ""
    if (pval <= 0.001)
      stars = "***"
    if (pval > 0.001 & pval <= 0.01)
      stars = "**"
    if (pval > 0.01 & pval <= 0.05)
      stars = "*"
    if (pval > 0.05 & pval <= 0.1)
      stars = ""
    stars
  }
  heatmap <- melt(cor) %>%
    mutate(pvalue = melt(pvalue)[, 3]) %>%
    mutate(signif = sapply(pvalue, function(x)
      myfun(x)))
  cor_char <- as.character(round(heatmap$value, 2))

  colnames(heatmap) <- c("Cell", "Gene", "cor", "pvalue", "signif")

  p1 <- ggplot(heatmap, aes(Gene, Cell, fill = cor)) +
    geom_tile() + theme_minimal() +
    geom_text(
      aes(label = signif),
      size = 5,
      color = "black",
      hjust = 0.5,
      vjust = 0.5
    ) +
    labs(x = NULL, y = NULL, color = NULL) +
    scale_fill_viridis_c() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(
      text = element_text(family = "Times"),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_rect(
        fill = NA,
        color = "grey70",
        size = 1,
        linetype = "solid"
      )
    ) +
    guides(color = guide_colorbar(
      direction = "vertical",
      reverse = F,
      barwidth = unit(.5, "cm"),
      barheight = unit(9, "cm")
    )) +
    scale_fill_gradientn(
      colours = col1,
      limits = c(-1, 1),
      name = paste0(
        "*    p < 0.05",
        "\n\n",
        "**  p < 0.01",
        "\n\n",
        "*** p < 0.001",
        "\n\n",
        "Correlation"
      ),
      guide = guide_colourbar(
        ticks = T,
        nbin = 50,
        barheight = 10,
        label = T,
        barwidth = 1
      )
    ) +
    theme(legend.position = "right") +
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
        size = 16,
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
    labs(x = "", y = "") + coord_flip()
  pdf(
    file = paste0("correlation.pdf"),
    width = 12,
    height = 8
  )
  a <- dev.cur()
  png(
    file = paste0("correlation.png"),
    width = 12,
    height = 8,
    units = "in",
    res = 300
  )
  dev.control("enable")
  print(p1)
  dev.copy(which = a)
  dev.off()
  dev.off()
}

library(tibble)
library(magrittr)
library(corrplot)
library(reshape2)
cor_dat <- cbind(risk.score, sigcell_exp)
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
  library("viridis")
  library(ggpubr)

  pdf(
    file = paste0("01.riskScore_Corr_Cells.pdf"),
    width = 8,
    height = 7
  )
  a <- dev.cur()
  png(
    file = paste0("01.riskScore_Corr_Cells.png"),
    width = 8,
    height = 7,
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

if (singlecor == 1) {
  risk.group <- read.csv("../06_model_validation/risk.csv", row.names = 1,check.names = F)
  sat.act <- cbind(sigcell_exp, risk.group)
  cellname <- colnames(sigcell_exp)
  colnames(sat.act)[c(1:ncol(sigcell_exp))] <-
    paste0('cell', 1:ncol(sigcell_exp))

  for (i in 1:ncol(sigcell_exp)) {
    p <- ggplot(sat.act, aes_string(paste0('cell', i), 'riskScore')) +
      geom_point(aes(colour = risk)) +
      stat_smooth(method = lm) +
      labs(x = cellname[i], y = "riskScore") +
      theme_bw() +
      stat_cor(data = sat.act, method = "pearson") +
      theme(legend.position = "top") +
      theme(
        axis.title.x = element_text(
          size = 24,
          family = "Times",
          face = "bold"
        ),
        axis.text.x = element_text(
          size = 17,
          family = "Times",
          face = "bold"
        ),
        axis.title.y = element_text(
          size = 24,
          family = "Times",
          face = "bold"
        ),
        axis.text.y = element_text(
          size = 17,
          family = "Times",
          face = "bold"
        ),
        legend.text = element_text(
          size = 20,
          family = "Times",
          face = "bold"
        ),
        legend.title = element_text(
          size = 20,
          family = "Times",
          face = "bold"
        )
      ) +
      ggsci::scale_fill_d3()
    assign(paste0('p', i), p)
  }
  library(cowplot)
  pall <- plot_grid(p1, p2, p3, p4, nrow = 2)
  ggsave(width = 5, height = 5, 'cor_ALL.pdf', p4)
  ggsave(width = 5, height = 5, 'cor_ALL.png', p4)
}

save.image("../ws_11_CIBERSORT.RData")
