.libPaths("/data/nas2/Software/miniconda3/envs/public_R/lib/R/library")

rm(list = ls())
setwd('/data/nas1/huangyifeng/project/10_TJZK-10322-3')
if (!dir.exists("./00_singlecell")) {
  dir.create("./00_singlecell")
}
setwd("./00_singlecell")

library(Seurat)
library(ggplot2)
library(GEOquery)
library(data.table)
library(tidyverse)

library(stringr)

folders = list.files('./GSE162454_RAW//')
folders
dir = paste0('GSE162454_RAW///', folders, '/')
dir

scList <- list()
for (i in c(1:6)) {
  count =  CreateSeuratObject(
    counts = Read10X(dir[i]),
    min.cells = 3,

    min.features = 200
  )
  scList[[i]] <- count
}

scRNA <- merge(
  scList[[1]],
  y = c(scList[[2]], scList[[3]], scList[[4]], scList[[5]], scList[[6]]),
  add.cell.ids = folders,
  project = "GSE162454"
)

table(scRNA$orig.ident)
head(scRNA@meta.data)
length(colnames(scRNA))
length(rownames(scRNA))
scRNA <- JoinLayers(scRNA)
head(scRNA@meta.data)

gset <- getGEO(
  "GSE162454",
  destdir = '.',
  GSEMatrix = T,
  getGPL = F
)

a = gset[[1]]
pd <- pData(a)
group <- data.frame(sample = pd$geo_accession,
                    group = pd$source_name_ch1)
library(stringr)
table(group$group)
group$group <- ifelse(grepl("sham", group$group),
                      'Control',
                      ifelse(grepl("Tumor tissue", group$group), 'Tumor', 'Other_Disease'))

table(group$group)

group$group
write.csv(group, file = "group_GSE162454.csv", row.names = F)

library(tidyverse)
smps = Cells(scRNA) %>% str_sub(1, 10)
table(smps)
grps = group
group.label = grps$group[match(smps, grps$sample)]
names(group.label) = Cells(scRNA)
group.label
scRNA = AddMetaData(scRNA, group.label, col.name = 'group')

Patients_ID <-
  tstrsplit(rownames(scRNA@meta.data), "_", fixed = TRUE)[[1]]
scRNA$Patients_ID <- Patients_ID
table(scRNA$group)

scRNA[["percent.mt"]] <-
  PercentageFeatureSet(scRNA, pattern = "^MT-")

scRNA[["percent.rb"]] <-
  PercentageFeatureSet(scRNA, pattern = "^RP[SL]")

scRNA[["percent.hb"]] <-
  PercentageFeatureSet(scRNA, pattern = "^HB[SL]")

head(scRNA@meta.data)

system.time(save(scRNA, file = "01.scRNA_orig.Rdata"))

rm(list = ls())
library(Seurat)
library(tidyverse)
setwd('/data/nas1/huangyifeng/project/10_TJZK-10322-3/00_singlecell')

if (!dir.exists("./01_QC")) {
  dir.create("./01_QC")
}
setwd("./01_QC")
load('/data/nas1/huangyifeng/project/10_TJZK-10322-3/00_singlecell/01.scRNA_orig.Rdata')

length(colnames(scRNA))
length(rownames(scRNA))

plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")

library(patchwork)
plots = list()
for (i in seq_along(plot.featrures)) {
  plots[[i]] = VlnPlot(
    scRNA,
    group.by = "Patients_ID",
    cols = c(
      '#74AED4',
      '#D3E2B7',
      '#CFAFD4',
      '#F7C97E',
      "#A6D854",
      "#FFD92F"
    ),

    pt.size = 0,
    features = plot.featrures[i]
  ) + NoLegend()
}
violin <-
  wrap_plots(
    plots = plots,
    nrow = 1,
    ncol = 3,
    cols = palettes(category = "random", 22, show_col = F)
  )
violin

ggsave(
  "01.vlnplot_before_qc.pdf",
  plot = violin,
  width = 15,
  height = 6
)
ggsave(
  "01.vlnplot_before_qc.png",
  plot = violin,
  width = 15,
  height = 6,
  dpi = 600
)

pctMT = 10

minGene = 200
maxGene = 6000

scRNA <- subset(scRNA,
                subset =

                  nFeature_RNA > minGene &
                  nFeature_RNA < maxGene &
                  percent.mt < pctMT)
length(colnames(scRNA))
length(rownames(scRNA))

plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")

plots = list()
for (i in seq_along(plot.featrures)) {
  plots[[i]] = VlnPlot(
    scRNA,
    group.by = "Patients_ID",
    cols = c(
      '#74AED4',
      '#D3E2B7',
      '#CFAFD4',
      '#F7C97E',
      "#A6D854",
      "#FFD92F"
    ),
    pt.size = 0,
    features = plot.featrures[i]
  ) + NoLegend()
}
violin <-
  wrap_plots(plots = plots, nrow = 1, ncol = 3)
violin
ggsave(
  "02.vlnplot_after_qc.pdf",
  plot = violin,
  width = 15,
  height = 6
)
ggsave(
  "02.vlnplot_after_qc.png",
  plot = violin,
  width = 15,
  height = 6,
  dpi = 600
)

system.time(save(scRNA, file = "02.scRNA_qc.Rdata"))

rm(list = ls())
setwd('/data/nas1/huangyifeng/project/10_TJZK-10322-3/00_singlecell')
if (!dir.exists("./02_PCA")) {
  dir.create("./02_PCA")
}
setwd("./02_PCA")
load('../01_QC/02.scRNA_qc.Rdata')
library(Seurat)

for (i in 1:6) {
  scRNA1 <- subset(scRNA, Patients_ID == unique(scRNA$Patients_ID)[i])
  scRNA2 <-
    NormalizeData(scRNA1,
                  normalization.method = "LogNormalize",
                  scale.factor = 10000)
  scRNA3 <-
    FindVariableFeatures(scRNA2, selection.method = "vst", nfeatures = 2000)
  assign(unique(scRNA$Patients_ID)[i],  scRNA3)
}

scRNA.anchors <-
  FindIntegrationAnchors(
    object.list = list(
      GSM5155200,
      GSM5155199,
      GSM5155198,
      GSM4952365,
      GSM4952364,
      GSM4952363
    ),
    dims = 1:20
  )
scRNA.integrated <-
  IntegrateData(anchorset = scRNA.anchors, dims = 1:20)
system.time(save(scRNA.integrated, file = "scRNA_integrated.Rdata"))

rm(list = ls())
setwd('/data/nas1/huangyifeng/project/10_TJZK-10322-3/00_singlecell')
if (!dir.exists("./02_PCA")) {
  dir.create("./02_PCA")
}
setwd("./02_PCA")
library(Seurat)
load('../01_QC/02.scRNA_qc.Rdata')
load('../02_PCA/scRNA_integrated.Rdata')

scRNA.norm <-
  NormalizeData(scRNA,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

all.genes <- rownames(scRNA.norm)
scRNA.norm <-
  FindVariableFeatures(scRNA.norm,
                       selection.method = "vst",
                       nfeatures = 2000)

top10 <- head(VariableFeatures(scRNA.norm), 10)

p <- VariableFeaturePlot(scRNA.norm)
p <- LabelPoints(plot = p,
                 points = top10,
                 repel = T)
p

ggsave(filename = "01.feature_selection.pdf", p, w = 6, h = 6)
ggsave(filename = "01.feature_selection.png", p, w = 6, h = 6)
dev.off()

scRNA.nor.sca <- ScaleData(scRNA.norm)

scRNA.norm.pca <-
  RunPCA(scRNA.nor.sca,
         features = VariableFeatures(object = scRNA.nor.sca),
         npcs = 50)
print(scRNA.norm.pca[["pca"]], dims = 1:50, nfeatures = 5)
VizDimLoadings(scRNA.norm.pca, dims = 1:2, reduction = "pca")

table(scRNA$Patients_ID)
table(scRNA$orig.ident)

pdf("02.PCA.pdf", w = 7, h = 5)
DimPlot(
  scRNA.norm.pca,
  reduction = "pca",
  group.by = "orig.ident",
  cols = c(
    '#D3E2B7',
    '#CFAFD4',
    '#F7C97E',
    "#A6D854",
    "#FFD92F",
    "#0070C0"
  )
)
dev.off()
png(
  "02.PCA.png",
  w = 7,
  h = 5,
  units = 'in',
  res = 600
)
DimPlot(
  scRNA.norm.pca,
  reduction = "pca",
  group.by = "orig.ident",
  cols = c(
    '#D3E2B7',
    '#CFAFD4',
    '#F7C97E',
    "#A6D854",
    "#FFD92F",
    "#0070C0"
  )
)
dev.off()
DimHeatmap(
  scRNA.norm.pca,
  dims = 1:5,
  cells = 500,
  balanced = TRUE
)

scRNA.norm.pca <-
  JackStraw(scRNA.norm.pca, num.replicate = 100, dims = 50)
scRNA.norm.pca2 <- ScoreJackStraw(scRNA.norm.pca, dims = 1:50)
plot_pca <- JackStrawPlot(scRNA.norm.pca2, dims = 1:50)
plot_pca
ggsave(filename = '03.pca_cluster.pdf',
       plot_pca,
       w = 13,
       h = 8)
ggsave(
  filename = '03.pca_cluster.png',
  plot_pca,
  w = 11,
  h = 7,
  dpi = 600
)
plot_elbow <- ElbowPlot(scRNA.norm.pca, ndims = 50)
plot_elbow
ggsave(filename = '04.pca_sd.pdf',
       plot_elbow,
       w = 6,
       h = 5)
ggsave(
  filename = '04.pca_sd.png',
  plot_elbow,
  w = 5,
  h = 4,
  dpi = 600
)

scRNA.norm.pca.c <-
  FindNeighbors(scRNA.norm.pca, dims = 1:50)
scRNA.norm.pca.c <-
  FindClusters(scRNA.norm.pca.c, resolution = 0.4)
UMAP <- RunUMAP(scRNA.norm.pca.c, dims = 1:50)
system.time(save(UMAP, file = "UMAP.Rdata"))

theme.set = theme(
  axis.title = element_text(
    size = 20,
    face = "bold",
    family = "Times"
  ),
  axis.text.x = element_text(
    size = 14,
    face = "bold",
    family = "Times"
  ),
  axis.text.y = element_text(
    size = 14,
    face = "bold",
    family = "Times"
  ),
  legend.text = element_text(
    size = 16,
    face = "bold",
    family = "Times"
  ),
  legend.title = element_blank(),
  text = element_text(family = "Times")
)

library(RColorBrewer)
p11 <- DimPlot(
  UMAP,
  reduction = "umap",
  group.by = "Patients_ID",
  cols = brewer.pal(n = 8, name = "Set2"),
  pt.size = 0.5,
  label = FALSE,
  repel = TRUE
) + theme.set
p22 <- DimPlot(
  UMAP,
  reduction = "umap",
  group.by = "group",
  cols = brewer.pal(n = 8, name = "Set2"),
  pt.size = 0.5,
  label = FALSE,
  repel = TRUE
) + theme.set

p11 + p22

pdf(file = "04.UMAP_after.pdf",
    width = 17,
    height = 5)
p11 + p22
dev.off()

png(
  file = "04.UMAP_after.png",
  width = 17,
  height = 5,
  units = 'in',
  res = 600
)
p11 + p22
dev.off()

pdf(file = "02.umap_final.pdf",
    width = 8,
    height = 7)
DimPlot(UMAP, reduction = "umap", label = T)
dev.off()
png(
  file = "02.umap_final.png",
  width = 8,
  height = 7,
  units = 'in',
  res = 600
)
DimPlot(UMAP, reduction = "umap", label = T)
dev.off()

rm(list = ls())
setwd('/data/nas1/huangyifeng/project/10_TJZK-10322-3/00_singlecell')
if (!dir.exists("./02_PCA")) {
  dir.create("./02_PCA")
}
setwd("./02_PCA")

load('UMAP.Rdata')

library(Seurat)

all.markers <-
  FindAllMarkers(
    UMAP,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.5,
    test.use = 'wilcox',
    return.thresh = 0.01
  )

write.csv(all.markers,
          file = '01.AllMarkers.csv',
          row.names = F,
          quote = F)
library(SingleR)
library(Seurat)

scRNA_nor_singleR <- GetAssayData(UMAP, layer = "data")

library(celldex)

hpca.se <- celldex::BlueprintEncodeData()
# yes

hpca.se
table(hpca.se@colData@listData$label.main)

table(UMAP@meta.data$seurat_clusters)

clusters = UMAP@meta.data$seurat_clusters
pred.hesc <- SingleR(
  test = scRNA_nor_singleR,
  ref = hpca.se,
  labels = hpca.se$label.main,
  method = "cluster",
  clusters = clusters,
  quantile = 0.2,

  assay.type.test = "logcounts",
  assay.type.ref = "logcounts"
)

table(pred.hesc$labels)

celltype = data.frame(
  ClusterID = rownames(pred.hesc),
  celltype = pred.hesc$labels,
  stringsAsFactors = F
)

write.csv(celltype, "celltype.csv")
library(tidyverse)

UMAP@meta.data$celltype = celltype[match(clusters, celltype$ClusterID), 'celltype']
anno_plot <-
  DimPlot(UMAP,
          reduction = "umap",
          group.by = "celltype",
          label = T)
anno_plot
library(ggplot2)
ggsave(filename = '01.celltype.pdf',
       anno_plot,
       w = 8,
       h = 6)
ggsave(
  filename = '01.celltype.png',
  anno_plot,
  w = 8,
  h = 6,
  dpi = 600
)
phe = UMAP@meta.data
table(phe$celltype)

cell_type_stat <- as.data.frame(sort(table(phe$celltype)))
write.csv(cell_type_stat,
          file = "02.cell_type_stat.csv",
          quote = F)

system.time(save(UMAP, file = "UMAP.celltype.Rdata"))

anno_plot2 <-
  DimPlot(
    UMAP,
    reduction = "umap",
    group.by = "celltype",
    split.by = 'group',
    label = T
  ) +
  labs(title = "")
anno_plot2
ggsave(filename = '02.celltype(split).pdf',
       anno_plot2,
       w = 10,
       h = 6)
ggsave(
  filename = '02.celltype(split).png',
  anno_plot2,
  w = 10,
  h = 6,
  dpi = 600
)

new.cluster.ids = as.character(celltype$celltype)
new.cluster.ids
names(new.cluster.ids) <- celltype$ClusterID
UMAP <- RenameIdents(UMAP, new.cluster.ids)

all.markers <-
  FindAllMarkers(
    UMAP,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.5,
    test.use = 'wilcox',
    return.thresh = 0.01
  )
write.csv(all.markers,
          file = '03.MarkersBetweenCells.csv',
          row.names = F,
          quote = F)

table(UMAP$celltype)
gene1 <-
  subset(all.markers,
         all.markers$avg_log2FC > 1 & all.markers$p_val_adj < 0.05)

EN1 <- gene1[gene1$cluster == "Astrocytes", ]
EN <-
  EN1[rownames(EN1) %in% rownames(head(EN1[order(EN1$avg_log2FC, decreasing = T), ], 5)), ]

B1 <- gene1[gene1$cluster == "B-cells", ]
B <-
  B1[rownames(B1) %in% rownames(head(B1[order(B1$avg_log2FC, decreasing = T), ], 5)), ]

CH1 <- gene1[gene1$cluster == "CD8+ T-cells", ]
CH <-
  CH1[rownames(CH1) %in% rownames(head(CH1[order(CH1$avg_log2FC, decreasing = T), ], 5)), ]

mo1 <- gene1[gene1$cluster == "Chondrocytes", ]
mo <-
  mo1[rownames(mo1) %in% rownames(head(mo1[order(mo1$avg_log2FC, decreasing = T), ], 5)), ]

MA1 <- gene1[gene1$cluster == "Endothelial cells", ]
MA <-
  MA1[rownames(MA1) %in% rownames(head(MA1[order(MA1$avg_log2FC, decreasing = T), ], 5)), ]

NK1 <- gene1[gene1$cluster == "Fibroblasts", ]
NK <-
  NK1[rownames(NK1) %in% rownames(head(NK1[order(NK1$avg_log2FC, decreasing = T), ], 5)), ]

SM1 <- gene1[gene1$cluster == "HSC", ]
SM <-
  SM1[rownames(SM1) %in% rownames(head(SM1[order(SM1$avg_log2FC, decreasing = T), ], 5)), ]

AD1 <- gene1[gene1$cluster == "Macrophages", ]
AD <-
  AD1[rownames(AD1) %in% rownames(head(AD1[order(AD1$avg_log2FC, decreasing = T), ], 5)), ]

C81 <- gene1[gene1$cluster == "Monocytes", ]
C8 <-
  C81[rownames(C81) %in% rownames(head(C81[order(C81$avg_log2FC, decreasing = T), ], 5)), ]

C41 <- gene1[gene1$cluster == "NK cells", ]
C4 <-
  C41[rownames(C41) %in% rownames(head(C41[order(C41$avg_log2FC, decreasing = T), ], 5)), ]

gene <- rbind(EN, B, CH, mo, MA, NK, SM, C8, C4, AD)

write.csv(gene, "03.Markergenes_top5.csv")

gene <- read.csv("03.Markergenes_top5.csv", row.names = 1)
library(ggplot2)
theme.set = theme(
  axis.title = element_text(
    size = 20,
    face = "bold",
    family = "Times"
  ),
  axis.text.x = element_text(
    size = 14,
    face = "bold",
    angle = 45,
    hjust = 1,
    family = "Times"
  ),
  axis.text.y = element_text(
    size = 14,
    face = "bold",
    family = "Times"
  ),
  legend.text = element_text(
    size = 16,
    face = "bold",
    family = "Times"
  ),
  legend.title = element_blank(),
  text = element_text(family = "Times")
)
library(Seurat)
DotPlot(UMAP, features = gene$gene) + theme.set + coord_flip() +
  scale_color_gradientn(colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))

pdf(
  file = "04.DotPlot.pdf",
  width = 9,
  height = 15,
  family = 'Times'
)
DotPlot(UMAP, features = gene$gene) + theme.set + coord_flip() +
  scale_color_gradientn(colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))
dev.off()

png(
  file = "04.DotPlot.png",
  width = 9,
  height = 15,
  family = 'Times',
  units = 'in',
  res = 600
)
DotPlot(UMAP, features = gene$gene) + theme.set + coord_flip() +
  scale_color_gradientn(colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))
dev.off()

system.time(save(UMAP, file = "UMAP.Rdata"))
load('UMAP.Rdata')
library(tidyverse)
dat_cell <- as.data.frame(table(UMAP$group, UMAP$celltype))

colnames(dat_cell) <- c('id', 'cell', 'number')
dat_cell2 <- spread(dat_cell, cell, number)

dat_cell2 <- dat_cell2 %>% column_to_rownames(., var = "id")
dat_cell2 <-
  dat_cell2 %>% mutate(rowsum = base::rowSums(.[1:dim(dat_cell2)[2]]))

dat_cell3 <- dat_cell2

for (i in 1:dim(dat_cell3)[1]) {
  for (j in 1:(dim(dat_cell3)[2] - 1)) {
    dat_cell3[i, j] <- dat_cell3[i, j] / dat_cell3[i, dim(dat_cell3)[2]]
  }
}
dat_cell3 <- dat_cell3[, -dim(dat_cell3)[2]]
dat_cell3 <- dat_cell3 %>% rownames_to_column('sample')
write.csv(dat_cell3, '05.cell_type_group.csv')

library(reshape2)
library(RColorBrewer)
mydata <-
  melt(
    dat_cell3,
    id.vars = "sample",
    variable.name = "Celltype",
    value.name = "Quantity"
  )
ggplot(mydata, aes(x = Quantity, y = sample,
                   fill = Celltype)) +
  geom_bar(position = "fill",
           stat = "identity",
           alpha = 0.7) +
  theme_bw() +
  scale_fill_manual(values = brewer.pal(11, "Set3")) +
  theme(
    axis.title.x = element_text(
      size = 22,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    axis.text.x = element_text(
      size = 18,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    axis.title.y = element_text(
      size = 22,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    axis.text.y = element_text(
      size = 18,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    legend.title = element_text(
      size = 20,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    legend.text = element_text(
      size = 14,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    title = element_text(
      size = 20,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    strip.text = element_text(
      size = 14,
      color = 'black',
      family = "Times",
      face = "bold"
    )
  ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Proportion", y = "Celltype", fill = "")

ggsave('05.Proportion.pdf', w = 14, h = 20)
ggsave('05.Proportion.png', w = 14, h = 20)


table(UMAP$group)
# Tumor
# 32737
scRNA_case <- subset(UMAP, group %like% 'Tumor')
phe = scRNA_case@meta.data

table(phe$celltype)
cell_type_case <- as.data.frame(sort(table(phe$celltype)))
colnames(cell_type_case) <- c('Celltype', 'Tumor')

cell_type <- cell_type_case

write.csv(cell_type, '06.cell_type_group.csv')
library(reshape2)
library(plyr)
library(RColorBrewer)
library(ggplot2)
mydata <-
  melt(
    cell_type,
    id.vars = "Celltype",
    variable.name = "Group",
    value.name = "Quantity"
  )
ggplot(mydata, aes(x = Group, y = Quantity, fill = Celltype)) +
  geom_bar(position = "fill",
           stat = "identity",
           alpha = 0.7) +
  theme_bw() +
  scale_fill_manual(values = brewer.pal(11, "Set3")) +
  theme(
    axis.title.x = element_text(
      size = 22,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    axis.text.x = element_text(
      size = 18,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    axis.title.y = element_text(
      size = 22,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    axis.text.y = element_text(
      size = 18,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    legend.title = element_text(
      size = 20,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    legend.text = element_text(
      size = 14,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    title = element_text(
      size = 20,
      color = 'black',
      face = "bold",
      family = 'Times'
    ),
    strip.text = element_text(
      size = 14,
      color = 'black',
      family = "Times",
      face = "bold"
    )
  ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Celltype", y = "Proportion", fill = "")

ggsave('06.Proportion_group.pdf', w = 10, h = 8)
ggsave('06.Proportion_group.png', w = 10, h = 8)

