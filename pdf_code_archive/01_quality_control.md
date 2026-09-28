# Ambient RNA, QC and doublet filtering

PDF transcription only. The PDF may wrap code or omit glyphs;
use the curated scripts for executable analysis.

## Source PDF page 5

```text
#ck_fcr SoupX
sc = load10X("/home/xiehangfei/scRNAseq/ck_fcr_10000cells/outs")
sc = autoEstCont(sc)
out = adjustCounts(sc)
dir.create("ck_fcr_matrix_soupx")
writeMM(out, file="ck_fcr_matrix_soupx/matrix.mtx")
#fcr SoupX
sc = load10X("/home/xiehangfei/scRNAseq/fcr_10000cells/outs")
sc = autoEstCont(sc)
out = adjustCounts(sc)
dir.create("fcr_matrix_soupx")
writeMM(out, file="fcr_matrix_soupx/matrix.mtx")
#ck_fcr_matrix_soupx
gzip matrix.mtx
cp ../ck_fcr_10000cells/outs/filtered_feature_bc_matrix/features.tsv.gz ./
cp ../ck_fcr_10000cells/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ./

#fcr_matrix_soupx
gzip matrix.mtx
cp ../fcr_10000cells/outs/filtered_feature_bc_matrix/features.tsv.gz ./
cp ../fcr_10000cells/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ./
```

## Source PDF page 6

```text
#ck_fcr filter
ck_fcr.data = Read10X("ck_fcr_matrix_soupx")
ck_fcr <- CreateSeuratObject(counts = ck_fcr.data, project = "ck_fcr",min.cells = 3, min.features = 200)
rownames(ck_fcr[["RNA"]]) <- sub("^gene:", "", rownames(ck_fcr[["RNA"]]))
keep_genes <- grepl("^TraesCS[1-8][ABD]", rownames(ck_fcr[["RNA"]]))
ck_fcr <- subset(ck_fcr, features = rownames(ck_fcr[["RNA"]])[keep_genes])

gene_data <- FetchData(ck_fcr, vars = "TraesCS8B02G990504")
expressed_cells <- sum(gene_data > 0)
total_cells <- nrow(gene_data)
print(paste("           :", expressed_cells))
print(paste("   :", total_cells))
print(paste("   :", round(expressed_cells / total_cells * 100, 2),"%"))
summary(gene_data)

gene_data <- FetchData(ck_fcr, vars = "TraesCS8D02G020913")
expressed_cells <- sum(gene_data > 0)
total_cells <- nrow(gene_data)
print(paste("           :", expressed_cells))
print(paste("   :", total_cells))
print(paste("   :", round(expressed_cells / total_cells * 100, 2),"%"))
summary(gene_data)

ck_fcr[["percent.mito"]] <- PercentageFeatureSet(ck_fcr, features = "TraesCS8B02G990504")
ck_fcr[["percent.chloro"]] <- PercentageFeatureSet(ck_fcr, features ="TraesCS8D02G020913")
head(ck_fcr@meta.data)
p <- VlnPlot(ck_fcr,
features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.chloro"),
ncol = 4,
pt.size = 0.1)

ggsave("ck_fcr_raw_QC_VlnPlot.pdf", p, width = 12, height = 8)

p1 <- FeatureScatter(ck_fcr, feature1 = "nCount_RNA", feature2 = "percent.mito") +
```

## Source PDF page 7

```text
geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
geom_hline(yintercept = 10, linetype = "dashed", color = "blue") +
geom_hline(yintercept = 15, linetype = "dashed", color = "green") +
annotate("text", x = Inf, y = 5, label = "5%", color = "red", hjust =1.1, vjust = -0.5, fontface = "bold") +
annotate("text", x = Inf, y = 10, label = "10%", color = "blue", hjust= 1.1, vjust = -0.5, fontface = "bold") +
annotate("text", x = Inf, y = 15, label = "15%", color = "green", hjust = 1.1, vjust = -0.5, fontface = "bold") +
ggtitle(" ")

p2 <- FeatureScatter(ck_fcr, feature1 = "nCount_RNA", feature2 = "percent.chloro") +
geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
geom_hline(yintercept = 30, linetype = "dashed", color = "blue") +
geom_hline(yintercept = 50, linetype = "dashed", color = "green") +
annotate("text", x = Inf, y = 10, label = "10%", color = "red", hjust= 1.1, vjust = -0.5, fontface = "bold") +
annotate("text", x = Inf, y = 30, label = "30%", color = "blue", hjust= 1.1, vjust = -0.5, fontface = "bold") +
annotate("text", x = Inf, y = 50, label = "50%", color = "green", hjust = 1.1, vjust = -0.5, fontface = "bold") +
ggtitle(" ")

p <- p1 | p2

ggsave("ck_fcr_FeatureScatter_mito_chloro.pdf", p, width = 12, height= 8)

ck_fcr <- subset(ck_fcr, subset =
nFeature_RNA > 500 &
nFeature_RNA < 10000 &
nCount_RNA < 60000 &
percent.mito < 5 &
percent.chloro < 20
)

p <- VlnPlot(ck_fcr,
features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.chloro"),
ncol = 4,
pt.size = 0.1)

ggsave("ck_fcr_filter_QC_VlnPlot.pdf", p, width = 12, height = 8)
```

## Source PDF page 10

```text
#fcr filter
fcr.data = Read10X("fcr_matrix_soupx")
fcr <- CreateSeuratObject(counts = fcr.data, project = "fcr", min.cells = 3, min.features = 200)
rownames(fcr[["RNA"]]) <- sub("^gene:", "", rownames(fcr[["RNA"]]))
keep_genes <- grepl("^TraesCS[1-8][ABD]", rownames(fcr[["RNA"]]))
fcr <- subset(fcr, features = rownames(fcr[["RNA"]])[keep_genes])

gene_data <- FetchData(fcr, vars = "TraesCS8B02G990504")
expressed_cells <- sum(gene_data > 0)
total_cells <- nrow(gene_data)
print(paste("           :", expressed_cells))
print(paste("   :", total_cells))
print(paste("   :", round(expressed_cells / total_cells * 100, 2),"%"))
summary(gene_data)

gene_data <- FetchData(fcr, vars = "TraesCS8D02G020913")
expressed_cells <- sum(gene_data > 0)
total_cells <- nrow(gene_data)
print(paste("           :", expressed_cells))
print(paste("   :", total_cells))
print(paste("   :", round(expressed_cells / total_cells * 100, 2),"%"))
summary(gene_data)

fcr[["percent.mito"]] <- PercentageFeatureSet(fcr, features = "TraesCS8B02G990504")
fcr[["percent.chloro"]] <- PercentageFeatureSet(fcr, features = "TraesCS8D02G020913")
head(fcr@meta.data)
p <- VlnPlot(fcr,
features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.chloro"),
ncol = 4,
pt.size = 0.1)

ggsave("fcr_raw_QC_VlnPlot.pdf", p, width = 12, height = 8)

p1 <- FeatureScatter(fcr, feature1 = "nCount_RNA", feature2 = "percent.mito") +
geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
geom_hline(yintercept = 10, linetype = "dashed", color = "blue") +
```

## Source PDF page 11

```text
geom_hline(yintercept = 15, linetype = "dashed", color = "green") +
annotate("text", x = Inf, y = 5, label = "5%", color = "red", hjust =1.1, vjust = -0.5, fontface = "bold") +
annotate("text", x = Inf, y = 10, label = "10%", color = "blue", hjust= 1.1, vjust = -0.5, fontface = "bold") +
annotate("text", x = Inf, y = 15, label = "15%", color = "green", hjust = 1.1, vjust = -0.5, fontface = "bold") +
ggtitle("percent.mito vs nCount_RNA")

p2 <- FeatureScatter(fcr, feature1 = "nCount_RNA", feature2 = "percent.chloro") +
geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
geom_hline(yintercept = 30, linetype = "dashed", color = "blue") +
geom_hline(yintercept = 50, linetype = "dashed", color = "green") +
annotate("text", x = Inf, y = 10, label = "10%", color = "red", hjust= 1.1, vjust = -0.5, fontface = "bold") +
annotate("text", x = Inf, y = 30, label = "30%", color = "blue", hjust= 1.1, vjust = -0.5, fontface = "bold") +
annotate("text", x = Inf, y = 50, label = "50%", color = "green", hjust = 1.1, vjust = -0.5, fontface = "bold") +
ggtitle("percent.chloro vs nCount_RNA")

p <- p1 | p2

ggsave("fcr_FeatureScatter_mito_chloro.pdf", p, width = 12, height = 8)

fcr <- subset(fcr, subset =
nFeature_RNA > 500 &
nFeature_RNA < 10000 &
nCount_RNA < 60000 &
percent.mito < 5 &
percent.chloro < 20
)

p <- VlnPlot(fcr,
features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.chloro"),
ncol = 4,
pt.size = 0.1)

ggsave("fcr_filter_QC_VlnPlot.pdf", p, width = 12, height = 8)
```

## Source PDF page 14

```text
set.seed(54)
#ck_fcr normalization
ck_fcr <- NormalizeData(ck_fcr, normalization.method = "LogNormalize",scale.factor = 10000)
#fcr normalization
fcr <- NormalizeData(fcr, normalization.method = "LogNormalize", scale.factor = 10000)
save.image()
vim findgenes_ck_fcr.R
```

## Source PDF page 15

```text
load(".RData")
library(remotes)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(batchelor)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(clustree)
library(plotly)
library(rtracklayer)
library(tidyHeatmap)
library(future)
library(DoubletFinder)
library(SoupX)
library(glmGamPoi)
library(ggrepel)
library(tidydr)
library(RColorBrewer)
library(AUCell)
library(GSEABase)
library(pheatmap)
library(tibble)
library(clusterProfiler)
library(edgeR)
library(org.Ta.eg.db)
library(igraph)
library(ComplexHeatmap)
library(circlize)
ck_fcr <- FindVariableFeatures(ck_fcr, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ck_fcr)
ck_fcr <- ScaleData(ck_fcr, features = all.genes)
save.image()
nohup Rscript findgenes_ck_fcr.R &
```

## Source PDF page 16

```text
vim findgenes_fcr.R
load(".RData")
library(remotes)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(batchelor)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(clustree)
library(plotly)
library(rtracklayer)
library(tidyHeatmap)
library(future)
library(DoubletFinder)
library(SoupX)
library(glmGamPoi)
library(ggrepel)
library(tidydr)
library(RColorBrewer)
library(AUCell)
library(GSEABase)
library(pheatmap)
library(tibble)
library(clusterProfiler)
library(edgeR)
library(org.Ta.eg.db)
library(igraph)
library(ComplexHeatmap)
library(circlize)
fcr <- FindVariableFeatures(fcr, selection.method = "vst", nfeatures =2000)
all.genes <- rownames(fcr)
fcr <- ScaleData(fcr, features = all.genes)
save.image()
```

## Source PDF page 17

```text
nohup Rscript findgenes_fcr.R &
```

## Source PDF page 18

```text
#ck_fcr dobletfinder
ck_fcr <- SCTransform(ck_fcr, vst.flavor = "v2")
ck_fcr <- RunPCA(ck_fcr, features = VariableFeatures(object = ck_fcr),npcs = 100)
ck_fcr <- FindNeighbors(ck_fcr, dims = 1:50)
ck_fcr <- FindClusters(ck_fcr, resolution = 0.6)
ck_fcr <- RunUMAP(ck_fcr, dims = 1:50)
sweep.res.list <- paramSweep(ck_fcr, PCs = 1:50, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
sweep.stats[order(sweep.stats$BCreal),]
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
homotypic.prop <- modelHomotypic(ck_fcr$seurat_clusters)
nExp_poi <- round(0.08 *nrow(ck_fcr@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
ck_fcr <- doubletFinder(ck_fcr, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, sct = TRUE)
df_col <- grep("^DF.classifications", colnames(ck_fcr@meta.data), value = TRUE)
print(df_col)
p <- DimPlot(
ck_fcr,
reduction = "umap",
group.by = "DF.classifications_0.25_0.14_631",
raster = FALSE
) +
theme_dr(
xlength = 0.2,
ylength = 0.2,
arrow = arrow(length = unit(0.2, "inches"), type = "closed")
) +
theme(
panel.grid = element_blank(),
axis.title = element_text(face = 2, hjust = 0.03),
plot.title = element_text(face = "bold", hjust = 0.5)
) +
labs(title = " ")
ggsave("doublet_singlet_ck_fcr.pdf", p, width = 8, height = 8)
ck_fcr <- subset(ck_fcr, subset = DF.classifications_0.25_0.14_631 =="Singlet")

#fcr dobletfinder
```

## Source PDF page 19

```text
fcr <- SCTransform(fcr)
fcr <- RunPCA(fcr, features = VariableFeatures(object = fcr), npcs = 100)
fcr <- FindNeighbors(fcr, dims = 1:50)
fcr <- FindClusters(fcr, resolution = 0.6)
fcr <- RunUMAP(fcr, dims = 1:50)
fcr <- RunTSNE(fcr, dims = 1:50)
sweep.res.list <- paramSweep(fcr, PCs = 1:50, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
sweep.stats[order(sweep.stats$BCreal),]
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
homotypic.prop <- modelHomotypic(fcr$seurat_clusters)
nExp_poi <- round(0.08 *nrow(fcr@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
fcr <- doubletFinder(fcr, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, nExp =nExp_poi.adj, sct = TRUE)
df_col <- grep("^DF.classifications", colnames(fcr@meta.data), value =TRUE)
print(df_col)
p <- DimPlot(
fcr,
reduction = "umap",
group.by = "DF.classifications_0.25_0.3_728",
raster = FALSE
) +
theme_dr(
xlength = 0.2,
ylength = 0.2,
arrow = arrow(length = unit(0.2, "inches"), type = "closed")
) +
theme(
panel.grid = element_blank(),
axis.title = element_text(face = 2, hjust = 0.03),
plot.title = element_text(face = "bold", hjust = 0.5)
) +
labs(title = " ")
ggsave("doublet_singlet_fcr.pdf", p, width = 8, height = 8)
fcr <- subset(fcr, subset = DF.classifications_0.25_0.3_728 == "Singlet")
```
