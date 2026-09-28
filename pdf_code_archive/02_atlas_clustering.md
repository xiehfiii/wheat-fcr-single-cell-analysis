# Atlas integration, clustering and marker discovery

PDF transcription only. The PDF may wrap code or omit glyphs;
use the curated scripts for executable analysis.

## Source PDF page 21

```text
merged <- merge(x = ck_fcr,
y = fcr,
add.cell.ids = c("ck_fcr", "fcr"),
project = "merged")
DefaultAssay(merged) <- "RNA"
merged <- JoinLayers(merged)
merged <- NormalizeData(merged, normalization.method = "LogNormalize",scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merged)
merged <- ScaleData(merged)
merged <- SCTransform(merged, vst.flavor = "v2")
merged <- RunPCA(merged, features = VariableFeatures(object = merged))
merged <- FindNeighbors(merged, dims = 1:50)
merged <- FindClusters(merged, resolution = 0.6)
merged <- RunUMAP(merged, dims = 1:50)

p1 <- DimPlot(merged, reduction = "umap", group.by = "orig.ident") +
ggtitle(" ")
p2 <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
ggtitle(" ")
p <- p1 + p2
ggsave("merged_umap_sample_vs_cluster.pdf", p, width = 16, height = 8)
save.image()
```

## Source PDF page 22

```text
merged <- FindClusters(
object = merged,
resolution = c(seq(0,1.5,.1))
)
mycols <- c(
"#E0A19A","#F2C8D6","#F7EAE5","#BFABCF",
"#BEDEEE","#B1D9D4","#E0DFA8","#FFD193",
"#F8B8C3","#C4BF18","#FFA512","#F45082",
"#63C163","#25A030","#FECDCF","#FE8D80"
)
fmt <- function(x) ifelse(abs(x - round(x)) < 1e-8, as.character(as.integer(round(x))), sprintf("%.1f", x))
res_vals <- seq(0, 1.5, 0.1)
names(mycols) <- sapply(res_vals, fmt)
p <- clustree(merged@meta.data, prefix = "SCT_snn_res.") +
scale_colour_manual(values = mycols, breaks = names(mycols))
ggsave("clustree_merged.pdf", p, width = 10, height = 12)

res_seq <- seq(0.1, 1.5, 0.1)
meta_names <- paste0("SCT_snn_res.", sapply(res_seq, fmt))
meta_names <- intersect(meta_names, colnames(merged@meta.data))
p1 <- lapply(meta_names, function(m) {
DimPlot(merged,
reduction = "umap",
group.by = m,
raster = FALSE,
pt.size = 0.5) +
ggtitle(sub("SCT_snn_res\\.", "cluster_", m)) +
theme(plot.title = element_text(hjust = 0.5),
legend.position = "right")
})

p <- wrap_plots(p1, ncol = 3, nrow = 5)
ggsave("cluster_resolution_merged.pdf", p, width = 15, height = 21)

merged <- FindClusters(merged, resolution = 0.6)
merged <- RunUMAP(merged, dims = 1:50)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(grid)

mycols <- c(
```

## Source PDF page 23

```text
"#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
"#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
"#767171", "#136952", "#E18727", "#20854E", "#7876B1"
)
cluster_levels <- levels(merged$seurat_clusters)
cluster_names <- sapply(cluster_levels, function(x) {
unique(as.character(merged@meta.data$celltype[merged$seurat_clusters== x]))
})

draw_key_circled_num <- function(data, params, size) {
idx <- which(toupper(mycols) == toupper(data$colour)) - 1
if(length(idx) == 0) idx <- ""

grobTree(
circleGrob(r = 0.42, gp = gpar(fill = data$colour, col = NA)),
textGrob(idx, x = 0.5, y = 0.5,
gp = gpar(col = "white", fontsize = 8, fontface = "bold"))
)
}
```

## Source PDF page 24

```text
DefaultAssay(merged) <- "RNA"
merged <- JoinLayers(merged)
merged.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
top100 <- merged.markers %>% group_by(cluster) %>% top_n(100, avg_log2FC) %>% arrange(cluster, desc(avg_log2FC)) %>% group_split()
top10 <- merged.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% arrange(cluster, desc(avg_log2FC)) %>% group_split()
top3 <- merged.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC) %>% arrange(cluster, desc(avg_log2FC)) %>% group_split()
top2 <- merged.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC) %>% arrange(cluster, desc(avg_log2FC)) %>% group_split()
top1 <- merged.markers %>% group_by(cluster) %>% top_n(1, avg_log2FC) %>% arrange(cluster, desc(avg_log2FC)) %>% group_split()
#####vlnplot
top1_genes <- pull(arrange(top_n(group_by(merged.markers, cluster), 1,avg_log2FC), cluster, desc(avg_log2FC)), gene)
p <- VlnPlot(merged, features = top1_genes, group.by = "seurat_clusters", pt.size = 0)
ggsave("VlnPlot_top1_merged.pdf", p, width = 20, height = 15)
#save markers
clusters <- sort(unique(merged.markers$cluster))
outfile <- "cluster0_14_top100_merged.tsv"
con <- file(outfile, open = "w")
for (i in seq_along(clusters)) {
clust <- clusters[i]
df <- top100[[i]]$gene
writeLines(paste0("cluster", clust), con)
writeLines(df, con)
writeLines("", con)
}
close(con)

write.table(merged.markers, "merged_all_markers.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

top10_df <- merged.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% arrange(cluster, desc(avg_log2FC))
write.table(top10_df, "merged_top10_markers.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

## Source PDF page 25

```text
#####dotplot
top3_df <- merged.markers %>%
group_by(cluster) %>%
top_n(3, avg_log2FC) %>%
arrange(cluster, desc(avg_log2FC))
top3_unique <- top3_df$gene
mycols <- c("#8FB4BE", "#AFC9CF", "#D5E1E3", "#EBBFC2", "#E28187", "#D93F49")
p <- DotPlot(merged, features = top3_unique, group.by = "seurat_clusters") +
coord_flip() +
theme(axis.text.y = element_text(size = 8)) +
scale_color_gradientn(colors = mycols)
ggsave("dotplot_top3_merged.pdf", p, width = 15, height = 10)
```

## Source PDF page 27

```text
#####heatmap
merged <- JoinLayers(merged)
top10_df <- merged.markers %>%
group_by(cluster) %>%
top_n(10, avg_log2FC)

top10_sorted <- top10_df %>%
arrange(cluster, desc(avg_log2FC)) %>%
pull(gene)

merged <- ScaleData(merged, features = top10_sorted)

p <- DoHeatmap(
merged,
features = top10_sorted,
group.by = "seurat_clusters",
group.bar = TRUE,
size = 4,
group.colors = c(
"#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
"#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
"#767171", "#136952", "#E18727", "#20854E", "#7876B1"
)
) +
scale_fill_gradientn(
colours = c("#61AACF", "#98CADD", "#EAEFF6", "#F9EFEF", "#E9C6C6","#DA9599"),
na.value = "grey90"
)

ggsave("heatmap_top10_merged.pdf", p, width = 20, height = 15)
```
