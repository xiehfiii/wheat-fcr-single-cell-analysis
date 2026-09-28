# Pseudobulk contrast and enrichment

PDF transcription only. The PDF may wrap code or omit glyphs;
use the curated scripts for executable analysis.

## Source PDF page 53

```text
ggsave(
"merged_umap_celltype.pdf",
p,
width = 12,
height = 10
)
```

## Source PDF page 54

```text
merged@meta.data$celltype = Idents(merged)
table(merged$celltype)
sample_table = as.data.frame(table(merged@meta.data$orig.ident,
merged@meta.data$celltype))
names(sample_table ) = c("Samples","celltype","CellNumber")
sample_table$Samples = factor(sample_table$Samples,levels = c("ck_fcr","fcr"))
colors = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E78AC3", "#E5C494", "#A6D854")

#       ,
p1 = ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype)) +
geom_bar(position="fill",width = 0.7, linewidth = 0.5,colour = '#222222') +
scale_fill_manual(values = colors) +
theme(panel.grid = element_blank(),
panel.background = element_rect(fill = "transparent",colour = NA),
axis.line.x = element_line(colour = "black"),
axis.line.y = element_line(colour = "black"),
plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size = 15)) +
labs(y="Percentage") +
RotatedAxis()

ggsave("merged_sample_cells_ratio.pdf", p1, width = 8, height = 8)

#####
names(sample_table) = c("celltype","Samples","CellNumber")
sample_colors = c("#FF6F61", "#6B5B95")
p2 = ggplot(sample_table,aes(x=celltype,weight=CellNumber,fill=Samples))+
geom_bar(position="fill",width = 0.7,size = 0.5,colour = '#222222')+
scale_fill_manual(values = sample_colors) +
theme(panel.grid = element_blank(),
panel.background = element_rect(fill = "transparent",colour = NA),
axis.line.x = element_line(colour = "black"),
axis.line.y = element_line(colour = "black"),
plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)) +
```

## Source PDF page 55

```text
labs(y="Percentage")+
RotatedAxis() +
coord_flip()
ggsave("merged_cells_ratio_cells_type.pdf", p2, width = 8, height = 8)
```

## Source PDF page 57

```text
cluster_order <- levels(merged)
my_cols <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
"#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
"#767171", "#136952", "#E18727", "#20854E", "#7876B1")

cell_type_colors <- my_cols[1:length(cluster_order)]
names(cell_type_colors) <- cluster_order
col_fun = colorRamp2(c(0, 0.5, 1), c("#4575B4", "#FFFFBF", "#D73027"))
av1 <- AggregateExpression(merged, group.by = c("orig.ident", "celltype"), assays = "RNA", slot = "counts", return.seurat = FALSE)
av1 <- as.data.frame(av1[[1]])
cg1 <- names(tail(sort(apply(log(av1+1), 1, sd)), 1000))
df1 <- cor(as.matrix(log(av1[cg1,]+1)))
target_rows <- paste0("fcr_", cluster_order)
target_cols <- paste0("ck-fcr_", cluster_order)
use_rows <- target_rows[target_rows %in% rownames(df1)]
use_cols <- target_cols[target_cols %in% colnames(df1)]
df1_clean <- df1[use_rows, use_cols]
df1_clean[is.infinite(df1_clean)] <- 0
df1_clean[is.na(df1_clean)] <- 0
col_anno_data <- gsub("^ck-fcr_", "", colnames(df1_clean))
ha_top = HeatmapAnnotation(
CellType = factor(col_anno_data, levels = cluster_order),
col = list(CellType = cell_type_colors),
show_legend = TRUE,
annotation_name_side = "left"
)
row_anno_data <- gsub("^fcr_", "", rownames(df1_clean))
ha_left = HeatmapAnnotation(
CellType = factor(row_anno_data, levels = cluster_order),
col = list(CellType = cell_type_colors),
show_legend = FALSE,
show_annotation_name = FALSE,
which = "row"
)

pdf("heatmap_1_split_Complex.pdf", width = 13, height = 9)
ht1 <- Heatmap(df1_clean,
name = "Correlation",
col = col_fun,
cluster_rows = FALSE,
cluster_columns = FALSE,
show_row_names = FALSE,
```

## Source PDF page 58

```text
show_column_names = FALSE,
top_annotation = ha_top,
left_annotation = ha_left,
rect_gp = gpar(col = "white", lwd = 1),
row_title = "fcr",
column_title = "ck_fcr"
)

draw(ht1, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

#####
library(ComplexHeatmap)
library(circlize)

# 1.
cluster_order <- levels(merged)
my_cols <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
"#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
"#767171", "#136952", "#E18727", "#20854E", "#7876B1")
cell_type_colors <- my_cols[1:length(cluster_order)]
names(cell_type_colors) <- cluster_order

#
col_fun = colorRamp2(c(0, 0.5, 1), c("#4575B4", "#FFFFBF", "#D73027"))

# 2.                     celltype
av2 <- AggregateExpression(merged, group.by = "celltype", assays = "RNA", slot = "counts", return.seurat = FALSE)
av2 <- as.data.frame(av2[[1]])

cg2 <- names(tail(sort(apply(log(av2+1), 1, sd)), 1000))
df2 <- cor(as.matrix(log(av2[cg2,]+1)))

#       cluster_order
use_types <- cluster_order[cluster_order %in% rownames(df2)]
df2_clean <- df2[use_types, use_types]
df2_clean[is.infinite(df2_clean)] <- 0
df2_clean[is.na(df2_clean)] <- 0

# 3.           Annotation
#
ha_top_2 = HeatmapAnnotation(
CellType = factor(colnames(df2_clean), levels = cluster_order),
col = list(CellType = cell_type_colors),
```

## Source PDF page 59

```text
show_legend = TRUE,         #
show_annotation_name = TRUE, #
annotation_name_side = "left", #         "CellType"
annotation_name_rot = 0     #
)

#
ha_left_2 = HeatmapAnnotation(
CellType = factor(rownames(df2_clean), levels = cluster_order),
col = list(CellType = cell_type_colors),
show_legend = FALSE,
show_annotation_name = FALSE, #
which = "row"
)

# 4.                            CellType
pdf("heatmap_2_merged_Complex.pdf", width = 11, height = 7)

ht2 <- Heatmap(df2_clean,
name = "Correlation",
col = col_fun,

cluster_rows = FALSE,
cluster_columns = FALSE,

#
show_row_names = FALSE,
show_column_names = FALSE,

#
top_annotation = ha_top_2,
left_annotation = ha_left_2,

#
rect_gp = gpar(col = "white", lwd = 1),

#
row_title = "",
column_title = ""
)

#                     "CellType"                 padding
draw(ht2,
heatmap_legend_side = "right",
```

## Source PDF page 60

```text
annotation_legend_side = "right",
padding = unit(c(2, 8, 2, 2), "mm") # , , ,
)

dev.off()
```

## Source PDF page 61

```text
dir.create("All_Clusters_DEG_Results", showWarnings = FALSE)

merged$Samples <- ifelse(grepl("ck", merged$orig.ident, ignore.case = TRUE), "ck", "fcr")

all_clusters <- levels(merged)

for (cluster_name in all_clusters) {

print(paste(">>> Processing cluster:", cluster_name))

safe_name <- str_replace_all(cluster_name, " ", "_")

tryCatch({

sub_seurat <- subset(merged, idents = cluster_name)

sample_counts <- table(sub_seurat$Samples)
if (length(sample_counts) < 2) {
print(paste(" [Skip] Insufficient samples for:", cluster_name))
next
}

agg_data <- AggregateExpression(sub_seurat,
group.by = "Samples",
assays = "RNA",
slot = "counts",
return.seurat = FALSE)
counts_matrix <- agg_data$RNA

y <- DGEList(counts = counts_matrix, group = colnames(counts_matrix))
keep <- rowSums(cpm(y) > 1) >= 1
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

bcv <- 0.2
dispersion_val <- bcv ^ 2

et <- exactTest(y, dispersion = dispersion_val, pair = c("ck", "fcr"))

top_genes <- topTags(et, n = Inf)$table
```

## Source PDF page 62

```text
DEG_edgeR <- top_genes %>%
rownames_to_column("Symbol") %>%
rename(padj = FDR)

fc_cutoff <- 1
p_cutoff <- 0.05

DEG_edgeR <- DEG_edgeR %>%
mutate(Type = case_when(
padj < p_cutoff & logFC >= fc_cutoff ~ "up",
padj < p_cutoff & logFC <= -fc_cutoff ~ "down",
TRUE ~ "stable"
)) %>%
arrange(desc(abs(logFC)))

top_labels <- head(DEG_edgeR[DEG_edgeR$Type != "stable", ], 10)$Symbol

p <- ggplot(DEG_edgeR, aes(x = logFC, y = -log10(padj))) +
geom_point(aes(color = Type), size = 2.5, alpha = 0.8) +
scale_color_manual(values = c("down"="#00468B", "stable"="gray","up"="#E64B35")) +
geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = 2, color = 'black', linewidth = 0.5) +
geom_hline(yintercept = -log10(p_cutoff), linetype = 2, color ='black', linewidth = 0.5) +

labs(x = "Log2(Fold Change)",
y = "-Log10(FDR)",
title = " ") +

theme_bw() +
theme(
panel.grid = element_blank(),
plot.title = element_blank(),
legend.position = "top"
) +
geom_text_repel(data = subset(DEG_edgeR, Symbol %in% top_labels),
aes(label = Symbol),
box.padding = 0.5,
max.overlaps = Inf)

ggsave(filename = paste0("All_Clusters_DEG_Results/", safe_name, "_Volcano.pdf"),
plot = p, width = 8, height = 8)

```

## Source PDF page 63

```text
sig_genes <- subset(DEG_edgeR, Type != "stable")
write.csv(sig_genes,
file = paste0("All_Clusters_DEG_Results/", safe_name, "_Significant_DEGs.csv"),
row.names = FALSE)

up_genes <- subset(DEG_edgeR, Type == "up")
write.csv(up_genes,
file = paste0("All_Clusters_DEG_Results/", safe_name, "_Up_Genes.csv"),
row.names = FALSE)

down_genes <- subset(DEG_edgeR, Type == "down")
write.csv(down_genes,
file = paste0("All_Clusters_DEG_Results/", safe_name, "_Down_Genes.csv"),
row.names = FALSE)

print(paste(" [Done] Saved results for:", cluster_name))

}, error = function(e) {
print(paste(" [Error]", cluster_name, ":", e$message))
})
}
```

## Source PDF page 64

```text
library(clusterProfiler)
library(org.Ta.eg.db)
library(ggplot2)
library(dplyr)
library(stringr)

# 1.
dir.create("All_Clusters_GO_Results", showWarnings = FALSE)
deg_dir <- "All_Clusters_DEG_Results"

# 2.
COLS <- c("#CECEE8", "#FFE7AD", "#90BFA1")
dot_cols <- c("#61AACF", "#98CADD", "#EAEFF6", "#F9EFEF", "#E9C6C6", "#DA9599")

display_number_bar <- c(15, 15, 15)
display_number_dot <- 21

directions <- c("up", "down")
all_clusters <- levels(merged) #  merged

# 3.
for (cluster_name in all_clusters) {

print(paste(">>>      :", cluster_name))
safe_name <- str_replace_all(cluster_name, " ", "_")

for (dir in directions) {

file_name <- paste0(safe_name, "_", tools::toTitleCase(dir), "_Genes.csv")
file_path <- file.path(deg_dir, file_name)

if (!file.exists(file_path)) {
print(paste(" [  ]        :", file_name))
next
}

deg_data <- read.csv(file_path)

if (nrow(deg_data) == 0) {
print(paste(" [  ]          :", dir))
next
```

## Source PDF page 65

```text
}

deg_genes <- as.character(deg_data$Symbol)
deg_genes <- deg_genes[deg_genes != "" & !is.na(deg_genes)]

print(paste("       GO       :", dir, "(   :", length(deg_genes), ")"))

tryCatch({
ego_BP <- enrichGO(gene = deg_genes, OrgDb = org.Ta.eg.db, keyType = "GENEID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ego_CC <- enrichGO(gene = deg_genes, OrgDb = org.Ta.eg.db, keyType = "GENEID", ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ego_MF <- enrichGO(gene = deg_genes, OrgDb = org.Ta.eg.db, keyType = "GENEID", ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

res_BP <- if(is.null(ego_BP)) data.frame() else as.data.frame(ego_BP)
res_CC <- if(is.null(ego_CC)) data.frame() else as.data.frame(ego_CC)
res_MF <- if(is.null(ego_MF)) data.frame() else as.data.frame(ego_MF)

count_BP <- nrow(res_BP)
count_CC <- nrow(res_CC)
count_MF <- nrow(res_MF)
total_count <- count_BP + count_CC + count_MF

print(paste0(" [  ] ", cluster_name, " (", dir, "): BP=", count_BP, " , CC=", count_CC, " , MF=", count_MF, " , =", total_count, " "))

if(count_BP > 0) res_BP$Ontology <- "BP"
if(count_CC > 0) res_CC$Ontology <- "CC"
if(count_MF > 0) res_MF$Ontology <- "MF"

df_all_GO <- rbind(res_BP, res_CC, res_MF)

if (nrow(df_all_GO) > 0) {
# *****           dplyr::select     BiocGenerics::select*****
df_all_GO <- df_all_GO %>% dplyr::select(Ontology, everything())
```

## Source PDF page 66

```text
write.csv(df_all_GO,
file = file.path("All_Clusters_GO_Results", paste0(safe_name, "_", dir, "_GO_All_Results.csv")),
row.names = FALSE)
}

df_BP_bar <- if(count_BP > 0) res_BP %>% slice_head(n = display_number_bar[1]) else NULL
df_CC_bar <- if(count_CC > 0) res_CC %>% slice_head(n = display_number_bar[2]) else NULL
df_MF_bar <- if(count_MF > 0) res_MF %>% slice_head(n = display_number_bar[3]) else NULL

go_bar_df <- bind_rows(
if(!is.null(df_BP_bar)) df_BP_bar %>% mutate(type = "biological process") else NULL,
if(!is.null(df_CC_bar)) df_CC_bar %>% mutate(type = "cellular component") else NULL,
if(!is.null(df_MF_bar)) df_MF_bar %>% mutate(type = "molecularfunction") else NULL
)

if (nrow(go_bar_df) > 0) {
go_bar_df$Description <- sapply(go_bar_df$Description, function(x) paste(head(strsplit(x, " ")[[1]], 5), collapse = " "))
go_bar_df$Description <- factor(go_bar_df$Description, levels= rev(unique(go_bar_df$Description)))

p_go_bar <- ggplot(go_bar_df, aes(x = Description, y = Count, fill = type)) +
geom_bar(stat = "identity", width = 0.8) +
coord_flip() +
scale_fill_manual(values = COLS) +
theme_bw() +
theme(panel.grid = element_blank(), plot.title = element_blank()) +
labs(x = NULL, y = "Gene Number")

ggsave(filename = file.path("All_Clusters_GO_Results", paste0(safe_name, "_", dir, "_GO_barplot.pdf")), plot = p_go_bar, width = 12, height = 8)
}

if (count_BP > 0) {
df_BP_dot <- res_BP %>% slice_head(n = display_number_dot)

```

## Source PDF page 67

```text
df_BP_dot$GeneRatio_Num <- sapply(strsplit(df_BP_dot$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
df_BP_dot$Description <- sapply(df_BP_dot$Description, function(x) paste(head(strsplit(x, " ")[[1]], 5), collapse = " "))

p_go_dot <- ggplot(df_BP_dot, aes(x = GeneRatio_Num, y = reorder(Description, GeneRatio_Num), size = Count, color = p.adjust)) +
geom_point() +
scale_color_gradientn(colors = dot_cols, trans = "reverse", name = "adj.P") +
scale_size(range = c(3, 8)) +
theme_bw() +
theme(panel.grid = element_blank(), plot.title = element_blank()) +
labs(x = "Gene Ratio", y = NULL)

ggsave(filename = file.path("All_Clusters_GO_Results", paste0(safe_name, "_", dir, "_GO_BP_dotplot.pdf")), plot = p_go_dot, width = 10, height = 8)
}

print(paste(" [  ]", dir, "         "))

}, error = function(e) {
print(paste(" [  ]", cluster_name, ":", e$message))
})
}
}
```

## Source PDF page 68

```text
dir.create("All_Clusters_Overlap_Results", showWarnings = FALSE)
deg_dir <- "All_Clusters_DEG_Results"

if (!exists("merged.markers")) stop("           merged.markers")

all_clusters <- levels(merged)
directions <- c("up", "down")

# 2.      (       )
for (cluster_name in all_clusters) {

safe_name <- str_replace_all(cluster_name, " ", "_")
message(paste0("\n========================================================"))
message(paste0(">>>      : ", cluster_name))

# --- Step 1:      Marker    (       ) ---
#            merged.markers$cluster         "Mesophyll cells"
#     as.character                      (factor)
cluster_marker_df <- subset(merged.markers, as.character(cluster) ==as.character(cluster_name))
marker_genes_list <- cluster_marker_df$gene
n_marker <- length(marker_genes_list)

if (n_marker == 0) {
message(paste0(" [  ]   merged.markers         '", cluster_name, "'     "))
message(" (     merged.markers$clusterID)")
next
}

# --- Step 2:                   ---

#     Up
up_file <- file.path(deg_dir, paste0(safe_name, "_Up_Genes.csv"))
if (file.exists(up_file)) {
up_df <- read.csv(up_file)
n_deg_up <- nrow(up_df)
up_genes <- if(n_deg_up > 0) up_df$Symbol else character(0)
} else {
```

## Source PDF page 69

```text
n_deg_up <- 0
up_genes <- character(0)
}

#     Down
down_file <- file.path(deg_dir, paste0(safe_name, "_Down_Genes.csv"))
if (file.exists(down_file)) {
down_df <- read.csv(down_file)
n_deg_down <- nrow(down_df)
down_genes <- if(n_deg_down > 0) down_df$Symbol else character(0)
} else {
n_deg_down <- 0
down_genes <- character(0)
}

#
n_overlap_up <- length(intersect(up_genes, marker_genes_list))
n_overlap_down <- length(intersect(down_genes, marker_genes_list))
n_deg_total <- n_deg_up + n_deg_down

# ---                ---
cat(paste0("          :\n"))
cat(paste0("    - Marker    : ", n_marker, "\n"))
cat(paste0("    -         :   ", n_deg_total, " (Up: ", n_deg_up, " | Down: ", n_deg_down, ")\n"))🔥
cat(paste0("    -   Marker ∩ Up: ", n_overlap_up, "\n"))
cat(paste0("    -  Marker ∩ Down: ", n_overlap_down, "\n"))❄
cat(paste0(" --------------------------------------\n"))


# --- Step 3:             ---
process_list <- list(
list(dir="up", df=if(exists("up_df")) up_df else NULL, genes=up_genes),
list(dir="down", df=if(exists("down_df")) down_df else NULL, genes=down_genes)
)

for (item in process_list) {
dir <- item$dir
deg_df <- item$df
deg_genes_list <- item$genes

if (length(deg_genes_list) == 0) next

#
```

## Source PDF page 70

```text
overlap_genes <- intersect(deg_genes_list, marker_genes_list)

if (length(overlap_genes) > 0) {
result_df <- data.frame(Gene = overlap_genes)

#    LogFC   FDR
deg_idx <- match(overlap_genes, deg_df$Symbol)
result_df$Induced_LogFC <- deg_df$logFC[deg_idx]
result_df$FDR <- deg_df$padj[deg_idx]

#    Marker LogFC
marker_idx <- match(overlap_genes, cluster_marker_df$gene)
result_df$Marker_LogFC <- cluster_marker_df$avg_log2FC[marker_idx]

#
if (dir == "up") {
result_df <- result_df %>% arrange(desc(Induced_LogFC))
} else {
result_df <- result_df %>% arrange(Induced_LogFC)
}

#
out_file_name <- paste0(safe_name, "_", tools::toTitleCase(dir),"_Overlap.csv")
out_file <- file.path("All_Clusters_Overlap_Results", out_file_name)

write.csv(result_df, out_file, row.names = FALSE)
}
}
}
```

## Source PDF page 71

```text
library(clusterProfiler)
library(org.Ta.eg.db)
library(ggplot2)
library(dplyr)
library(stringr)

# 1.
#
dir.create("All_Clusters_Overlap_GO_Results", showWarnings = FALSE)
#        (                )
input_dir <- "All_Clusters_Overlap_Results"

# 2.      (      )
COLS <- c("#CECEE8", "#FFE7AD", "#90BFA1")
dot_cols <- c("#61AACF", "#98CADD", "#EAEFF6", "#F9EFEF", "#E9C6C6", "#DA9599")

display_number_bar <- c(15, 15, 15)
display_number_dot <- 21

directions <- c("up", "down")
all_clusters <- levels(merged) #  merged

# 3.
for (cluster_name in all_clusters) {

print(paste(">>>      :", cluster_name))
safe_name <- str_replace_all(cluster_name, " ", "_")

for (dir in directions) {

# ---     A:          ---
#        : ClusterName_Up_Overlap.csv
file_name <- paste0(safe_name, "_", tools::toTitleCase(dir), "_Overlap.csv")
file_path <- file.path(input_dir, file_name)

if (!file.exists(file_path)) {
print(paste(" [  ]        :", file_name))
next
}

deg_data <- read.csv(file_path)
```

## Source PDF page 72

```text

if (nrow(deg_data) == 0) {
print(paste(" [  ]          :", dir))
next
}

# ---     B:        Gene ---
#                "Gene"
if (!"Gene" %in% colnames(deg_data)) {
print(paste(" [  ]         'Gene'  :", file_name))
next
}

deg_genes <- as.character(deg_data$Gene)
deg_genes <- deg_genes[deg_genes != "" & !is.na(deg_genes)]

print(paste("       GO       :", dir, "(   :", length(deg_genes), ")"))

tryCatch({
# --- 1.         ---
ego_BP <- enrichGO(gene = deg_genes, OrgDb = org.Ta.eg.db, keyType = "GENEID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ego_CC <- enrichGO(gene = deg_genes, OrgDb = org.Ta.eg.db, keyType = "GENEID", ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
ego_MF <- enrichGO(gene = deg_genes, OrgDb = org.Ta.eg.db, keyType = "GENEID", ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

# --- 2.      ---
res_BP <- if(is.null(ego_BP)) data.frame() else as.data.frame(ego_BP)
res_CC <- if(is.null(ego_CC)) data.frame() else as.data.frame(ego_CC)
res_MF <- if(is.null(ego_MF)) data.frame() else as.data.frame(ego_MF)

count_BP <- nrow(res_BP)
count_CC <- nrow(res_CC)
count_MF <- nrow(res_MF)
total_count <- count_BP + count_CC + count_MF

print(paste0(" [  ] ", cluster_name, " (", dir, "): BP=", count_BP, " , CC=", count_CC, " , MF=", count_MF, " , =", total_coun
```

## Source PDF page 73

```text

# --- 3.   CSV ---
if(count_BP > 0) res_BP$Ontology <- "BP"
if(count_CC > 0) res_CC$Ontology <- "CC"
if(count_MF > 0) res_MF$Ontology <- "MF"

df_all_GO <- rbind(res_BP, res_CC, res_MF)

if (nrow(df_all_GO) > 0) {
# *****             dplyr::select    BiocGenerics::select*****
df_all_GO <- df_all_GO %>% dplyr::select(Ontology, everything())

#         up/down
write.csv(df_all_GO,
file = file.path("All_Clusters_Overlap_GO_Results", paste0(safe_name, "_", dir, "_Overlap_GO_All_Results.csv")),
row.names = FALSE)
}

# --- 4.  : Barplot ---
df_BP_bar <- if(count_BP > 0) res_BP %>% slice_head(n = display_number_bar[1]) else NULL
df_CC_bar <- if(count_CC > 0) res_CC %>% slice_head(n = display_number_bar[2]) else NULL
df_MF_bar <- if(count_MF > 0) res_MF %>% slice_head(n = display_number_bar[3]) else NULL

go_bar_df <- bind_rows(
if(!is.null(df_BP_bar)) df_BP_bar %>% mutate(type = "biological process") else NULL,
if(!is.null(df_CC_bar)) df_CC_bar %>% mutate(type = "cellular component") else NULL,
if(!is.null(df_MF_bar)) df_MF_bar %>% mutate(type = "molecularfunction") else NULL
)

if (nrow(go_bar_df) > 0) {
go_bar_df$Description <- sapply(go_bar_df$Description, function(x) paste(head(strsplit(x, " ")[[1]], 5), collapse = " "))
go_bar_df$Description <- factor(go_bar_df$Description, levels= rev(unique(go_bar_df$Description)))

p_go_bar <- ggplot(go_bar_df, aes(x = Description, y = Count, f
```

## Source PDF page 74

```text
geom_bar(stat = "identity", width = 0.8) +
coord_flip() +
scale_fill_manual(values = COLS) +
theme_bw() +
theme(panel.grid = element_blank(), plot.title = element_blank()) +
labs(x = NULL, y = "Gene Number")

ggsave(filename = file.path("All_Clusters_Overlap_GO_Results",paste0(safe_name, "_", dir, "_Overlap_GO_barplot.pdf")), plot = p_go_bar, width = 12, height = 8)
}

# --- 5.  : Dotplot ( BP) ---
if (count_BP > 0) {
df_BP_dot <- res_BP %>% slice_head(n = display_number_dot)

df_BP_dot$GeneRatio_Num <- sapply(strsplit(df_BP_dot$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
df_BP_dot$Description <- sapply(df_BP_dot$Description, function(x) paste(head(strsplit(x, " ")[[1]], 5), collapse = " "))

p_go_dot <- ggplot(df_BP_dot, aes(x = GeneRatio_Num, y = reorder(Description, GeneRatio_Num), size = Count, color = p.adjust)) +
geom_point() +
scale_color_gradientn(colors = dot_cols, trans = "reverse", name = "adj.P") +
scale_size(range = c(3, 8)) +
theme_bw() +
theme(panel.grid = element_blank(), plot.title = element_blank()) +
labs(x = "Gene Ratio", y = NULL)

ggsave(filename = file.path("All_Clusters_Overlap_GO_Results",paste0(safe_name, "_", dir, "_Overlap_GO_BP_dotplot.pdf")), plot = p_go_dot, width = 10, height = 8)
}

print(paste(" [  ]", dir, "         "))

}, error = function(e) {
print(paste(" [  ]", cluster_name, ":", e$message))
})
}
}
```

## Source PDF page 76

```text
library(clusterProfiler)
library(ggplot2)
library(biomaRt)
library(enrichplot)
library(stringr)
library(dplyr)

# ==============================================================================
# 1.
# ==============================================================================
options(timeout = 600)
message(">>>  KEGG    ")
gene_data <- read.table("up.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(gene_data)[1] <- "Symbol"
my_genes <- gene_data$Symbol
message(paste0("              : ", length(my_genes)))

# ==============================================================================
# 2. ID   (BiomaRt)
# ==============================================================================
message(">>> [1/4]  Ensembl Plants  ID    ...")

ensembl_plants <- useMart(biomart = "plants_mart",
dataset = "taestivum_eg_gene",
host = "https://plants.ensembl.org")

#    Entrez ID
gene_map <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id'),
filters = 'ensembl_gene_id',
values = my_genes,
mart = ensembl_plants)

#
gene_map <- gene_map[!is.na(gene_map$entrezgene_id), ]
gene_data <- merge(gene_data, gene_map, by.x = "Symbol", by.y = "ensembl_gene_id")
gene_data <- gene_data[grepl("^\\d+$", gene_data$entrezgene_id), ]
message(paste0("   ID                 : ", nrow(gene_data)))
```

## Source PDF page 77

```text
if (nrow(gene_data) == 0) stop("     ID          Entrez IDKEGG     ")

# ==============================================================================
# 3.   ID      (Entrez -> Symbol)
# ==============================================================================
id_lookup <- setNames(gene_data$Symbol, gene_data$entrezgene_id)

# ==============================================================================
# 4.   KEGG     (     )
# ==============================================================================
message(">>> [2/4]  KEGG      ...")
target_entrez <- gene_data$entrezgene_id
kk <- tryCatch({
enrichKEGG(gene = target_entrez,
organism = 'taes',
pvalueCutoff = 0.05,
qvalueCutoff = 0.05)
}, error = function(e) return(NULL))

if (is.null(kk) || nrow(kk) == 0) {
stop("              KEGG     ")
}

# ==============================================================================
# 5.
# ==============================================================================
message(">>> [3/4]    ...")

#
kk@result$Description <- gsub(" - Triticum aestivum \\(bread wheat\\)", "", kk@result$Description)

# ID     (Entrez -> Symbol)
kk_readable <- kk
kk_readable@result$geneID <- sapply(kk@result$geneID, function(x) {
ids <- unlist(strsplit(x, "/"))
symbols <- id_lookup[ids]
symbols <- symbols[!is.na(symbols)]
paste(symbols, collapse = "/")
```

## Source PDF page 78

```text
})

message(paste0("            : ", nrow(kk)))

# ==============================================================================
# 6.
# ==============================================================================
message(">>> [4/4]    ...")

# ---     ---
p1 <- dotplot(kk_readable, showCategory = 15) +
ggtitle(" ") +
theme(plot.title = element_text(hjust = 0.5))

ggsave("up_KEGG_Dotplot.pdf", p1, width = 10, height = 8)

#
write.csv(as.data.frame(kk_readable), "up_KEGG_Results.csv", row.names= FALSE)

message("                    ")
target_file <- "All_Clusters_Overlap_Results/Mesophyll_cells_Up_Overlap.csv"

if (file.exists(target_file)) {
cluster_data <- read.csv(target_file)
print(">>> Mesophyll_cells Top 10        ")
print(head(cluster_data, 10))
} else {
print("                   ")
}
```

## Source PDF page 79

```text
genes_file <- "All_Clusters_Overlap_Results/Mesophyll_cells_Up_Overlap.csv"
go_result_file <- "All_Clusters_GO_Results/Mesophyll_cells_up_GO_All_Results.csv"
if (!file.exists(genes_file)) stop("         ")
gene_df <- read.csv(genes_file)
target_col <- if("Gene" %in% colnames(gene_df)) "Gene" else "Symbol"
top10_genes <- head(gene_df[[target_col]], 10)
if (!file.exists(go_result_file)) stop(" GO")
go_df <- read.csv(go_result_file)
message(paste0(">>>    Top 10      GO         ...\n"))

for (gene in top10_genes) {
cat(paste0("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
cat(paste0("     ID: ", gene, "\n"))
matched_rows <- go_df[grep(gene, go_df$geneID), ]
if (nrow(matched_rows) > 0) {
cat(paste0("           ", nrow(matched_rows), "\n"))
out_info <- matched_rows[, c("Ontology", "ID", "Description", "p.adjust")]
out_info <- out_info[order(out_info$p.adjust), ]
print(head(out_info, 10))
if (nrow(out_info) > 10) {
cat(paste0(" ... (      ", nrow(out_info) - 10, "  )\n"))
}
} else {
cat("                  GO     \n")
cat("  (                   geneID      )\n")
}
cat("\n")
}
```
