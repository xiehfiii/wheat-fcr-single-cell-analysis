# scTenifoldKnk virtual knockout

PDF transcription only. The PDF may wrap code or omit glyphs;
use the curated scripts for executable analysis.

## Source PDF page 122

```text
options(timeout = 600)
devtools::install_github("cailab-tamu/scTenifoldKnk")
library(dplyr)
library(Seurat)
library(scTenifoldKnk)
library(ggplot2)
library(ggrepel)
merged <- readRDS("merged.rds")
merged_sub <- subset(merged, idents = c("Epidermis cells I"))
merged_sub
#An object of class Seurat
#134222 features across 1246 samples within 2 assays
#Active assay: RNA (68949 features, 2000 variable features)
# 3 layers present: scale.data, data, counts
# 1 other assay present: SCT
counts <- GetAssayData(merged_sub, layer = "counts", assay = "RNA")
percent_expressed <- rowSums(counts > 0) / ncol(counts)
keep_genes <- names(which(percent_expressed > 0.1))
merged_sub <- subset(merged_sub, features = keep_genes)
merged_sub
#An object of class Seurat
#16464 features across 1246 samples within 2 assays
#Active assay: RNA (8232 features, 370 variable features)
# 3 layers present: scale.data, data, counts
# 1 other assay present: SCT
# 2 dimensional reductions calculated: pca, umap
countMatrix <- GetAssayData(merged_sub, layer = "counts", assay = "RNA")
dim(countMatrix)
#[1] 8232 1246
"TraesCS4A02G001300" %in% rownames(countMatrix)
#[1] TRUE
"TraesCS2B02G369100" %in% rownames(countMatrix)
#[1] TRUE

result <- scTenifoldKnk(countMatrix = countMatrix,
gKO = 'TraesCS4A02G001300',
qc = FALSE,
qc_mtThreshold = 0.1,
qc_minLSize = 1000,
nc_nNet = 10,
nc_nCells = 500,
nc_nComp = 3
```

## Source PDF page 123

```text
)

df_save <- result$diffRegulation
df_save$Gene <- rownames(df_save)
df_save <- df_save %>% dplyr::select(Gene, everything())
write.csv(df_save, file = "Gene_knockout_result_diffRegulation.csv", row.names = FALSE)

ordered_data <- result$diffRegulation[order(-result$diffRegulation$FC), ]
top_genes <- ordered_data[2:21, ]
p <- ggplot(top_genes, aes(x=reorder(gene, FC), y=FC)) +
geom_bar(stat='identity', fill='steelblue') +
coord_flip() +
labs(title=" ",
x="Gene", y="FC") +
theme_minimal()
ggsave("gene_knock_out_barplot.pdf", p, width = 10, height = 10)

df <- result$diffRegulation
df <- df[df$gene != "TraesCS4A02G001300", ]
df$log_pval <- -log10(df$p.adj)
label_genes <- subset(df, abs(Z) > 2 & p.adj < 0.01)

p <- ggplot(df, aes(x=Z, y=log_pval)) +
geom_point(alpha=0.5) +
geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
geom_vline(xintercept=c(2), linetype="dashed", color="blue") +
geom_text_repel(data=label_genes, aes(label=gene),
size=2, max.overlaps=20) +
labs(title=" ",
x="Z-score", y="-log10(p-value)") +
theme_classic()
ggsave("gene_knock_out_dotplot.pdf", p, width = 10, height = 10)



################GO
library(clusterProfiler)
library(org.Ta.eg.db)
library(ggplot2)
library(dplyr)
library(stringr)
df <- result$diffRegulation
target_data <- subset(df, abs(Z) > 2 & p.adj < 0.01 & df$gene != "TraesCS4A02G001300")
```

## Source PDF page 124

```text
deg_genes <- as.character(target_data$gene)
if (length(deg_genes) == 0) {
stop("                            ")
}
print(paste(">>>                :", length(deg_genes)))
out_dir <- "Gene_Knockout_GO_Results"
dir.create(out_dir, showWarnings = FALSE)
COLS <- c("#CECEE8", "#FFE7AD", "#90BFA1")
dot_cols <- c("#61AACF", "#98CADD", "#EAEFF6", "#F9EFEF", "#E9C6C6", "#DA9599")
display_number_bar <- c(15, 15, 15)
display_number_dot <- 21
print(">>>      GO       (BP, CC, MF)...")

run_enrich <- function(ont_type) {
enrichGO(gene = deg_genes,
OrgDb = org.Ta.eg.db,
keyType = "GENEID",
ont = ont_type,
pAdjustMethod = "BH",
pvalueCutoff = 0.05,
qvalueCutoff = 0.05)
}

ego_BP <- run_enrich("BP")
ego_CC <- run_enrich("CC")
ego_MF <- run_enrich("MF")

res_BP <- if(is.null(ego_BP)) data.frame() else as.data.frame(ego_BP)
res_CC <- if(is.null(ego_CC)) data.frame() else as.data.frame(ego_CC)
res_MF <- if(is.null(ego_MF)) data.frame() else as.data.frame(ego_MF)

if(nrow(res_BP) > 0) res_BP$Ontology <- "BP"
if(nrow(res_CC) > 0) res_CC$Ontology <- "CC"
if(nrow(res_MF) > 0) res_MF$Ontology <- "MF"

#
print(paste0(" [  ] BP=", nrow(res_BP), ", CC=", nrow(res_CC), ", MF=", nrow(res_MF)))
df_all_GO <- rbind(res_BP, res_CC, res_MF)

if (nrow(df_all_GO) > 0) {
df_all_GO <- df_all_GO %>% dplyr::select(Ontology, everything())
write.csv(df_all_GO, file = file.path(out_dir, "Knockout_GO_All_Results.csv"), row.names = FALSE)
print(">>>        GO      CSV")
```

## Source PDF page 125

```text
} else {
stop("             GO            ")
}

df_BP_bar <- if(nrow(res_BP) > 0) res_BP %>% slice_head(n = display_number_bar[1]) else NULL
df_CC_bar <- if(nrow(res_CC) > 0) res_CC %>% slice_head(n = display_number_bar[2]) else NULL
df_MF_bar <- if(nrow(res_MF) > 0) res_MF %>% slice_head(n = display_number_bar[3]) else NULL

go_bar_df <- bind_rows(
if(!is.null(df_BP_bar)) df_BP_bar %>% mutate(type = "biological process") else NULL,
if(!is.null(df_CC_bar)) df_CC_bar %>% mutate(type = "cellular component") else NULL,
if(!is.null(df_MF_bar)) df_MF_bar %>% mutate(type = "molecular function") else NULL
)

if (nrow(go_bar_df) > 0) {
go_bar_df$Description <- sapply(go_bar_df$Description, function(x) paste(head(strsplit(x, " ")[[1]], 5), collapse = " "))
go_bar_df$Description <- factor(go_bar_df$Description, levels = rev(unique(go_bar_df$Description)))
go_bar_df$type <- factor(go_bar_df$type, levels = c("biological process", "cellular component", "molecular function"))

p_go_bar <- ggplot(go_bar_df, aes(x = Description, y = Count, fill =type)) +
geom_bar(stat = "identity", width = 0.8) +
coord_flip() +
scale_fill_manual(values = COLS) +
theme_bw() +
theme(panel.grid = element_blank(),
plot.title = element_text(hjust = 0.5),
axis.text.y = element_text(size = 10)) +
labs(x = NULL, y = "Gene Number", title = "GO Enrichment Barplot")

ggsave(filename = file.path(out_dir, "Knockout_GO_barplot.pdf"), plot = p_go_bar, width = 12, height = 8)
print(">>>     Barplot")
}

if (nrow(res_BP) > 0) {
df_BP_dot <- res_BP %>% slice_head(n = display_number_dot)
```
