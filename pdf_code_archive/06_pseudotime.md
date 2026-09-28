# Monocle2 pseudotime and BEAM

PDF transcription only. The PDF may wrap code or omit glyphs;
use the curated scripts for executable analysis.

## Source PDF page 87

```text
nohup python3 scvelo_pipeline.py --loom /home/xiehangfei/scRNAseq/merged_15000cells/velocyto/merged_15000cells.loom --outdir scvelo_result --n_top_genes 2000 --plot --umap seurat_umap_coords.csv --clusters seurat_clusters.csv &
# scrna
id <- c("8")
scRNAsub <- subset(merged, seurat_clusters %in% id)
data <- GetAssayData(scRNAsub, assay = "RNA", layer = "counts")
data <- as(data, 'sparseMatrix')
pd <- scRNAsub@meta.data
fData <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
saveRDS(list(data = data, pd = pd, fData = fData), "Epidermis_cells.rds")
```

## Source PDF page 88

```text
#####scRNAseq R
epidermis_markers <- subset(merged.markers, cluster == "Epidermis cells I")
epidermis_genes <- epidermis_markers$gene
write.csv(epidermis_genes, "Epidermis_cells_I_markers.csv",
row.names = FALSE, quote = FALSE)

#####monocle2 R
library(monocle)
library(ggplot2)
library(igraph)
library(gridExtra)
library(dplyr)
input_list <- readRDS("Epidermis_cells.rds")
data <- input_list$data
pd_df <- input_list$pd
fd_df <- input_list$fData
pd <- new('AnnotatedDataFrame', data = pd_df)
fd <- new('AnnotatedDataFrame', data = fd_df)
mycds <- newCellDataSet(data,
phenoData = pd,
featureData = fd,
expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores = 10, relative_expr = TRUE)
deg_df <- read.csv("Epidermis_cells_I_Significant_DEGs.csv",
header = TRUE,
stringsAsFactors = FALSE)
ordering_genes <- deg_df$Symbol
mycds <- setOrderingFilter(mycds, ordering_genes)

p <- plot_pc_variance_explained(mycds, return_all = F)
ggsave("elbowplot_epidermis_pseudotime.pdf", p, width = 10, height = 8)

mycds <- reduceDimension(mycds,
max_components = 2,
num_dim = 15,
ncenter = 80,
reduction_method = 'DDRTree',
verbose = F)
mycds <- orderCells(mycds)

########select root
```

## Source PDF page 89

```text
meta_data <- pData(mycds)
control_group_name <- "ck"
state_stats <- meta_data %>%
group_by(State) %>%
summarise(
Total_Cells = n(),
CK_Count = sum(Samples == control_group_name),
FCR_Count = sum(Samples != control_group_name),
CK_Percent = round((CK_Count / Total_Cells) * 100, 2)
) %>%
arrange(desc(CK_Percent))
print(">>> State       (CK)            ")
print(as.data.frame(state_stats))
candidate <- state_stats %>% filter(Total_Cells > 10) %>% slice(1)
root_suggestion <- candidate$State
cat(paste0("\n>>>  Root State: ", root_suggestion,
" (CK   : ", candidate$CK_Percent, "%)\n"))

mycds <- orderCells(mycds, root_state = 3)
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
plot2 <- plot_cell_trajectory(mycds, color_by = "Samples")
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
plotc <- plot1 | plot2 | plot3
ggsave("Epidermis_cells_pseudotime.pdf", plotc, width = 18, height = 6)

plot_gene_monocle <- function(cds, gene_id, group_by = NULL, ncol = 2){

if (!gene_id %in% rownames(cds)) {
message(paste("                ", gene_id))
return(NULL)
}

expr_val <- as.numeric(exprs(cds)[gene_id, ])
expr_val[expr_val == 0] <- NA
cds_temp <- cds
pData(cds_temp)[[gene_id]] <- expr_val

do_facet <- !is.null(group_by) && nzchar(group_by)
if (do_facet) {
if (!group_by %in% colnames(pData(cds_temp))) {
message(paste("     pData         ", group_by))
message("    colnames(pData(cds))")
return(NULL)
}
}
```

## Source PDF page 90

```text

p <- plot_cell_trajectory(cds_temp, color_by = gene_id) +
scale_color_gradient(
low = "yellow",
high = "red",
na.value = "grey90"
) +
theme(legend.position = "right")

if (do_facet) {
p <- p + facet_wrap(as.formula(paste("~", group_by)), ncol = ncol)
}

return(p)
}

p <- plot_gene_monocle(mycds, "TraesCS4A02G001300")
ggsave("TraesCS4A02G001300_pseudotime.pdf", p, width = 12, height = 8)

p <- plot_gene_monocle(mycds, "TraesCS4A02G001300", group_by = "Samples")
ggsave("TraesCS4A02G001300_pseudotime_facet.pdf", p, width = 18, height = 8)

#####BEAM
library(pheatmap)
BEAM_res <- BEAM(mycds, branch_point = 1, cores = 4, progenitor_method= "duplicate")
BEAM_genes <- row.names(subset(BEAM_res, qval < 1e-4))
print(paste("   ", length(BEAM_genes), "       "))

pdf("BEAM_Heatmap_epidermis.pdf", width = 8, height = 10)
hm <- plot_genes_branched_heatmap(mycds[BEAM_genes,],
branch_point = 1,
num_clusters = 5,
cores = 1,
use_gene_short_name = FALSE,
show_rownames = FALSE,
return_heatmap = TRUE)
dev.off()

gene_clusters <- cutree(hm$ph_res$tree_row, k = 5)
gene_cluster_df <- data.frame(GeneID = names(gene_clusters),
Cluster = gene_clusters)

my_tf <- "TraesCS4A02G001300"
```
