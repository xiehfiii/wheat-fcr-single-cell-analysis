suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else Sys.getenv("SCRNASEQ_PROJECT_DIR", unset = ".")
out_dir <- if (length(args) >= 2) args[[2]] else file.path(project_dir, "submission_source_20260922")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsvgz <- function(x, file) {
  con <- gzfile(file, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

merged <- readRDS(file.path(project_dir, "merged.rds"))
DefaultAssay(merged) <- "RNA"

cell_id <- colnames(merged)
md <- merged@meta.data
um <- Embeddings(merged, "umap")
stopifnot(identical(rownames(md), rownames(um)), identical(rownames(md), cell_id))

get_expr <- function(gene) {
  if (!gene %in% rownames(merged)) return(rep(NA_real_, ncol(merged)))
  as.numeric(GetAssayData(merged, assay = "RNA", layer = "data")[gene, ])
}

cell_source <- data.frame(
  cell_id = cell_id,
  sample = as.character(md$Samples),
  original_library = as.character(md$orig.ident),
  cluster = as.character(md$seurat_clusters),
  cell_type = as.character(md$celltype_final),
  annotation_basis = as.character(md$annotation_basis),
  nCount_RNA = md$nCount_RNA,
  nFeature_RNA = md$nFeature_RNA,
  percent_mito = md$percent.mito,
  percent_chloro = md$percent.chloro,
  UMAP_1 = um[, 1],
  UMAP_2 = um[, 2],
  TaERF87_expression = get_expr("TraesCS4A02G001300"),
  TaMGBP1_expression = get_expr("TraesCS3D02G094200"),
  stringsAsFactors = FALSE
)
write_tsvgz(cell_source, file.path(out_dir, "main_Fig1_Fig2_cell_metadata_UMAP_expression.tsv.gz"))

cluster_key <- unique(cell_source[, c("cluster", "cell_type", "annotation_basis")])
cluster_key <- cluster_key[order(as.integer(cluster_key$cluster)), ]
write.csv(cluster_key, file.path(out_dir, "main_Fig1_cluster_annotation_key.csv"), row.names = FALSE, quote = TRUE)

cell_counts <- as.data.frame(table(sample = cell_source$sample, cluster = cell_source$cluster,
                                   cell_type = cell_source$cell_type), stringsAsFactors = FALSE)
cell_counts <- cell_counts[cell_counts$Freq > 0, ]
cell_counts$total_by_sample <- ave(cell_counts$Freq, cell_counts$sample, FUN = sum)
cell_counts$proportion <- cell_counts$Freq / cell_counts$total_by_sample
write.csv(cell_counts, file.path(out_dir, "main_Fig1_cell_counts_and_proportions.csv"), row.names = FALSE, quote = TRUE)

row_map <- data.frame(
  row_pos = 1:15,
  cluster = c("8", "12", "2", "14", "1", "6", "13", "7", "10", "3", "9", "0", "5", "11", "4"),
  stringsAsFactors = FALSE
)
block_map <- data.frame(
  block_order = 1:13,
  marker_block = c(
    "Epidermal I/II lineage", "Epidermal III identity and subtype",
    "Guard-cell identity", "Mesophyll I/II lineage", "Procambial identity",
    "Vascular-parenchyma evidence", "Phloem-parenchyma identity",
    "Cell-wall-remodelling state I", "Cell-wall-remodelling state II",
    "Defense-detoxification state", "Defense-JA-responsive state",
    "Proliferating-S-phase state", "Translation-active state"
  ),
  clusters_in_block = c("8,12", "2", "14", "1,6", "13", "7", "10", "3", "9", "0", "5", "11", "4"),
  stringsAsFactors = FALSE
)
genes_by_block <- list(
  c("TraesCS5B02G145900", "TraesCS3D02G140300", "TraesCS4A02G068400"),
  c("TraesCS4A02G357400", "TraesCS1A02G334100", "TraesCS3D02G293300"),
  c("TraesCS4A02G235600", "TraesCS4B02G079300", "TraesCS2D02G302800"),
  c("TraesCS2A02G187200", "TraesCS2A02G252600", "TraesCS5A02G482800"),
  c("TraesCS6A02G036100", "TraesCS6B02G050700", "TraesCS6D02G041700"),
  c("TraesCS4B02G243500", "TraesCS6A02G404500", "TraesCS7B02G213000"),
  c("TraesCS7A02G261100", "TraesCS7B02G160000", "TraesCS7D02G263100"),
  c("TraesCS1A02G001100", "TraesCS1A02G000900", "TraesCS1A02G000400"),
  c("TraesCS5D02G459400", "TraesCS6B02G207500", "TraesCS1A02G157200"),
  c("TraesCS5A02G236700", "TraesCS5D02G243700", "TraesCS4A02G202200"),
  c("TraesCS1D02G034800", "TraesCS1B02G042200", "TraesCS2A02G395000"),
  c("TraesCS6A02G034500", "TraesCS6A02G086200", "TraesCS1D02G039800"),
  c("TraesCS4A02G245100", "TraesCS5B02G088100", "TraesCS4D02G069100")
)
feature_map <- do.call(rbind, lapply(seq_along(genes_by_block), function(i) {
  data.frame(block_order = i, marker_block = block_map$marker_block[i],
             gene = genes_by_block[[i]], stringsAsFactors = FALSE)
}))
feature_map$gene_order <- seq_len(nrow(feature_map))
dot_raw <- DotPlot(merged, features = feature_map$gene, assay = "RNA",
                   group.by = "seurat_clusters", scale = TRUE,
                   col.min = -1, col.max = 2.5)$data
dot_raw$cluster <- as.character(dot_raw$id)
dot_raw$gene <- as.character(dot_raw$features.plot)
dot_source <- merge(dot_raw[, c("cluster", "gene", "pct.exp", "avg.exp", "avg.exp.scaled")],
                    feature_map, by = "gene", all.x = TRUE, sort = FALSE)
dot_source <- merge(dot_source, row_map, by = "cluster", all.x = TRUE, sort = FALSE)
dot_source <- dot_source[order(dot_source$gene_order, dot_source$row_pos), ]
names(dot_source)[names(dot_source) == "pct.exp"] <- "percent_expressed"
names(dot_source)[names(dot_source) == "avg.exp"] <- "average_expression_unscaled"
names(dot_source)[names(dot_source) == "avg.exp.scaled"] <- "average_expression_scaled"
write.csv(dot_source, file.path(out_dir, "main_Fig1C_marker_dotplot_source.csv"), row.names = FALSE, quote = TRUE)
write.csv(block_map, file.path(out_dir, "main_Fig1C_marker_blocks.csv"), row.names = FALSE, quote = TRUE)

# Post-filter QC values used in the violin/scatter panels.
qc_filtered <- cell_source[, c("cell_id", "sample", "nCount_RNA", "nFeature_RNA", "percent_mito", "percent_chloro")]
write_tsvgz(qc_filtered, file.path(out_dir, "supp_Fig1_filtered_cell_QC.tsv.gz"))

# Raw FCR object retained on the server.
fcr <- readRDS(file.path(project_dir, "fcr.rds"))
raw_fcr <- cbind(cell_id = rownames(fcr@meta.data), sample = "FCR", fcr@meta.data)
write_tsvgz(raw_fcr, file.path(out_dir, "supp_Fig1_raw_FCR_cell_QC.tsv.gz"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo_source_export.txt"))
cat("WROTE", out_dir, "\n")
