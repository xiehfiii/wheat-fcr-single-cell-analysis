suppressPackageStartupMessages({
  library(monocle)
  library(dplyr)
  library(cluster)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript 08_build_BEAM_heatmap_matrix.R SELECTED_CDS_RDS BEAM_CSV BRANCH_POINT OUTPUT_DIRECTORY")
}

set.seed(54)
cds_file <- normalizePath(args[[1]], mustWork = TRUE)
beam_file <- normalizePath(args[[2]], mustWork = TRUE)
branch_point <- as.integer(args[[3]])
outdir <- args[[4]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cds <- readRDS(cds_file)
beam <- read.csv(beam_file, stringsAsFactors = FALSE, check.names = FALSE)
sig <- beam %>% filter(status == "OK", !is.na(qval), qval < 1e-4) %>% arrange(qval, pval)
genes <- intersect(sig$gene_id, rownames(cds))
if (length(genes) < 50L) stop("Too few BEAM-significant genes")

branch_labels <- if (branch_point == 1L) {
  c("State 4 terminal", "State 5 terminal")
} else {
  c("State 2 terminal", "Other downstream states")
}

# Generate the standardized, smoothed branch-expression matrix once. The
# temporary PDF is retained as an audit rendering, not used as the final panel.
grDevices::cairo_pdf(
  file.path(outdir, sprintf("BEAM_branch_point_%d_audit_heatmap.pdf", branch_point)),
  width = 7.2, height = 9, family = "Arial"
)
hm <- plot_genes_branched_heatmap(
  cds[genes, ],
  branch_point = branch_point,
  branch_labels = branch_labels,
  num_clusters = 4,
  cores = 10,
  use_gene_short_name = FALSE,
  show_rownames = FALSE,
  return_heatmap = TRUE
)
dev.off()

mat <- hm$heatmap_matrix
tree <- hm$ph_res$tree_row
if (is.null(tree)) stop("The heatmap row tree was not returned")

distance <- dist(mat)
k_grid <- 2:min(8L, nrow(mat) - 1L)
k_audit <- lapply(k_grid, function(k) {
  cl <- cutree(tree, k = k)
  sil <- silhouette(cl, distance)
  tab <- table(cl)
  data.frame(
    branch_point = branch_point,
    k = k,
    mean_silhouette = mean(sil[, "sil_width"]),
    min_module_n = min(tab),
    min_module_fraction = min(tab) / length(cl),
    max_module_fraction = max(tab) / length(cl)
  )
}) %>% bind_rows()

# Prefer the highest silhouette among solutions without tiny modules. If no
# solution has every module >=5%, fall back to the highest silhouette overall.
eligible <- k_audit %>% filter(min_module_fraction >= 0.05)
if (nrow(eligible) == 0L) eligible <- k_audit
selected_k <- eligible %>% arrange(desc(mean_silhouette), k) %>% slice(1) %>% pull(k)
modules <- cutree(tree, k = selected_k)

module_df <- data.frame(
  gene_id = names(modules),
  module = paste0("M", modules),
  stringsAsFactors = FALSE
) %>%
  left_join(sig %>% select(gene_id, pval, qval, num_cells_expressed), by = "gene_id")

module_summary <- module_df %>%
  count(module, name = "genes") %>%
  mutate(branch_point = branch_point, selected_k = selected_k, .before = 1)

write.csv(k_audit, file.path(outdir, sprintf("BEAM_branch_point_%d_module_k_audit.csv", branch_point)), row.names = FALSE)
write.csv(module_df, file.path(outdir, sprintf("BEAM_branch_point_%d_gene_modules.csv", branch_point)), row.names = FALSE)
write.csv(module_summary, file.path(outdir, sprintf("BEAM_branch_point_%d_module_summary.csv", branch_point)), row.names = FALSE)
write.csv(data.frame(gene_id = rownames(mat), mat, check.names = FALSE),
          file.path(outdir, sprintf("BEAM_branch_point_%d_smoothed_z_matrix.csv", branch_point)), row.names = FALSE)
saveRDS(list(matrix = mat, tree = tree, selected_k = selected_k,
             modules = modules, branch_labels = branch_labels,
             annotation_col = hm$annotation_col),
        file.path(outdir, sprintf("BEAM_branch_point_%d_heatmap_components.rds", branch_point)))
capture.output(sessionInfo(), file = file.path(outdir, sprintf("sessionInfo_heatmap_branch_point_%d.txt", branch_point)))
print(k_audit)
print(module_summary)
