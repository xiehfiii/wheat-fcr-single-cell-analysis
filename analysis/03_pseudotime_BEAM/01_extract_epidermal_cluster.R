suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 01_extract_epidermal_cluster.R MERGED_RDS OUTPUT_RDS")
}

set.seed(54)
merged_file <- normalizePath(args[[1]], mustWork = TRUE)
output_file <- args[[2]]

merged <- readRDS(merged_file)
keep <- as.character(merged$seurat_clusters) == "8"
if (!any(keep)) stop("Cluster 8 is absent from merged.")
scRNAsub <- subset(merged, cells = colnames(merged)[keep])
counts <- GetAssayData(scRNAsub, assay = "RNA", layer = "counts")
counts <- as(counts, "sparseMatrix")
pd <- scRNAsub@meta.data
fData <- data.frame(
  gene_short_name = rownames(counts),
  row.names = rownames(counts),
  stringsAsFactors = FALSE
)

stopifnot(identical(colnames(counts), rownames(pd)))
stopifnot(all(as.character(pd$seurat_clusters) == "8"))

saveRDS(
  list(data = counts, pd = pd, fData = fData),
  output_file,
  compress = FALSE
)

cat("cells\t", ncol(counts), "\n", sep = "")
cat("genes\t", nrow(counts), "\n", sep = "")
cat("nonzero_genes\t", sum(Matrix::rowSums(counts > 0) > 0), "\n", sep = "")
cat("CK_cells\t", sum(pd$Samples == "ck"), "\n", sep = "")
cat("FCR_cells\t", sum(pd$Samples == "fcr"), "\n", sep = "")
cat("final_celltype\t", paste(unique(as.character(pd$celltype)), collapse = ";"), "\n", sep = "")
