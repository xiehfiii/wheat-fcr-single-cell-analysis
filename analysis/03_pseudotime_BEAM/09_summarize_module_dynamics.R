suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript 09_summarize_module_dynamics.R COMPONENT_RDS OUTPUT_CSV")

x <- readRDS(normalizePath(args[[1]], mustWork = TRUE))
mat <- x$matrix
module <- paste0("M", x$modules[rownames(mat)])
if (anyNA(module)) stop("Module labels do not cover every matrix row")

curve <- data.frame(gene_id = rownames(mat), module = module, mat, check.names = FALSE) %>%
  pivot_longer(-c(gene_id, module), names_to = "matrix_column", values_to = "z") %>%
  group_by(module, matrix_column) %>%
  summarise(mean_z = mean(z), median_z = median(z), .groups = "drop")

column_order <- data.frame(matrix_column = colnames(mat), matrix_index = seq_len(ncol(mat)))
curve <- curve %>% left_join(column_order, by = "matrix_column") %>% arrange(module, matrix_index)
write.csv(curve, args[[2]], row.names = FALSE)

info_file <- sub("\\.csv$", "_structure.txt", args[[2]])
capture.output({
  cat("matrix dimensions:", paste(dim(mat), collapse = " x "), "\n")
  cat("branch labels:", paste(x$branch_labels, collapse = " | "), "\n")
  cat("first columns:\n")
  print(head(colnames(mat), 15))
  cat("last columns:\n")
  print(tail(colnames(mat), 15))
  cat("annotation_col:\n")
  str(x$annotation_col)
  print(x$annotation_col)
}, file = info_file)

print(curve %>% group_by(module) %>% slice_max(mean_z, n = 3, with_ties = FALSE))
