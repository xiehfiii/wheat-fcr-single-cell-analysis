suppressPackageStartupMessages({
  library(monocle)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript 07_run_BEAM_one_branch.R SELECTED_CDS_RDS BRANCH_POINT OUTPUT_DIRECTORY")
}

set.seed(54)
input_file <- normalizePath(args[[1]], mustWork = TRUE)
branch_point <- as.integer(args[[2]])
outdir <- args[[3]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cds <- readRDS(input_file)
if (!"num_cells_expressed" %in% colnames(fData(cds))) {
  cds <- detectGenes(cds, min_expr = 0.1)
}
tested_genes <- rownames(fData(cds))[fData(cds)$num_cells_expressed >= 10]
if (length(tested_genes) < 100L) stop("Unexpectedly few expressed genes")

message(sprintf("Branch point %d: testing %d genes", branch_point, length(tested_genes)))
t0 <- Sys.time()
beam <- BEAM(
  cds[tested_genes, ],
  branch_point = branch_point,
  progenitor_method = "duplicate",
  relative_expr = TRUE,
  cores = 10,
  verbose = TRUE
)
elapsed_minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
beam$gene_id <- rownames(beam)
beam <- beam %>% arrange(qval, pval)

summary_df <- data.frame(
  branch_point = branch_point,
  tested_genes = length(tested_genes),
  q_lt_1e_2 = sum(beam$qval < 1e-2, na.rm = TRUE),
  q_lt_1e_3 = sum(beam$qval < 1e-3, na.rm = TRUE),
  q_lt_1e_4 = sum(beam$qval < 1e-4, na.rm = TRUE),
  elapsed_minutes = elapsed_minutes
)

write.csv(beam, file.path(outdir, sprintf("BEAM_branch_point_%d_full.csv", branch_point)), row.names = FALSE)
write.csv(data.frame(gene_id = tested_genes),
          file.path(outdir, sprintf("BEAM_branch_point_%d_tested_genes.csv", branch_point)), row.names = FALSE)
write.csv(summary_df,
          file.path(outdir, sprintf("BEAM_branch_point_%d_summary.csv", branch_point)), row.names = FALSE)
saveRDS(beam, file.path(outdir, sprintf("BEAM_branch_point_%d_full.rds", branch_point)))
capture.output(sessionInfo(), file = file.path(outdir, sprintf("sessionInfo_branch_point_%d.txt", branch_point)))
print(summary_df)
