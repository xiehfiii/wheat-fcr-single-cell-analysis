suppressPackageStartupMessages({
  library(monocle)
  library(dplyr)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: Rscript 02_audit_current_monocle2.R",
    "EPIDERMAL_RDS DEG_CSV OUTPUT_DIRECTORY ROOT_STATE"
  ))
}

set.seed(54)
input_file <- normalizePath(args[[1]], mustWork = TRUE)
deg_file <- normalizePath(args[[2]], mustWork = TRUE)
outdir <- args[[3]]
root_state_requested <- as.numeric(args[[4]])
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

input_list <- readRDS(input_file)
counts <- input_list$data
pd_df <- input_list$pd
fd_df <- input_list$fData

pd <- new("AnnotatedDataFrame", data = pd_df)
fd <- new("AnnotatedDataFrame", data = fd_df)
cds <- newCellDataSet(
  counts,
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size()
)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, cores = 10, relative_expr = TRUE)

deg_df <- read.csv(deg_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"Symbol" %in% names(deg_df)) stop("DEG table has no Symbol column.")
ordering_genes_input <- unique(na.omit(as.character(deg_df$Symbol)))
ordering_genes <- intersect(ordering_genes_input, rownames(cds))
ordering_genes_expressed <- ordering_genes[
  Matrix::rowSums(exprs(cds)[ordering_genes, , drop = FALSE] > 0) >= 10
]

cds <- setOrderingFilter(cds, ordering_genes_expressed)
cds <- reduceDimension(
  cds,
  max_components = 2,
  num_dim = 15,
  ncenter = 80,
  reduction_method = "DDRTree",
  verbose = FALSE
)
cds <- orderCells(cds)

state_stats_initial <- pData(cds) %>%
  as.data.frame() %>%
  count(State, Samples, name = "n") %>%
  tidyr::complete(State, Samples = c("ck", "fcr"), fill = list(n = 0)) %>%
  group_by(State) %>%
  mutate(total = sum(n), fraction = n / total) %>%
  ungroup()
write.csv(
  state_stats_initial,
  file.path(outdir, "state_by_sample_before_reroot.csv"), row.names = FALSE
)

root_summary <- state_stats_initial %>%
  filter(Samples == "ck", total > 10) %>%
  arrange(desc(fraction), desc(total))
recommended_root <- root_summary$State[[1]]

if (!root_state_requested %in% unique(pData(cds)$State)) {
  stop("Requested root state is absent from the reconstructed trajectory.")
}
cds <- orderCells(cds, root_state = root_state_requested)

meta <- pData(cds) %>% as.data.frame()
meta$cell_id <- rownames(meta)
write.csv(meta, file.path(outdir, "cell_metadata_with_pseudotime.csv"), row.names = FALSE)

state_sample <- meta %>%
  count(State, Samples, name = "n") %>%
  tidyr::complete(State, Samples = c("ck", "fcr"), fill = list(n = 0)) %>%
  group_by(State) %>%
  mutate(total = sum(n), fraction = n / total) %>%
  ungroup()
write.csv(state_sample, file.path(outdir, "state_by_sample_after_reroot.csv"), row.names = FALSE)

numeric_covariates <- intersect(
  c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.chloro"),
  names(meta)
)
cor_summary <- data.frame(
  covariate = numeric_covariates,
  spearman_rho = vapply(
    numeric_covariates,
    function(x) cor(meta$Pseudotime, meta[[x]], method = "spearman", use = "complete.obs"),
    numeric(1)
  )
)
write.csv(cor_summary, file.path(outdir, "pseudotime_qc_correlations.csv"), row.names = FALSE)

target_gene <- "TraesCS4A02G001300"
target_summary <- data.frame()
if (target_gene %in% rownames(cds)) {
  value <- as.numeric(exprs(cds)[target_gene, ])
  target_summary <- data.frame(
    gene = target_gene,
    cells = length(value),
    expressing_cells = sum(value > 0),
    fraction_expressing = mean(value > 0),
    mean_expression = mean(value),
    median_expression = median(value),
    spearman_pseudotime = cor(value, meta$Pseudotime, method = "spearman")
  )
  write.csv(target_summary, file.path(outdir, "target_gene_summary.csv"), row.names = FALSE)
}

summary <- data.frame(
  metric = c(
    "cells", "genes", "ordering_genes_input", "ordering_genes_in_matrix",
    "ordering_genes_expressed_in_at_least_10_cells", "states",
    "recommended_root_state_from_CK_fraction", "requested_root_state",
    "root_state_matches_recommendation", "branch_points"
  ),
  value = c(
    ncol(cds), nrow(cds), length(ordering_genes_input), length(ordering_genes),
    length(ordering_genes_expressed), length(unique(meta$State)),
    recommended_root, root_state_requested,
    identical(as.numeric(recommended_root), root_state_requested),
    length(cds@auxOrderingData[["DDRTree"]]$branch_points)
  )
)
write.csv(summary, file.path(outdir, "analysis_audit_summary.csv"), row.names = FALSE)
saveRDS(cds, file.path(outdir, "mycds_reconstructed_current_definition.rds"), compress = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))

print(summary)
print(state_sample)
print(cor_summary)
print(target_summary)
