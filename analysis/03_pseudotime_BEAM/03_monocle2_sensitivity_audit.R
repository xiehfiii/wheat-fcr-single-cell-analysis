suppressPackageStartupMessages({
  library(monocle)
  library(dplyr)
  library(tidyr)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript 03_monocle2_sensitivity_audit.R EPIDERMAL_RDS DEG_CSV OUTPUT_DIRECTORY")
}

set.seed(54)
input_file <- normalizePath(args[[1]], mustWork = TRUE)
deg_file <- normalizePath(args[[2]], mustWork = TRUE)
outdir <- args[[3]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

input_list <- readRDS(input_file)
pd <- new("AnnotatedDataFrame", data = input_list$pd)
fd <- new("AnnotatedDataFrame", data = input_list$fData)
base_cds <- newCellDataSet(
  input_list$data,
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size()
)
base_cds <- estimateSizeFactors(base_cds)
base_cds <- estimateDispersions(base_cds, cores = 10, relative_expr = TRUE)

deg <- read.csv(deg_file, stringsAsFactors = FALSE)
deg_genes <- unique(intersect(as.character(deg$Symbol), rownames(base_cds)))
detected_in <- Matrix::rowSums(exprs(base_cds) > 0)
deg_genes <- deg_genes[detected_in[deg_genes] >= 10]

disp <- dispersionTable(base_cds) %>%
  filter(
    gene_id %in% rownames(base_cds),
    mean_expression >= 0.1,
    dispersion_fit > 0,
    dispersion_empirical >= dispersion_fit,
    detected_in[gene_id] >= 10
  ) %>%
  mutate(dispersion_ratio = dispersion_empirical / dispersion_fit) %>%
  arrange(desc(dispersion_ratio), desc(mean_expression))

hvg_1000 <- head(unique(disp$gene_id), 1000)
hvg_2000 <- head(unique(disp$gene_id), 2000)
write.csv(disp, file.path(outdir, "unsupervised_dispersion_candidates.csv"), row.names = FALSE)

rank_auc <- function(value, group) {
  x <- value[group == "fcr"]
  y <- value[group == "ck"]
  as.numeric(wilcox.test(x, y, exact = FALSE)$statistic) / (length(x) * length(y))
}

run_trajectory <- function(label, ordering_genes, num_dim, ncenter) {
  set.seed(54)
  cds <- base_cds
  cds <- setOrderingFilter(cds, ordering_genes)
  cds <- reduceDimension(
    cds, max_components = 2, num_dim = num_dim, ncenter = ncenter,
    reduction_method = "DDRTree", verbose = FALSE
  )
  cds <- orderCells(cds)
  meta0 <- as.data.frame(pData(cds))
  state_stats <- meta0 %>%
    count(State, Samples, name = "n") %>%
    complete(State, Samples = c("ck", "fcr"), fill = list(n = 0)) %>%
    group_by(State) %>%
    mutate(total = sum(n), fraction = n / total) %>%
    ungroup()
  root_state <- state_stats %>%
    filter(Samples == "ck", total > 10) %>%
    arrange(desc(fraction), desc(total)) %>%
    slice(1) %>%
    pull(State) %>%
    as.character()
  cds <- orderCells(cds, root_state = as.numeric(root_state))
  meta <- as.data.frame(pData(cds))
  meta$cell_id <- rownames(meta)
  target <- "TraesCS4A02G001300"
  target_expr <- as.numeric(exprs(cds)[target, ])
  metrics <- data.frame(
    run = label,
    ordering_genes = length(ordering_genes),
    num_dim = num_dim,
    ncenter = ncenter,
    states = length(unique(meta$State)),
    branch_points = length(cds@auxOrderingData[["DDRTree"]]$branch_points),
    root_state = root_state,
    root_ck_fraction = max(state_stats$fraction[state_stats$Samples == "ck"]),
    median_pseudotime_ck = median(meta$Pseudotime[meta$Samples == "ck"]),
    median_pseudotime_fcr = median(meta$Pseudotime[meta$Samples == "fcr"]),
    fcr_vs_ck_pseudotime_auc = rank_auc(meta$Pseudotime, meta$Samples),
    rho_nCount_RNA = cor(meta$Pseudotime, meta$nCount_RNA, method = "spearman"),
    rho_nFeature_RNA = cor(meta$Pseudotime, meta$nFeature_RNA, method = "spearman"),
    rho_percent_mito = cor(meta$Pseudotime, meta$percent.mito, method = "spearman"),
    rho_percent_chloro = cor(meta$Pseudotime, meta$percent.chloro, method = "spearman"),
    target_gene_rho_pseudotime = cor(target_expr, meta$Pseudotime, method = "spearman")
  )
  state_stats$run <- label
  write.csv(meta, file.path(outdir, paste0(label, "_cell_metadata.csv")), row.names = FALSE)
  saveRDS(cds, file.path(outdir, paste0(label, "_cds.rds")), compress = FALSE)
  list(metrics = metrics, state = state_stats)
}

runs <- list(
  treatment_DEG = list(genes = deg_genes, num_dim = 15, ncenter = 80),
  unsupervised_HVG2000 = list(genes = hvg_2000, num_dim = 15, ncenter = 80),
  unsupervised_HVG1000 = list(genes = hvg_1000, num_dim = 15, ncenter = 80),
  unsupervised_HVG2000_dim10 = list(genes = hvg_2000, num_dim = 10, ncenter = 80)
)

results <- lapply(names(runs), function(nm) {
  x <- runs[[nm]]
  run_trajectory(nm, x$genes, x$num_dim, x$ncenter)
})

metrics <- bind_rows(lapply(results, `[[`, "metrics"))
states <- bind_rows(lapply(results, `[[`, "state"))
write.csv(metrics, file.path(outdir, "trajectory_sensitivity_metrics.csv"), row.names = FALSE)
write.csv(states, file.path(outdir, "trajectory_sensitivity_state_composition.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
print(metrics)
print(states)
