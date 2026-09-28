suppressPackageStartupMessages({
  library(monocle)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 04_parameter_grid_monocle2.R HVG1000_CDS_RDS OUTPUT_DIRECTORY")
}

set.seed(54)
input_file <- normalizePath(args[[1]], mustWork = TRUE)
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
base_cds <- readRDS(input_file)
ordering_genes <- rownames(subset(fData(base_cds), use_for_ordering %in% TRUE))
stopifnot(length(ordering_genes) == 1000L)

settings <- tidyr::crossing(num_dim = 15L, ncenter = c(50L, 55L, 60L, 65L, 70L))
plots <- vector("list", nrow(settings))
metrics <- vector("list", nrow(settings))

for (i in seq_len(nrow(settings))) {
  set.seed(54)
  nd <- settings$num_dim[[i]]
  nc <- settings$ncenter[[i]]
  cds <- setOrderingFilter(base_cds, ordering_genes)
  cds <- reduceDimension(
    cds, max_components = 2, num_dim = nd, ncenter = nc,
    reduction_method = "DDRTree", verbose = FALSE
  )
  cds <- orderCells(cds)
  state_stats <- as.data.frame(pData(cds)) %>%
    count(State, Samples, name = "n") %>%
    complete(State, Samples = c("ck", "fcr"), fill = list(n = 0)) %>%
    group_by(State) %>%
    mutate(total = sum(n), fraction = n / total) %>%
    ungroup()
  root <- state_stats %>%
    filter(Samples == "ck", total > 10) %>%
    arrange(desc(fraction), desc(total)) %>%
    slice(1) %>%
    pull(State) %>%
    as.character()
  cds <- orderCells(cds, root_state = as.numeric(root))
  meta <- as.data.frame(pData(cds))
  bp <- length(cds@auxOrderingData[["DDRTree"]]$branch_points)
  metrics[[i]] <- data.frame(
    num_dim = nd, ncenter = nc,
    states = length(unique(meta$State)), branch_points = bp,
    root_state = root,
    root_ck_fraction = max(state_stats$fraction[state_stats$Samples == "ck"]),
    median_ck = median(meta$Pseudotime[meta$Samples == "ck"]),
    median_fcr = median(meta$Pseudotime[meta$Samples == "fcr"]),
    rho_chloro = cor(meta$Pseudotime, meta$percent.chloro, method = "spearman"),
    rho_mito = cor(meta$Pseudotime, meta$percent.mito, method = "spearman")
  )
  plots[[i]] <- plot_cell_trajectory(cds, color_by = "Samples", cell_size = 0.55) +
    scale_color_manual(values = c(ck = "#4C78A8", fcr = "#D07A5F")) +
    labs(title = sprintf("PCs %d | centers %d | branch points %d", nd, nc, bp)) +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 7.5, face = "bold"),
      axis.title = element_text(size = 6.5), axis.text = element_text(size = 5.5),
      legend.position = "none", panel.grid = element_blank()
    )
}

metrics_df <- bind_rows(metrics)
write.csv(metrics_df, file.path(outdir, "parameter_grid_metrics.csv"), row.names = FALSE)

comparison <- wrap_plots(plots, ncol = 5)
ragg::agg_png(
  file.path(outdir, "parameter_grid_preview.png"),
  width = 183 / 25.4, height = 62 / 25.4, units = "in", res = 300
)
print(comparison)
dev.off()

print(metrics_df)
