suppressPackageStartupMessages({
  library(monocle)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 05_build_selected_pseudotime_preview.R HVG1000_CDS_RDS OUTPUT_DIRECTORY")
}

set.seed(54)
input_file <- normalizePath(args[[1]], mustWork = TRUE)
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
tag_family <- "Times New Roman"
col_ck <- "#4C78A8"
col_fcr <- "#D07A5F"
state_colours <- c("#456A95", "#63A58C", "#D2A24C", "#B46A87", "#7B6FA8")

cds <- readRDS(input_file)
ordering_genes <- rownames(subset(fData(cds), use_for_ordering %in% TRUE))
stopifnot(length(ordering_genes) == 1000L)
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(
  cds,
  max_components = 2,
  num_dim = 10,
  ncenter = 40,
  reduction_method = "DDRTree",
  verbose = FALSE
)
cds <- orderCells(cds)

state_composition <- as.data.frame(pData(cds)) %>%
  count(State, Samples, name = "n") %>%
  complete(State, Samples = c("ck", "fcr"), fill = list(n = 0)) %>%
  group_by(State) %>%
  mutate(total = sum(n), fraction = n / total) %>%
  ungroup()
root_state <- state_composition %>%
  filter(Samples == "ck", total > 10) %>%
  arrange(desc(fraction), desc(total)) %>%
  slice(1) %>%
  pull(State) %>%
  as.character()
cds <- orderCells(cds, root_state = as.numeric(root_state))

target_gene <- "TraesCS4A02G001300"
normalized_target <- log1p(
  as.numeric(exprs(cds)[target_gene, ]) / as.numeric(sizeFactors(cds))
)
pData(cds)$target_expression <- normalized_target
pData(cds)$Samples <- factor(pData(cds)$Samples, levels = c("ck", "fcr"), labels = c("CK", "FCR"))

meta <- as.data.frame(pData(cds))
meta$cell_id <- rownames(meta)
state_composition_final <- meta %>%
  count(State, Samples, name = "n") %>%
  group_by(State) %>%
  mutate(total = sum(n), fraction = n / total) %>%
  ungroup()

rank_auc <- function(value, group) {
  x <- value[group == "FCR"]
  y <- value[group == "CK"]
  as.numeric(wilcox.test(x, y, exact = FALSE)$statistic) / (length(x) * length(y))
}

qa <- data.frame(
  metric = c(
    "cells", "ordering_genes", "num_dim", "ncenter", "states",
    "branch_points", "root_state", "root_CK_fraction",
    "median_pseudotime_CK", "median_pseudotime_FCR",
    "FCR_vs_CK_pseudotime_AUC_descriptive",
    "rho_pseudotime_percent_mito", "rho_pseudotime_percent_chloro",
    "target_fraction_expressing", "target_rho_pseudotime"
  ),
  value = c(
    ncol(cds), length(ordering_genes), 10, 40,
    length(unique(meta$State)),
    length(cds@auxOrderingData[["DDRTree"]]$branch_points),
    root_state,
    max(state_composition$fraction[state_composition$Samples == "ck"], na.rm = TRUE),
    median(meta$Pseudotime[meta$Samples == "CK"]),
    median(meta$Pseudotime[meta$Samples == "FCR"]),
    rank_auc(meta$Pseudotime, meta$Samples),
    cor(meta$Pseudotime, meta$percent.mito, method = "spearman"),
    cor(meta$Pseudotime, meta$percent.chloro, method = "spearman"),
    mean(normalized_target > 0),
    cor(meta$Pseudotime, normalized_target, method = "spearman")
  )
)

trajectory_theme <- theme_classic(base_size = 7.2, base_family = font_family) +
  theme(
    axis.line = element_line(linewidth = 0.32, colour = "#282828"),
    axis.ticks = element_line(linewidth = 0.28, colour = "#282828"),
    axis.text = element_text(size = 5.5, colour = "#282828"),
    axis.title = element_text(size = 6.4, colour = "#282828"),
    plot.title = element_text(size = 7.4, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 6.0),
    legend.text = element_text(size = 5.5),
    legend.key.height = grid::unit(3.0, "mm"),
    legend.position = "top",
    panel.grid = element_blank(),
    plot.margin = margin(2, 3, 2, 3)
  )

p_state <- plot_cell_trajectory(
  cds, color_by = "State", cell_size = 0.62,
  show_branch_points = FALSE, show_backbone = TRUE
) +
  scale_colour_manual(values = state_colours[seq_along(unique(meta$State))]) +
  labs(title = "Transcriptional state") + trajectory_theme

p_sample <- plot_cell_trajectory(
  cds, color_by = "Samples", cell_size = 0.62,
  show_branch_points = FALSE, show_backbone = TRUE
) +
  scale_colour_manual(values = c(CK = col_ck, FCR = col_fcr)) +
  labs(title = "Sample origin", colour = "Sample") + trajectory_theme

p_time <- plot_cell_trajectory(
  cds, color_by = "Pseudotime", cell_size = 0.62,
  show_branch_points = FALSE, show_backbone = TRUE
) +
  scale_colour_gradientn(
    colours = c("#E8EEF3", "#8FB8CF", "#3F7DA8", "#173B63"),
    name = "Pseudotime"
  ) +
  labs(title = "Pseudotime") + trajectory_theme

medians <- meta %>% group_by(Samples) %>% summarise(median = median(Pseudotime), .groups = "drop")
p_density <- ggplot(meta, aes(x = Pseudotime, fill = Samples, colour = Samples)) +
  geom_density(alpha = 0.22, linewidth = 0.65, adjust = 1.05, key_glyph = "path") +
  geom_vline(
    data = medians, aes(xintercept = median, colour = Samples),
    linewidth = 0.55, linetype = "22"
  ) +
  scale_fill_manual(values = c(CK = col_ck, FCR = col_fcr)) +
  scale_colour_manual(values = c(CK = col_ck, FCR = col_fcr)) +
  labs(x = "Pseudotime", y = "Cell density", fill = NULL, colour = NULL) +
  theme_classic(base_size = 7.2, base_family = font_family) +
  theme(
    axis.line = element_line(linewidth = 0.32),
    axis.text = element_text(size = 5.6), axis.title = element_text(size = 6.4),
    legend.position = "top", legend.text = element_text(size = 5.6),
    legend.key.width = grid::unit(5, "mm"),
    legend.key.height = grid::unit(2.3, "mm"),
    panel.grid = element_blank(), plot.margin = margin(4, 4, 3, 5)
  ) +
  guides(fill = "none", colour = guide_legend(override.aes = list(linewidth = 1.0)))

p_gene <- plot_cell_trajectory(
  cds, color_by = "target_expression", cell_size = 0.68,
  show_branch_points = FALSE, show_backbone = TRUE
) +
  scale_colour_gradientn(
    colours = c("#E3E5E4", "#F1D28A", "#D8834E", "#9D3F35"),
    limits = range(normalized_target), oob = squish,
    name = "log1p normalized\nexpression"
  ) +
  labs(title = target_gene) + trajectory_theme +
  theme(legend.position = "right")

p_gene_facet <- plot_cell_trajectory(
  cds, color_by = "target_expression", cell_size = 0.68,
  show_branch_points = FALSE, show_backbone = TRUE
) +
  facet_wrap(~Samples, nrow = 1) +
  scale_colour_gradientn(
    colours = c("#E3E5E4", "#F1D28A", "#D8834E", "#9D3F35"),
    limits = range(normalized_target), oob = squish,
    name = "log1p normalized\nexpression"
  ) +
  labs(title = paste0(target_gene, " by sample")) + trajectory_theme +
  theme(legend.position = "right", strip.text = element_text(size = 6.2, face = "bold"))

panel_a <- wrap_elements(full = p_state | p_sample | p_time)
preview <- (
  panel_a / (p_density | p_gene | p_gene_facet) +
    plot_layout(heights = c(1.0, 0.92), widths = c(0.36, 0.25, 0.39)) +
    plot_annotation(tag_levels = "A")
) &
  theme(
    plot.tag = element_text(
      family = tag_family, face = "plain", size = 12,
      colour = "#111111"
    ),
    plot.tag.position = c(0.006, 0.994)
  )

ragg::agg_png(
  file.path(outdir, "Figure6_pseudotime_preview.png"),
  width = 183 / 25.4, height = 125 / 25.4, units = "in", res = 300
)
print(preview)
dev.off()

write.csv(meta, file.path(outdir, "cell_pseudotime_source.csv"), row.names = FALSE)
write.csv(state_composition_final, file.path(outdir, "state_sample_composition.csv"), row.names = FALSE)
write.csv(qa, file.path(outdir, "pseudotime_QA_summary.csv"), row.names = FALSE)
write.csv(
  data.frame(gene_id = ordering_genes),
  file.path(outdir, "ordering_genes_HVG1000.csv"), row.names = FALSE
)
saveRDS(cds, file.path(outdir, "Epidermal_cells_I_monocle2_selected.rds"), compress = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
print(qa)
