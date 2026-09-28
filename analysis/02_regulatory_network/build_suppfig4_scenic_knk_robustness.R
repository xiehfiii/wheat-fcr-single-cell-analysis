suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
})

set.seed(54)
args <- commandArgs(trailingOnly = TRUE)
scenic_dir <- if (length(args) >= 1) args[[1]] else Sys.getenv("SCENIC_PROJECT_DIR", unset = ".")
out_dir <- if (length(args) >= 2) args[[2]] else file.path(".", "output_suppfig4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tf_id <- "TraesCS4A02G001300"
target_id <- "TraesCS3D02G094200"
col_main <- "#3E6F86"
col_warm <- "#BD684B"
col_target <- "#B24C4C"
col_neutral <- "#C8CFD3"
col_dark <- "#24323D"

theme_pc <- function(base_size = 6.6) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.30, colour = col_dark),
      axis.ticks = element_line(linewidth = 0.28, colour = col_dark),
      axis.text = element_text(size = 5.3, colour = "#303A42"),
      axis.title = element_text(size = 6.0, colour = col_dark),
      legend.title = element_text(size = 5.6),
      legend.text = element_text(size = 5.2),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

ranking <- read.csv(
  file.path(scenic_dir, "Target_Ranking_TraesCS4A02G001300.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
) %>% arrange(desc(importance)) %>% mutate(grn_rank = row_number())

ko <- read.csv(
  file.path(scenic_dir, "Gene_knockout_result_diffRegulation.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
) %>%
  filter(gene != tf_id, is.finite(Z), is.finite(p.adj)) %>%
  mutate(neglog10_fdr = -log10(pmax(p.adj, 1e-300)))

scenic_summary <- read.csv(
  file.path(scenic_dir, "SCENIC_Summary.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
tf_row <- scenic_summary %>% filter(TF_ID == tf_id) %>% slice(1)
stopifnot(nrow(tf_row) == 1L)
motif_targets <- unique(strsplit(tf_row$All_Targets, ";", fixed = TRUE)[[1]])
motif_targets <- motif_targets[nzchar(motif_targets)]
motif_downstream <- setdiff(motif_targets, tf_id)
stopifnot(length(motif_targets) == 44L, target_id %in% motif_targets)

ranking <- ranking %>%
  mutate(
    motif_supported = target %in% motif_downstream,
    target_gene = target == target_id
  )

standard_selected <- ko %>% filter(Z > 2, p.adj < 0.05)
standard_shared <- intersect(motif_downstream, standard_selected$gene)
stopifnot(length(standard_selected$gene) == 136L, length(standard_shared) == 4L,
          target_id %in% standard_shared)

# A: motif-pruned SCENIC targets in the complete GRNBoost2 ranking.
target_point <- ranking %>% filter(target_gene)
p_a <- ggplot(ranking, aes(grn_rank, importance)) +
  geom_line(linewidth = 0.42, colour = "#AEB8BE") +
  geom_point(
    data = ranking %>% filter(motif_supported),
    size = 0.78, colour = col_main, alpha = 0.86
  ) +
  geom_point(
    data = target_point, shape = 21, size = 2.25, stroke = 0.65,
    fill = "white", colour = col_target
  ) +
  annotate(
    "text", x = target_point$grn_rank + 42, y = target_point$importance + 0.16,
    label = "italic('TaMGBP1')~'(rank 45)'", parse = TRUE,
    family = "Arial", size = 2.05, colour = col_target, hjust = 0
  ) +
  annotate(
    "text", x = 1970, y = max(ranking$importance) * 0.96,
    label = "43 motif-pruned downstream targets", hjust = 1,
    family = "Arial", size = 2.05, colour = col_main
  ) +
  labs(x = "GRNBoost2 target rank", y = "Edge importance", tag = "A") +
  theme_pc()

# B: complete virtual-knockout response and the standard downstream call.
ko_plot <- ko %>%
  mutate(
    status = case_when(
      gene == target_id ~ "TaMGBP1",
      gene %in% standard_shared ~ "Other shared targets",
      Z > 2 & p.adj < 0.05 ~ "Perturbed genes",
      TRUE ~ "Other genes"
    ),
    status = factor(status, levels = c("Other genes", "Perturbed genes", "Other shared targets", "TaMGBP1"))
  )
target_ko <- ko_plot %>% filter(gene == target_id)
ko_ymax <- max(target_ko$neglog10_fdr * 1.22,
               quantile(ko_plot$neglog10_fdr, 0.995, na.rm = TRUE) * 1.05)
p_b <- ggplot(ko_plot, aes(Z, neglog10_fdr)) +
  geom_vline(xintercept = 2, linewidth = 0.32, linetype = "dashed", colour = "#75828A") +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.32, linetype = "dashed", colour = "#75828A") +
  geom_point(aes(colour = status), size = 0.48, alpha = 0.58) +
  geom_point(
    data = target_ko, shape = 21, size = 2.2, stroke = 0.65,
    fill = "white", colour = col_target
  ) +
  annotate(
    "text", x = target_ko$Z + 0.14, y = target_ko$neglog10_fdr + 0.30,
    label = "italic('TaMGBP1')", parse = TRUE, hjust = 0,
    family = "Arial", size = 2.05, colour = col_target
  ) +
  scale_colour_manual(values = c(
    "Other genes" = col_neutral,
    "Perturbed genes" = "#657985",
    "Other shared targets" = col_warm,
    "TaMGBP1" = col_target
  )) +
  coord_cartesian(ylim = c(0, ko_ymax)) +
  labs(x = "Virtual-knockout Z-score", y = expression(-log[10] * " adjusted p value"),
       colour = NULL, tag = "B") +
  theme_pc() +
  theme(legend.position = "top", legend.justification = "left",
        legend.key.width = grid::unit(3.8, "mm"))

# C: robustness of the SCENIC-knockout intersection to the calling thresholds.
z_cutoffs <- c(1.5, 1.75, 2.0, 2.25, 2.5, 3.0)
fdr_cutoffs <- c(1e-4, 1e-3, 1e-2, 0.05, 0.10)
grid_df <- expand.grid(z_cutoff = z_cutoffs, fdr_cutoff = fdr_cutoffs) %>%
  rowwise() %>%
  mutate(
    selected_n = sum(ko$Z > z_cutoff & ko$p.adj < fdr_cutoff),
    overlap_n = length(intersect(motif_downstream, ko$gene[ko$Z > z_cutoff & ko$p.adj < fdr_cutoff])),
    target_retained = target_id %in% ko$gene[ko$Z > z_cutoff & ko$p.adj < fdr_cutoff]
  ) %>% ungroup() %>%
  mutate(
    fdr_label = factor(
      ifelse(fdr_cutoff < 0.01, format(fdr_cutoff, scientific = TRUE), sprintf("%.2f", fdr_cutoff)),
      levels = rev(ifelse(fdr_cutoffs < 0.01, format(fdr_cutoffs, scientific = TRUE), sprintf("%.2f", fdr_cutoffs)))
    ),
    z_label = factor(sprintf("%.2f", z_cutoff), levels = sprintf("%.2f", z_cutoffs))
  )

p_c <- ggplot(grid_df, aes(z_label, fdr_label, fill = overlap_n)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_text(aes(label = overlap_n), family = "Arial", size = 2.25, colour = col_dark) +
  geom_point(
    data = grid_df %>% filter(target_retained),
    aes(shape = "TaMGBP1 retained"),
    size = 1.30, stroke = 0.42, fill = "white", colour = col_target,
    position = position_nudge(x = 0.28, y = 0.28)
  ) +
  geom_tile(
    data = grid_df %>% filter(z_cutoff == 2, abs(fdr_cutoff - 0.05) < 1e-10),
    fill = NA, colour = col_dark, linewidth = 0.8
  ) +
  scale_fill_gradient(low = "#F0EEE8", high = col_main, name = "Shared\ntargets") +
  scale_shape_manual(values = c("TaMGBP1 retained" = 23), name = NULL) +
  guides(shape = guide_legend(override.aes = list(fill = "white", colour = col_target))) +
  labs(x = "Z-score cutoff", y = "FDR cutoff", tag = "C") +
  theme_pc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "right")

# D: observed overlap versus random expectation across Z-score thresholds.
universe <- unique(ko$gene)
M <- length(intersect(motif_downstream, universe))
N <- length(universe)
sens_df <- bind_rows(lapply(z_cutoffs, function(zc) {
  selected <- ko %>% filter(Z > zc, p.adj < 0.05) %>% pull(gene) %>% unique()
  obs <- length(intersect(motif_downstream, selected))
  expected <- length(selected) * M / N
  pval <- phyper(obs - 1, M, N - M, length(selected), lower.tail = FALSE)
  data.frame(z_cutoff = zc, selected_n = length(selected), observed = obs,
             expected = expected, enrichment = ifelse(expected > 0, obs / expected, NA_real_),
             p_value = pval, target_retained = target_id %in% selected)
}))
sens_long <- sens_df %>% select(z_cutoff, observed, expected) %>%
  pivot_longer(c(observed, expected), names_to = "metric", values_to = "overlap") %>%
  mutate(metric = recode(metric, observed = "Observed overlap", expected = "Random expectation"))
std <- sens_df %>% filter(z_cutoff == 2)
p_d <- ggplot(sens_long, aes(z_cutoff, overlap, colour = metric, group = metric)) +
  geom_line(linewidth = 0.72) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = 2, linewidth = 0.30, linetype = "dashed", colour = "#7A858C") +
  scale_colour_manual(values = c("Observed overlap" = col_warm, "Random expectation" = "#7C8C95")) +
  scale_x_continuous(breaks = z_cutoffs) +
  labs(x = "Z-score cutoff (FDR < 0.05)", y = "SCENIC targets retained", colour = NULL, tag = "D") +
  theme_pc() +
  theme(legend.position = "top", legend.justification = "left")

figure <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(heights = c(1, 0.92), widths = c(1, 1.06)) &
  theme(
    plot.tag = element_text(family = "Arial", face = "bold", size = 8, colour = "#111111"),
    plot.tag.position = c(0.006, 0.994)
  )

save_three_formats <- function(plot_obj, stem, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  grDevices::cairo_pdf(
    file.path(out_dir, paste0(stem, ".pdf")),
    width = width_in, height = height_in, family = "Arial"
  )
  print(plot_obj)
  grDevices::dev.off()

  ragg::agg_png(
    file.path(out_dir, paste0(stem, ".png")),
    width = width_in, height = height_in, units = "in", res = dpi,
    background = "white"
  )
  print(plot_obj)
  grDevices::dev.off()

  ragg::agg_tiff(
    file.path(out_dir, paste0(stem, ".tiff")),
    width = width_in, height = height_in, units = "in", res = dpi,
    compression = "lzw", background = "white"
  )
  print(plot_obj)
  grDevices::dev.off()
}

# Submission-size composite and every constituent panel.
save_three_formats(
  figure, "Supplementary_Figure4_SCENIC_virtual_knockout_robustness",
  183, 150
)
save_three_formats(p_a, "Supplementary_Figure4A_SCENIC_target_ranking", 88, 72)
save_three_formats(p_b, "Supplementary_Figure4B_virtual_knockout_response", 92, 72)
save_three_formats(p_c, "Supplementary_Figure4C_threshold_sensitivity", 88, 72)
save_three_formats(p_d, "Supplementary_Figure4D_overlap_enrichment_sensitivity", 92, 72)

preview <- file.path(out_dir, "Supplementary_Figure4_SCENIC_virtual_knockout_robustness_preview.png")
agg_png(preview, width = 2161, height = 1772, res = 300, background = "white")
print(figure)
dev.off()

write.csv(ranking, file.path(out_dir, "panelA_SCENIC_target_ranking.csv"), row.names = FALSE)
write.csv(ko_plot, file.path(out_dir, "panelB_virtual_knockout_source.csv"), row.names = FALSE)
write.csv(grid_df, file.path(out_dir, "panelC_threshold_sensitivity.csv"), row.names = FALSE)
write.csv(sens_df, file.path(out_dir, "panelD_overlap_enrichment_sensitivity.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

cat("PREVIEW=", preview, "\n", sep = "")
cat("SCENIC_TARGETS=", length(motif_downstream), "\n", sep = "")
cat("STANDARD_KO_TARGETS=", nrow(standard_selected), "\n", sep = "")
cat("STANDARD_SHARED=", length(standard_shared), "\n", sep = "")
cat("STANDARD_ENRICHMENT=", round(std$enrichment, 3), "\n", sep = "")
cat("STANDARD_HYPERGEOM_P=", signif(std$p_value, 4), "\n", sep = "")
