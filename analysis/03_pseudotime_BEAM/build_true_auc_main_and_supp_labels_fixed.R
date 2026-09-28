suppressPackageStartupMessages({
  library(monocle)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ragg)
})

set.seed(54)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else Sys.getenv("PSEUDOTIME_PROJECT_DIR", unset = ".")

# Rebuild the established main-figure panels in the supplied temporary output
# directory, then replace panel C using the exact per-cell AUCell values.
base_script <- if (length(args) >= 2) args[[2]] else file.path(project_dir, "build_figure6_pseudotime_AE_final.R")
source(base_script, local = FALSE)

auc_source <- file.path(project_dir, "true_auc_preview", "true_TaERF87_AUCell_TaMGBP1_pseudotime.csv")
state_source <- file.path(project_dir, "selected_preview", "state_sample_composition.csv")
sens_source <- file.path(project_dir, "sensitivity", "trajectory_sensitivity_metrics.csv")

auc_df <- read.csv(auc_source, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(auc_df) == ncol(cds), !anyDuplicated(auc_df$cell_id))
ord <- match(colnames(cds), auc_df$cell_id)
stopifnot(!anyNA(ord))
auc_df <- auc_df[ord, , drop = FALSE]

theme_pc <- function(base_size = 6.8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.30, colour = "#292929"),
      axis.ticks = element_line(linewidth = 0.28, colour = "#292929"),
      axis.text = element_text(size = 5.3, colour = "#292929"),
      axis.title = element_text(size = 6.0, colour = "#292929"),
      plot.title = element_text(size = 7.0, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 5.6, colour = "#5B5B5B", hjust = 0),
      strip.background = element_rect(fill = "#F1F1EF", colour = NA),
      strip.text = element_text(size = 6.0, face = "plain", colour = "#292929"),
      legend.title = element_text(size = 5.6),
      legend.text = element_text(size = 5.2),
      panel.grid = element_blank(),
      plot.margin = margin(3, 5, 3, 5)
    )
}

# New main-panel C: true AUCell and target-gene expression versus pseudotime.
trend_df <- bind_rows(
  data.frame(
    cell_id = auc_df$cell_id,
    Pseudotime = auc_df$Pseudotime,
    value = auc_df$TaERF87_AUCell,
    feature = "TaERF87 regulon activity (AUCell)"
  ),
  data.frame(
    cell_id = auc_df$cell_id,
    Pseudotime = auc_df$Pseudotime,
    value = auc_df$TaMGBP1_expression,
    feature = "TaMGBP1 expression"
  )
) %>%
  mutate(feature = factor(
    feature,
    levels = c("TaERF87 regulon activity (AUCell)", "TaMGBP1 expression")
  ))

strip_labels <- c(
  "TaERF87 regulon activity (AUCell)" = "italic('TaERF87')~regulon~activity~'(AUCell)'",
  "TaMGBP1 expression" = "italic('TaMGBP1')~expression"
)

p_trend <- ggplot(trend_df, aes(Pseudotime, value)) +
  geom_point(size = 0.24, alpha = 0.16, colour = "#858C91") +
  geom_smooth(
    method = "gam", formula = y ~ s(x, k = 6), se = FALSE,
    linewidth = 0.90, colour = "#B75D3E"
  ) +
  facet_wrap(
    ~feature, scales = "free_y", nrow = 1,
    labeller = as_labeller(strip_labels, label_parsed)
  ) +
  labs(x = "Pseudotime", y = NULL, tag = "C") +
  theme_pc() +
  theme(
    plot.tag = element_text(family = "Arial", face = "bold", size = 8),
    plot.tag.position = c(0.006, 0.994)
  )

# Preserve the biological geometry of the DDRTree embedding: one unit on
# Component 1 now occupies the same physical length as one unit on Component 2.
# The larger top-row allocation prevents these fixed-aspect panels from looking
# vertically compressed in the assembled figure.
# A small optical correction is retained after patchwork assembly because the
# long horizontal labels and legends otherwise make an exactly 1:1 coordinate
# panel appear vertically compressed at final page size.
p_state <- p_state + coord_fixed(ratio = 1.14)
p_sample <- p_sample + coord_fixed(ratio = 1.14)
p_time <- p_time + coord_fixed(ratio = 1.14)
top_new <- p_state | p_sample | p_time
middle_new <- (p_density | p_trend) + plot_layout(widths = c(0.30, 0.70))
main_new <- top_new / middle_new / bp2_panel$panel / bp1_panel$panel +
  plot_layout(heights = c(1.02, 0.68, 0.94, 0.94)) &
  theme(
    plot.tag = element_text(family = "Arial", face = "bold", size = 8, colour = "#111111"),
    plot.tag.position = c(0.006, 0.994)
  )

# Supplement panels A-B: composition and sensitivity.
pal_sample <- c(CK = "#2F6C8E", FCR = "#D98C3F")
state_df <- read.csv(state_source, stringsAsFactors = FALSE) %>%
  mutate(
    Samples = toupper(Samples),
    State = factor(State, levels = sort(unique(as.integer(State)))),
    Samples = factor(Samples, levels = c("CK", "FCR")),
    label = paste0(n, "\n", sprintf("%.1f%%", 100 * fraction))
  ) %>%
  group_by(State) %>%
  arrange(factor(Samples, levels = c("FCR", "CK")), .by_group = TRUE) %>%
  mutate(
    y_mid = cumsum(fraction) - fraction / 2,
    outside_top = Samples == "CK" & fraction < 0.10,
    y_label = case_when(
      # Keep the two-line labels comfortably below the panel ceiling after
      # patchwork rescales panel A; this prevents the State 2 label from being
      # visually clipped in the assembled figure.
      outside_top ~ 1.095,
      TRUE ~ y_mid
    ),
    label_colour = ifelse(outside_top, "#2F6C8E", "white")
  ) %>%
  ungroup()

p_state_comp <- ggplot(state_df, aes(State, fraction, fill = Samples)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.30) +
  geom_segment(
    data = state_df %>% filter(outside_top),
    aes(x = as.numeric(State), xend = as.numeric(State), y = 1.005, yend = 1.025),
    inherit.aes = FALSE, colour = "#738A9A", linewidth = 0.28
  ) +
  geom_text(
    aes(y = y_label, label = label, colour = label_colour),
    size = 1.76, family = "Arial", lineheight = 0.72
  ) +
  scale_fill_manual(values = pal_sample, breaks = c("CK", "FCR")) +
  scale_colour_identity() +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    breaks = seq(0, 1, 0.25), limits = c(0, 1.19),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Transcriptional State", y = "Cell fraction", fill = NULL, tag = "A") +
  theme_pc(7.0) +
  theme(legend.position = "top", legend.justification = "left")

sens <- read.csv(sens_source, stringsAsFactors = FALSE) %>%
  filter(run %in% c("unsupervised_HVG1000", "unsupervised_HVG2000", "unsupervised_HVG2000_dim10")) %>%
  mutate(
    setting = recode(
      run,
      unsupervised_HVG1000 = "HVG 1,000\n15 dimensions",
      unsupervised_HVG2000 = "HVG 2,000\n15 dimensions",
      unsupervised_HVG2000_dim10 = "HVG 2,000\n10 dimensions"
    ),
    setting = factor(setting, levels = c(
      "HVG 1,000\n15 dimensions", "HVG 2,000\n15 dimensions", "HVG 2,000\n10 dimensions"
    ))
  )
sens_long <- sens %>%
  select(setting, states, branch_points) %>%
  pivot_longer(c(states, branch_points), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, states = "States", branch_points = "Branch points"))

p_sens <- ggplot(sens_long, aes(setting, value, colour = metric, group = metric)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.1) +
  geom_text(aes(label = value), vjust = -0.75, size = 2.15, family = "Arial", show.legend = FALSE) +
  geom_text(
    data = sens,
    aes(setting, 0.65, label = paste0("AUC ", sprintf("%.3f", fcr_vs_ck_pseudotime_auc))),
    inherit.aes = FALSE, size = 1.95, family = "Arial", colour = "#606060"
  ) +
  scale_colour_manual(values = c("States" = "#596F8A", "Branch points" = "#C26B45")) +
  scale_y_continuous(breaks = 0:7, limits = c(0.2, 7.8), expand = c(0, 0)) +
  labs(x = NULL, y = "Number", colour = NULL, tag = "B") +
  theme_pc(7.0) +
  theme(legend.position = "top", legend.justification = "left", axis.text.x = element_text(size = 5.8))

# Supplement panels C-D: expression projected onto the established trajectory.
gene_erf <- "TraesCS4A02G001300"
gene_mgbp <- "TraesCS3D02G094200"
meta2 <- as.data.frame(pData(cds))
meta2$Samples <- factor(toupper(as.character(meta2$Samples)), levels = c("CK", "FCR"))
pData(cds)$Samples <- meta2$Samples
lognorm_gene <- function(g) log1p(as.numeric(exprs(cds)[g, ]) / as.numeric(sizeFactors(cds)))
pData(cds)$TaERF87_expression <- lognorm_gene(gene_erf)
pData(cds)$TaMGBP1_expression <- lognorm_gene(gene_mgbp)

make_expression_maps <- function(column_name, symbol, panel_tag) {
  vals <- pData(cds)[[column_name]]
  sc <- scale_colour_gradientn(
    colours = c("#E5E7E6", "#E7C980", "#CA7548", "#873A37"),
    limits = range(vals, na.rm = TRUE), oob = scales::squish,
    name = "log1p normalized\nexpression"
  )
  p_all <- plot_cell_trajectory(
    cds, color_by = column_name, cell_size = 0.52,
    show_branch_points = FALSE, show_backbone = TRUE
  ) + sc +
    labs(title = bquote(italic(.(symbol))), tag = panel_tag) + theme_traj +
    theme(legend.position = "none")
  p_by <- plot_cell_trajectory(
    cds, color_by = column_name, cell_size = 0.52,
    show_branch_points = FALSE, show_backbone = TRUE
  ) + facet_wrap(~Samples, nrow = 1) + sc +
    labs(title = bquote(italic(.(symbol))~by~sample)) + theme_traj +
    theme(legend.position = "right", strip.text = element_text(size = 5.7, face = "bold"))
  (p_all | p_by) + plot_layout(widths = c(0.37, 0.63))
}

p_erf_maps <- make_expression_maps("TaERF87_expression", "TaERF87", "C")
p_mgbp_maps <- make_expression_maps("TaMGBP1_expression", "TaMGBP1", "D")

make_overall_expression_map <- function(column_name, symbol, panel_tag = NULL) {
  vals <- pData(cds)[[column_name]]
  sc <- scale_colour_gradientn(
    colours = c("#E5E7E6", "#E7C980", "#CA7548", "#873A37"),
    limits = range(vals, na.rm = TRUE), oob = scales::squish,
    name = "log1p normalized\nexpression"
  )
  plot_cell_trajectory(
    cds, color_by = column_name, cell_size = 0.54,
    show_branch_points = FALSE, show_backbone = TRUE
  ) + sc +
    labs(title = bquote(italic(.(symbol))), tag = panel_tag) + theme_traj +
    theme(legend.position = "right")
}

p_erf_overall <- make_overall_expression_map("TaERF87_expression", "TaERF87", "C")
p_mgbp_overall <- make_overall_expression_map("TaMGBP1_expression", "TaMGBP1")
p_expression_overall <- p_erf_overall | p_mgbp_overall

# Supplement panel E: fate-dependent expression trends from both branch points.
lognorm_object <- function(object, genes) {
  x <- as.matrix(exprs(object)[genes, , drop = FALSE])
  log1p(sweep(x, 2, as.numeric(sizeFactors(object)), "/"))
}
extract_branch <- function(bp) {
  bc <- buildBranchCellDataSet(
    cds, branch_point = bp, progenitor_method = "duplicate",
    branch_labels = c("Branch_A", "Branch_B")
  )
  bm <- as.data.frame(pData(bc))
  vals <- lognorm_object(bc, c(gene_erf, gene_mgbp))
  branch_names <- if (bp == 1L) {
    c(Branch_A = "State 4 fate", Branch_B = "State 5 fate")
  } else {
    c(Branch_A = "State 2 fate", Branch_B = "States 4/5 fate")
  }
  bind_rows(lapply(seq_len(nrow(vals)), function(i) {
    data.frame(
      branch_point = paste0("Branch point ", bp),
      gene = if (rownames(vals)[i] == gene_erf) "TaERF87" else "TaMGBP1",
      Pseudotime = bm$Pseudotime,
      expression = as.numeric(vals[i, ]),
      fate = unname(branch_names[as.character(bm$Branch)])
    )
  }))
}

branch_df <- bind_rows(extract_branch(1L), extract_branch(2L)) %>%
  mutate(
    branch_point = factor(branch_point, levels = c("Branch point 1", "Branch point 2")),
    gene = factor(gene, levels = c("TaERF87", "TaMGBP1"))
  )

p_branch <- ggplot(branch_df, aes(Pseudotime, expression, colour = fate)) +
  geom_point(size = 0.18, alpha = 0.07) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = FALSE, linewidth = 0.76) +
  facet_grid(
    branch_point ~ gene, scales = "free_y",
    labeller = labeller(gene = as_labeller(
      c(TaERF87 = "italic('TaERF87')", TaMGBP1 = "italic('TaMGBP1')"), label_parsed
    ))
  ) +
  scale_colour_manual(values = c(
    "State 2 fate" = "#C06A3D", "State 4 fate" = "#A94F4F",
    "State 5 fate" = "#467A9B", "States 4/5 fate" = "#548E85"
  )) +
  labs(x = "Branched pseudotime", y = "log1p normalized expression", colour = NULL, tag = "D") +
  theme_pc(6.8) +
  theme(legend.position = "top", legend.justification = "left", strip.text.x = element_text(face = "plain"))

supp_new <- ((p_state_comp | p_sens) + plot_layout(widths = c(1, 1.10))) /
  p_expression_overall / p_branch +
  plot_layout(heights = c(0.92, 1.22, 1.40)) &
  theme(
    plot.tag = element_text(family = "Arial", face = "bold", size = 8, colour = "#111111"),
    plot.tag.position = c(0.006, 0.994)
  )

main_preview <- file.path(outdir, "Figure6_pseudotime_true_AUCell_main_preview.png")
agg_png(main_preview, width = 2161, height = 2717, res = 300, background = "white")
print(main_new)
dev.off()

supp_preview <- file.path(outdir, "Supplementary_pseudotime_streamlined_labels_fixed_preview.png")
agg_png(supp_preview, width = 2161, height = 2540, res = 300, background = "white")
print(supp_new)
dev.off()

panel_a_preview <- file.path(outdir, "Supplementary_pseudotime_panelA_labels_fixed_preview.png")
agg_png(panel_a_preview, width = 1063, height = 827, res = 300, background = "white")
print(p_state_comp)
dev.off()

save_three_formats <- function(plot_obj, stem, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  grDevices::cairo_pdf(
    file.path(outdir, paste0(stem, ".pdf")),
    width = width_in, height = height_in, family = "Arial", onefile = TRUE
  )
  print(plot_obj)
  grDevices::dev.off()

  ragg::agg_png(
    file.path(outdir, paste0(stem, ".png")),
    width = width_in, height = height_in, units = "in", res = dpi,
    background = "white"
  )
  print(plot_obj)
  grDevices::dev.off()

  ragg::agg_tiff(
    file.path(outdir, paste0(stem, ".tiff")),
    width = width_in, height = height_in, units = "in", res = dpi,
    compression = "lzw", background = "white"
  )
  print(plot_obj)
  grDevices::dev.off()
}

# Formal export of the complete supplementary figure and every constituent panel.
save_three_formats(
  supp_new, "Supplementary_pseudotime_robustness_and_gene_dynamics",
  183, 215
)
save_three_formats(
  p_state_comp, "Supplementary_pseudotime_panelA_state_composition",
  88, 80
)
save_three_formats(
  p_sens, "Supplementary_pseudotime_panelB_parameter_sensitivity",
  95, 80
)
save_three_formats(
  p_expression_overall, "Supplementary_pseudotime_panelC_gene_expression_trajectories",
  183, 96
)
save_three_formats(
  p_branch, "Supplementary_pseudotime_panelD_branch_dynamics",
  183, 105
)

write.csv(trend_df, file.path(outdir, "main_panelC_true_AUCell_source.csv"), row.names = FALSE)
write.csv(branch_df, file.path(outdir, "supp_panelE_branch_dynamics_source.csv"), row.names = FALSE)

# Export only the modified main-figure constituent (A) and the reassembled
# main figure. Unchanged panels are deliberately not replaced downstream.
save_three_formats(
  top_new, "Figure6A_pseudotime_trajectory",
  183, 68
)

# Formal full-size export of the approved main figure.
if (TRUE) {
formal_stem <- file.path(outdir, "Figure6_pseudotime_true_AUCell")
formal_width <- 183 / 25.4
formal_height <- 230 / 25.4

agg_png(
  paste0(formal_stem, ".png"), width = formal_width, height = formal_height,
  units = "in", res = 600, background = "white"
)
print(main_new)
dev.off()

agg_tiff(
  paste0(formal_stem, ".tiff"), width = formal_width, height = formal_height,
  units = "in", res = 600, compression = "lzw", background = "white"
)
print(main_new)
dev.off()

grDevices::cairo_pdf(
  paste0(formal_stem, ".pdf"), width = formal_width, height = formal_height,
  family = "Arial", onefile = TRUE
)
print(main_new)
dev.off()
}

cat("MAIN_PREVIEW=", main_preview, "\n", sep = "")
cat("SUPP_PREVIEW=", supp_preview, "\n", sep = "")
cat("PANEL_A_PREVIEW=", panel_a_preview, "\n", sep = "")
