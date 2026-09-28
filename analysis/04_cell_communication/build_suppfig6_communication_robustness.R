suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
})

set.seed(54)
args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(args) >= 1) args[[1]] else Sys.getenv("COMMUNICATION_ROBUSTNESS_DIR", unset = ".")
data_dir <- file.path(work_dir, "data")
out_dir <- file.path(work_dir, "output_suppfig6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

col_dcc <- "#B9654A"
col_ralf <- "#46768D"
col_dark <- "#24323D"
col_neutral <- "#D7DBDC"

theme_pc <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.30, colour = col_dark),
      axis.ticks = element_line(linewidth = 0.28, colour = col_dark),
      axis.text = element_text(size = 5.2, colour = "#303A42"),
      axis.title = element_text(size = 5.9, colour = col_dark),
      legend.title = element_text(size = 5.5),
      legend.text = element_text(size = 5.1),
      strip.background = element_rect(fill = "#F0F0ED", colour = NA),
      strip.text = element_text(size = 5.8, face = "bold", colour = col_dark),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

mapping <- read.delim(file.path(data_dir, "celltype_annotation_mapping.tsv"), check.names = FALSE)
map_vec <- setNames(mapping$final_celltype, mapping$old_celltype)
celltype_levels <- mapping$final_celltype

map_celltype <- function(x) {
  out <- unname(map_vec[x])
  ifelse(is.na(out), x, out)
}

core <- read.delim(
  gzfile(file.path(data_dir, "Figure6_core_route_audit.tsv.gz")),
  check.names = FALSE, stringsAsFactors = FALSE
) %>%
  filter(core_route %in% c(TRUE, "TRUE", "True", 1)) %>%
  mutate(
    module = recode(pathway_name, `DCC-DCCR` = "DCC1-DCCR1", `RALF-FER` = "RALF-FER"),
    source_final = map_celltype(source),
    target_final = map_celltype(target),
    route_index = row_number()
  )
stopifnot(sum(core$module == "DCC1-DCCR1") == 48L,
          sum(core$module == "RALF-FER") == 62L)

# A: distribution of directional stability for all 110 core routes.
stability_stats <- core %>% group_by(module) %>% summarise(
  n_routes = n(), minimum = min(direction_stability), median = median(direction_stability),
  .groups = "drop"
)
p_a <- ggplot(core, aes(module, direction_stability, fill = module)) +
  geom_violin(width = 0.72, trim = TRUE, alpha = 0.28, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.38, fill = "white") +
  geom_jitter(width = 0.11, height = 0, size = 0.62, alpha = 0.40, colour = col_dark) +
  geom_hline(yintercept = 0.90, linetype = "dashed", linewidth = 0.32, colour = "#8A9297") +
  geom_text(
    data = stability_stats,
    aes(module, 0.876, label = paste0("min ", sprintf("%.3f", minimum), "\nn = ", n_routes)),
    inherit.aes = FALSE, family = "Arial", size = 1.95, colour = col_dark
  ) +
  scale_fill_manual(values = c("DCC1-DCCR1" = col_dcc, "RALF-FER" = col_ralf)) +
  scale_y_continuous(limits = c(0.87, 1.005), breaks = seq(0.90, 1.00, 0.025)) +
  labs(x = NULL, y = "Direction stability", tag = "A") +
  theme_pc() + theme(legend.position = "none")

# B: route-wise effect direction and cell-sampling intervals.
route_df <- core %>%
  group_by(module) %>%
  arrange(normalized_delta_median, .by_group = TRUE) %>%
  mutate(route_rank = row_number()) %>% ungroup()
p_b <- ggplot(route_df, aes(route_rank, normalized_delta_median, colour = module)) +
  geom_hline(yintercept = 0, linewidth = 0.34, colour = "#7D858A") +
  geom_linerange(
    aes(ymin = normalized_delta_q025, ymax = normalized_delta_q975),
    linewidth = 0.32, alpha = 0.52
  ) +
  geom_point(size = 0.82) +
  facet_wrap(~module, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c("DCC1-DCCR1" = col_dcc, "RALF-FER" = col_ralf)) +
  labs(x = "Core-route rank", y = "Normalized FCR - CK change", tag = "B") +
  theme_pc() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

# C: full DCC1-A/B/D and DCCR1 expression basis across all final cell types.
dcc <- read.delim(
  gzfile(file.path(data_dir, "Figure6B_expression_source.tsv.gz")),
  check.names = FALSE, stringsAsFactors = FALSE
) %>%
  mutate(
    condition = toupper(condition),
    celltype = factor(map_celltype(celltype), levels = rev(celltype_levels)),
    gene = factor(gene, levels = c("DCC1-A", "DCC1-B", "DCC1-D", "DCCR1"))
  )

# D: RALF-FER candidate expression basis for the same cell-type order.
expr <- read.delim(
  gzfile(file.path(data_dir, "gene_celltype_condition_expression_summary.tsv.gz")),
  check.names = FALSE, stringsAsFactors = FALSE
)
ralf_genes <- c(
  "TraesCS5A02G093000", "TraesCS5B02G099100",
  "TraesCS5D02G105300", "TraesCS4A02G133800"
)
ralf_labels <- c(
  "TraesCS5A02G093000" = "RALF candidate-A",
  "TraesCS5B02G099100" = "RALF candidate-B",
  "TraesCS5D02G105300" = "RALF candidate-D",
  "TraesCS4A02G133800" = "FER candidate"
)
ralf <- expr %>% filter(gene_id %in% ralf_genes) %>%
  mutate(
    condition = toupper(condition),
    celltype = factor(map_celltype(celltype), levels = rev(celltype_levels)),
    gene = factor(unname(ralf_labels[gene_id]), levels = unname(ralf_labels[ralf_genes]))
  )
stopifnot(nrow(dcc) == 120L, nrow(ralf) == 120L)

make_dot <- function(df, tag, show_y = TRUE, mean_limit, fraction_limit,
                     mean_breaks, fraction_breaks) {
  plot_df <- df %>% mutate(fraction_display = pmin(fraction_expressing, fraction_limit))
  ggplot(plot_df, aes(gene, celltype)) +
    geom_point(aes(size = fraction_display, colour = mean_LogNormalize)) +
    facet_grid(. ~ condition) +
    scale_size_continuous(
      range = c(0.15, 3.30), limits = c(0, fraction_limit),
      breaks = fraction_breaks, labels = as.character(round(100 * fraction_breaks)),
      name = "Cells expressing (%)"
    ) +
    scale_colour_gradientn(
      colours = c("#E9E8E2", "#D3B37B", "#B76A4F", "#713A43"),
      limits = c(0, mean_limit), breaks = mean_breaks, oob = scales::squish,
      name = "Mean expression"
    ) +
    labs(x = NULL, y = NULL, tag = tag) +
    theme_pc() +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 34, hjust = 1, vjust = 1, size = 5.0),
      axis.text.y = if (show_y) element_text(size = 5.0) else element_blank(),
      legend.position = "bottom", legend.box = "vertical",
      legend.key.width = grid::unit(5.2, "mm"),
      panel.grid.major = element_line(colour = "#ECEDEA", linewidth = 0.25)
    )
}

p_c <- make_dot(
  dcc, "C", TRUE,
  mean_limit = 0.50, fraction_limit = 0.60,
  mean_breaks = seq(0, 0.5, 0.1), fraction_breaks = c(0.15, 0.30, 0.45, 0.60)
)
p_d <- make_dot(
  ralf, "D", FALSE,
  mean_limit = 1.60, fraction_limit = 0.90,
  mean_breaks = seq(0, 1.6, 0.4), fraction_breaks = c(0.25, 0.50, 0.75, 0.90)
)

tag_style <- element_text(family = "Arial", face = "bold", size = 8, colour = "#111111")
p_a <- p_a + theme(plot.tag = tag_style, plot.tag.position = c(0.006, 0.994))
p_b <- p_b + theme(plot.tag = tag_style, plot.tag.position = c(0.006, 0.994))
p_c <- p_c + theme(plot.tag = tag_style, plot.tag.position = c(0.006, 0.994))
p_d <- p_d + theme(
  plot.tag = tag_style, plot.tag.position = c(-0.045, 1.01),
  plot.margin = margin(4, 5, 4, 9)
)

figure <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(heights = c(0.72, 1.55), widths = c(1.05, 1))

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
  figure, "Supplementary_Figure6_communication_stability_expression",
  183, 200
)
save_three_formats(p_a, "Supplementary_Figure6A_direction_stability", 88, 72)
save_three_formats(p_b, "Supplementary_Figure6B_core_route_changes", 92, 72)
save_three_formats(p_c, "Supplementary_Figure6C_DCC_expression_basis", 100, 135)
save_three_formats(p_d, "Supplementary_Figure6D_RALF_FER_expression_basis", 100, 135)

preview <- file.path(out_dir, "Supplementary_Figure6_communication_stability_expression_preview.png")
agg_png(preview, width = 2161, height = 2362, res = 300, background = "white")
print(figure)
dev.off()

write.csv(stability_stats, file.path(out_dir, "panelA_direction_stability_summary.csv"), row.names = FALSE)
write.csv(route_df, file.path(out_dir, "panelB_core_route_sampling_intervals.csv"), row.names = FALSE)
write.csv(dcc, file.path(out_dir, "panelC_DCC_expression_basis.csv"), row.names = FALSE)
write.csv(ralf, file.path(out_dir, "panelD_RALF_FER_expression_basis.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

cat("PREVIEW=", preview, "\n", sep = "")
cat("DCC_ROUTES=", sum(core$module == "DCC1-DCCR1"), "\n", sep = "")
cat("RALF_ROUTES=", sum(core$module == "RALF-FER"), "\n", sep = "")
cat("DCC_MIN_STABILITY=", min(core$direction_stability[core$module == "DCC1-DCCR1"]), "\n", sep = "")
cat("RALF_MIN_STABILITY=", min(core$direction_stability[core$module == "RALF-FER"]), "\n", sep = "")
