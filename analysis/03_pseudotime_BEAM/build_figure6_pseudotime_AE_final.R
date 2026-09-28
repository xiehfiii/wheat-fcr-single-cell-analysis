suppressPackageStartupMessages({
  library(monocle)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ragg)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Usage: Rscript build_figure6_pseudotime_AE_preview.R SELECTED_CDS_RDS BP2_HEATMAP_RDS BP2_GO_CSV BP1_HEATMAP_RDS BP1_GO_CSV OUTPUT_DIRECTORY")
}

set.seed(54)
cds_file <- normalizePath(args[[1]], mustWork = TRUE)
bp2_heatmap_file <- normalizePath(args[[2]], mustWork = TRUE)
bp2_go_file <- normalizePath(args[[3]], mustWork = TRUE)
bp1_heatmap_file <- normalizePath(args[[4]], mustWork = TRUE)
bp1_go_file <- normalizePath(args[[5]], mustWork = TRUE)
outdir <- args[[6]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

font_body <- "Arial"
font_tag <- "Arial"
col_ck <- "#2F6C8E"
col_fcr <- "#D98C3F"
col_state4 <- "#B8664F"
col_state5 <- "#4F7698"
col_pre <- "#C9C9C7"
state_colours <- c("#315F78", "#4C956C", "#D6A84B", "#B56576", "#66578A")

cds <- readRDS(cds_file)
bp2_comp <- readRDS(bp2_heatmap_file)
bp2_go <- read.csv(bp2_go_file, stringsAsFactors = FALSE, check.names = FALSE)
bp1_comp <- readRDS(bp1_heatmap_file)
bp1_go <- read.csv(bp1_go_file, stringsAsFactors = FALSE, check.names = FALSE)
target_gene <- "TraesCS4A02G001300"

meta <- as.data.frame(pData(cds))
meta$cell_id <- rownames(meta)
meta$Samples <- factor(tolower(as.character(meta$Samples)), levels = c("ck", "fcr"), labels = c("CK", "FCR"))
meta$target_expression <- log1p(as.numeric(exprs(cds)[target_gene, ]) / as.numeric(sizeFactors(cds)))
pData(cds)$Samples <- meta$Samples
pData(cds)$target_expression <- meta$target_expression
coords <- as.data.frame(t(reducedDimS(cds)))
colnames(coords)[1:2] <- c("Component_1", "Component_2")
coords$cell_id <- rownames(coords)
source_cells <- left_join(meta %>% select(cell_id, Samples, State, Pseudotime,
                                          percent.mito, percent.chloro, target_expression),
                          coords, by = "cell_id")

theme_traj <- theme_classic(base_size = 6.8, base_family = font_body) +
  theme(
    axis.line = element_line(linewidth = 0.30, colour = "#292929"),
    axis.ticks = element_line(linewidth = 0.28, colour = "#292929"),
    axis.text = element_text(size = 5.2, colour = "#292929"),
    axis.title = element_text(size = 6.0, colour = "#292929"),
    plot.title = element_text(size = 7.0, face = "bold", hjust = 0.5, margin = margin(b = 1.5)),
    legend.title = element_text(size = 5.6),
    legend.text = element_text(size = 5.2),
    legend.key.height = grid::unit(2.7, "mm"),
    legend.key.width = grid::unit(3.6, "mm"),
    legend.position = "top",
    panel.grid = element_blank(),
    plot.margin = margin(2, 3, 2, 3)
  )

p_state <- plot_cell_trajectory(
  cds, color_by = "State", cell_size = 0.56,
  show_branch_points = FALSE, show_backbone = TRUE
) +
  scale_colour_manual(values = state_colours) +
  labs(title = "Transcriptional state", tag = "A") + theme_traj

p_sample <- plot_cell_trajectory(
  cds, color_by = "Samples", cell_size = 0.56,
  show_branch_points = FALSE, show_backbone = TRUE
) +
  scale_colour_manual(values = c(CK = col_ck, FCR = col_fcr)) +
  labs(title = "Sample origin", colour = "Sample") + theme_traj

p_time <- plot_cell_trajectory(
  cds, color_by = "Pseudotime", cell_size = 0.56,
  show_branch_points = FALSE, show_backbone = TRUE
) +
  scale_colour_gradientn(
    colours = c("#F3F2E9", "#BFD7C8", "#6EA69A", "#235E63"),
    name = "Pseudotime"
  ) +
  labs(title = "Pseudotime") + theme_traj

medians <- meta %>% group_by(Samples) %>% summarise(median = median(Pseudotime), .groups = "drop")
p_density <- ggplot(meta, aes(Pseudotime, fill = Samples, colour = Samples)) +
  geom_density(alpha = 0.22, linewidth = 0.62, adjust = 1.05, key_glyph = "path") +
  geom_vline(data = medians, aes(xintercept = median, colour = Samples),
             linewidth = 0.50, linetype = "22", show.legend = FALSE) +
  scale_fill_manual(values = c(CK = col_ck, FCR = col_fcr)) +
  scale_colour_manual(values = c(CK = col_ck, FCR = col_fcr)) +
  labs(x = "Pseudotime", y = "Cell density", colour = NULL, tag = "B") +
  theme_classic(base_size = 6.8, base_family = font_body) +
  theme(
    axis.line = element_line(linewidth = 0.30), axis.ticks = element_line(linewidth = 0.28),
    axis.text = element_text(size = 5.3), axis.title = element_text(size = 6.0),
    legend.position = "top", legend.text = element_text(size = 5.3),
    legend.key.width = grid::unit(4.8, "mm"), panel.grid = element_blank(),
    plot.margin = margin(4, 5, 3, 4)
  ) + guides(fill = "none", colour = guide_legend(override.aes = list(fill = NA, alpha = 1, linewidth = 1.1, linetype = 1)))

expr_scale <- scale_colour_gradientn(
  colours = c("#E4E6E5", "#EBCB82", "#CF7B4B", "#913B36"),
  limits = range(meta$target_expression), oob = squish,
  name = "log1p normalized\nexpression"
)
p_gene <- plot_cell_trajectory(
  cds, color_by = "target_expression", cell_size = 0.58,
  show_branch_points = FALSE, show_backbone = TRUE
) + expr_scale + labs(title = target_gene, tag = "C") + theme_traj + theme(legend.position = "none")

p_gene_facet <- plot_cell_trajectory(
  cds, color_by = "target_expression", cell_size = 0.58,
  show_branch_points = FALSE, show_backbone = TRUE
) +
  facet_wrap(~Samples, nrow = 1) + expr_scale +
  labs(title = paste0(target_gene, " by sample")) + theme_traj +
  theme(legend.position = "right", strip.text = element_text(size = 5.7, face = "bold"))

ratio_to_num <- function(x) vapply(strsplit(x, "/", fixed = TRUE), function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1))
module_labels <- c(M1 = "M1", M2 = "M2")

prepare_go <- function(go) {
  go %>%
    mutate(
      gene_ratio = ratio_to_num(GeneRatio),
      neglog10_fdr = -log10(pmax(p.adjust, .Machine$double.xmin)),
      description_wrapped = vapply(
        Description,
        function(x) paste(strwrap(x, width = 31), collapse = "\n"),
        character(1)
      ),
      module_label = factor(module_labels[module], levels = module_labels)
    ) %>%
    arrange(module_label, desc(neglog10_fdr)) %>%
    mutate(term_order = factor(paste(module, ID), levels = rev(paste(module, ID))))
}

bp2_go <- prepare_go(bp2_go)
bp1_go <- prepare_go(bp1_go)
all_go <- bind_rows(bp2_go, bp1_go)
go_ratio_limits <- c(0, max(all_go$gene_ratio, na.rm = TRUE) * 1.06)
go_count_limits <- range(all_go$Count, na.rm = TRUE)
go_fdr_limits <- range(all_go$neglog10_fdr, na.rm = TRUE)

build_branch_panel <- function(comp, go, panel_tag, branch_title, segment_levels) {
  mat <- comp$matrix
  tree_order <- rownames(mat)[comp$tree$order]
  module_key <- data.frame(
    gene_id = names(comp$modules),
    module = paste0("M", as.integer(comp$modules)),
    stringsAsFactors = FALSE
  )
  heat <- data.frame(gene_id = rownames(mat), mat, check.names = FALSE) %>%
    left_join(module_key, by = "gene_id") %>%
    pivot_longer(-c(gene_id, module), names_to = "matrix_column", values_to = "z") %>%
    mutate(
      matrix_index = as.integer(matrix_column),
      gene_id = factor(gene_id, levels = rev(tree_order)),
      module_label = factor(module_labels[module], levels = module_labels)
    )

  col_anno <- data.frame(
    matrix_index = seq_len(ncol(mat)),
    segment = factor(comp$annotation_col[[1]], levels = segment_levels)
  )
  segment_centres <- col_anno %>%
    group_by(segment) %>%
    summarise(x = mean(matrix_index), .groups = "drop") %>%
    filter(!is.na(segment))
  segment_values <- setNames(c(col_state4, col_pre, col_state5), segment_levels)

  p_colanno <- ggplot(col_anno, aes(matrix_index, 1, fill = segment)) +
    geom_raster() +
    scale_fill_manual(values = segment_values, guide = "none", na.value = col_pre) +
    scale_x_continuous(
      position = "top", breaks = segment_centres$x,
      labels = as.character(segment_centres$segment), expand = c(0, 0)
    ) +
    coord_cartesian(expand = FALSE) +
    labs(tag = panel_tag, title = branch_title) +
    theme_void(base_family = font_body) +
    theme(
      plot.title = element_text(size = 6.8, face = "bold", hjust = 0.5, margin = margin(b = 2)),
      axis.text.x.top = element_text(size = 5.2, colour = "#222222", margin = margin(b = 1)),
      plot.margin = margin(2, 2, 0, 1)
    )

  p_heat <- ggplot(heat, aes(matrix_index, gene_id, fill = z)) +
    geom_raster(interpolate = FALSE) +
    facet_grid(module_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2(
      low = "#315E7A", mid = "#F4F2EC", high = "#B65C4B",
      midpoint = 0, limits = c(-3, 3), oob = squish, name = "Row z-score"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 6.4, base_family = font_body) +
    theme(
      panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
      strip.placement = "outside", strip.background = element_blank(),
      strip.text.y.left = element_text(size = 5.4, angle = 0, face = "bold"),
      panel.spacing.y = grid::unit(0.7, "mm"),
      legend.position = "right", legend.title = element_text(size = 5.2),
      legend.text = element_text(size = 5.0), legend.key.height = grid::unit(7.2, "mm"),
      plot.margin = margin(0, 2, 2, 1)
    )
  p_heat_full <- p_colanno / p_heat + plot_layout(heights = c(0.12, 0.88))
  p_heat_full_clean <- (p_colanno + labs(tag = NULL)) / p_heat +
    plot_layout(heights = c(0.12, 0.88))

  p_go <- ggplot(go, aes(gene_ratio, term_order, size = Count, colour = neglog10_fdr)) +
    geom_point(alpha = 0.92) +
    facet_grid(module_label ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = setNames(go$description_wrapped, paste(go$module, go$ID))) +
    scale_x_continuous(limits = go_ratio_limits, expand = expansion(mult = c(0, 0.03))) +
    scale_size_continuous(
      range = c(1.6, 3.8), limits = go_count_limits,
      breaks = pretty(go_count_limits, n = 3), name = "Gene count"
    ) +
    scale_colour_gradient(
      low = "#B9C7CE", high = "#304F68", limits = go_fdr_limits,
      oob = squish, name = "-log10 FDR"
    ) +
    guides(
      colour = guide_colourbar(
        order = 1, title.position = "top",
        barheight = grid::unit(10.5, "mm"), barwidth = grid::unit(2.5, "mm")
      ),
      size = guide_legend(
        order = 2, title.position = "top",
        keyheight = grid::unit(2.6, "mm"), keywidth = grid::unit(4.0, "mm")
      )
    ) +
    labs(x = "Gene ratio", y = NULL) +
    theme_classic(base_size = 6.4, base_family = font_body) +
    theme(
      axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.x = element_line(linewidth = 0.30), axis.ticks.x = element_line(linewidth = 0.28),
      axis.text.x = element_text(size = 5.1), axis.text.y = element_text(size = 5.0, lineheight = 0.94),
      axis.title.x = element_text(size = 5.8),
      strip.background = element_blank(), strip.text = element_text(size = 5.3, angle = 0, face = "bold"),
      panel.spacing.y = grid::unit(1.0, "mm"),
      legend.position = "right", legend.title = element_text(size = 5.1),
      legend.text = element_text(size = 5.0),
      legend.spacing.y = grid::unit(1.0, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(4, 3, 2, 4)
    )

  list(
    panel = (p_heat_full | p_go) + plot_layout(widths = c(0.51, 0.49)),
    panel_clean = (p_heat_full_clean | p_go) + plot_layout(widths = c(0.51, 0.49)),
    heat_source = heat,
    go_source = go
  )
}

bp2_panel <- build_branch_panel(
  bp2_comp, bp2_go, "D", "Branch point 2",
  c("State 2 terminal", "Pre-branch", "Other downstream states")
)
bp1_panel <- build_branch_panel(
  bp1_comp, bp1_go, "E", "Branch point 1",
  c("State 4 terminal", "Pre-branch", "State 5 terminal")
)

top <- p_state | p_sample | p_time
middle_gene <- (p_gene | p_gene_facet) + plot_layout(widths = c(0.38, 0.62))
middle <- (p_density | middle_gene) + plot_layout(widths = c(0.30, 0.70))
figure <- top / middle / bp2_panel$panel / bp1_panel$panel +
  plot_layout(heights = c(0.72, 0.68, 0.94, 0.94)) &
  theme(
    plot.tag = element_text(family = font_tag, face = "bold", size = 8, colour = "#111111"),
    plot.tag.position = c(0.006, 0.994)
  )

figure_width_in <- 183 / 25.4
figure_height_in <- 230 / 25.4

panel_a <- (p_state + labs(tag = NULL)) | p_sample | p_time
panel_b <- p_density + labs(tag = NULL)
panel_c <- ((p_gene + labs(tag = NULL)) | p_gene_facet) + plot_layout(widths = c(0.38, 0.62))
panel_d <- bp2_panel$panel_clean
panel_e <- bp1_panel$panel_clean

save_plot_bundle <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  ragg::agg_png(
    file.path(outdir, paste0(stem, ".png")),
    width = width_in, height = height_in, units = "in", res = 600,
    background = "white"
  )
  print(plot)
  dev.off()

  ragg::agg_tiff(
    file.path(outdir, paste0(stem, ".tiff")),
    width = width_in, height = height_in, units = "in", res = 600,
    compression = "lzw", background = "white"
  )
  print(plot)
  dev.off()

  grDevices::cairo_pdf(
    file.path(outdir, paste0(stem, ".pdf")),
    width = width_in, height = height_in, family = font_body,
    onefile = TRUE
  )
  print(plot)
  dev.off()
}

save_plot_bundle(figure, "Figure6_pseudotime_revised", 183, 230)
save_plot_bundle(panel_a, "Figure6A_pseudotime_trajectory", 183, 68)
save_plot_bundle(panel_b, "Figure6B_pseudotime_density", 58, 58)
save_plot_bundle(panel_c, "Figure6C_target_gene_trajectory", 128, 58)
save_plot_bundle(panel_d, "Figure6D_branch_point2_BEAM_GO", 183, 59)
save_plot_bundle(panel_e, "Figure6E_branch_point1_BEAM_GO", 183, 59)
