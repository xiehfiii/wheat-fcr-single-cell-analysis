suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(magick)
  library(grid)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: script panelA_png source_dir output_dir")
panel_a_path <- args[[1]]
source_dir <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theme_pc <- function(base_size = 7.2) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(
      axis.line = element_line(colour = "#222222", linewidth = 0.35),
      axis.ticks = element_line(colour = "#222222", linewidth = 0.35),
      axis.text = element_text(colour = "#222222", size = base_size - 0.5),
      axis.title = element_text(colour = "#222222", size = base_size),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      plot.title = element_text(
        size = base_size + 0.5, face = "bold", colour = "#111111",
        hjust = 0.5, margin = margin(b = 3)
      ),
      plot.subtitle = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )
}

save_panel_preview <- function(plot, filename, width_mm, height_mm) {
  agg_png(file.path(output_dir, filename), width = width_mm, height = height_mm,
          units = "mm", res = 300, background = "white")
  print(plot)
  dev.off()
}

save_panel_formal <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  grDevices::cairo_pdf(file.path(output_dir, paste0(stem, ".pdf")),
                       width = width_in, height = height_in,
                       family = "Arial", onefile = TRUE)
  print(plot); dev.off()
  ragg::agg_png(file.path(output_dir, paste0(stem, ".png")),
                width = width_in, height = height_in, units = "in",
                res = 600, background = "white")
  print(plot); dev.off()
  ragg::agg_tiff(file.path(output_dir, paste0(stem, ".tiff")),
                 width = width_in, height = height_in, units = "in",
                 res = 600, compression = "lzw", background = "white")
  print(plot); dev.off()
}

# B: QC association. All secondary grey annotations are deliberately removed.
qc_long <- read.delim(file.path(source_dir, "01_qc_correlation_source.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE)
# The original audit table stored the intended two-line display label with a
# literal newline. Reconstruct it from the intact condition/model columns so
# the four QC columns remain unambiguous.
qc_long <- qc_long[is.finite(qc_long$rho), , drop = FALSE]
qc_long$column <- paste(qc_long$condition,
                        ifelse(qc_long$model == "Raw Spearman", "Raw", "Adjusted"),
                        sep = "\n")
qc_long$column <- factor(qc_long$column,
  levels = c("CK\nRaw", "FCR\nRaw", "CK\nAdjusted", "FCR\nAdjusted"))
qc_long$metric <- factor(qc_long$metric,
  levels = rev(c("RNA counts", "Detected genes", "Mitochondrial RNA (%)", "Chloroplast RNA (%)")))
p_qc <- ggplot(qc_long, aes(column, metric, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", rho)), family = "Arial", size = 2.45,
            colour = ifelse(abs(qc_long$rho) >= 0.48, "white", "#222222")) +
  scale_fill_gradient2(low = "#355C8A", mid = "#F5F2EA", high = "#9C3D4A",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation\ncoefficient") +
  labs(x = NULL, y = NULL, title = "Association with sequencing-quality covariates") +
  theme_pc() +
  theme(axis.text.x = element_text(size = 6.3, lineheight = 0.95),
        axis.text.y = element_text(size = 6.5), axis.ticks = element_blank(),
        axis.line = element_blank(), legend.position = "right") +
  coord_fixed(ratio = 0.78)

# C: stress-overlap audit.
stress_df <- read.delim(file.path(source_dir, "02_stress_overlap_source.tsv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
stress_df$category <- factor(stress_df$category, levels = rev(stress_df$category))
stress_cols <- c(
  "Wounding / mechanical" = "#B05A4A", "Heat / temperature" = "#D58B3D",
  "Hypoxia / low oxygen" = "#5C78A8", "Oxidative / ROS" = "#8D568C",
  "Pathogen / immunity" = "#3F8775", "No selected stress GO annotation" = "#B8B8B8"
)
p_stress <- ggplot(stress_df, aes(fraction, category, fill = category)) +
  geom_col(width = 0.64, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d  (%.1f%%)", n_genes, 100 * fraction)),
            hjust = -0.08, family = "Arial", size = 2.35, colour = "#222222") +
  scale_fill_manual(values = stress_cols, guide = "none") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.28))) +
  labs(x = "Fraction of conserved Arabidopsis orthogroups", y = NULL,
       title = "Overlap with stress-associated GO processes") +
  theme_pc() +
  theme(axis.text.y = element_text(size = 6.5), axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

# D: matched-resolution clustering stability.
flow <- read.delim(file.path(source_dir, "03_cluster_transition_source.tsv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
tab <- xtabs(Cells ~ Original + Core_filtered, data = flow)
row_pct <- 100 * prop.table(tab, margin = 1)
old_index <- seq_len(nrow(row_pct))
new_centroid <- colSums(row_pct * old_index) / pmax(colSums(row_pct), 1e-9)
new_order <- colnames(row_pct)[order(new_centroid)]
stability <- as.data.frame(row_pct, stringsAsFactors = FALSE)
colnames(stability) <- c("Original", "Core_filtered", "Percent")
stability$Original <- factor(stability$Original, levels = rev(rownames(row_pct)))
stability$Core_filtered <- factor(stability$Core_filtered, levels = new_order)
stability$label <- ifelse(stability$Percent >= 5, sprintf("%.0f", stability$Percent), "")
p_flow <- ggplot(stability, aes(Core_filtered, Original, fill = Percent)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = label), family = "Arial", size = 2.15,
            colour = ifelse(stability$Percent >= 55, "white", "#222222")) +
  scale_fill_gradientn(colours = c("#F4F2EC", "#B9D4D0", "#4A8F86", "#174F55"),
                       limits = c(0, 100), name = "Cells from each\noriginal cluster (%)") +
  labs(x = "Core-filtered cluster", y = "Original cluster",
       title = "Clustering stability after core-signature exclusion") +
  theme_pc() +
  theme(axis.ticks = element_blank(), axis.line = element_blank(),
        axis.text = element_text(size = 6.2), legend.position = "right") +
  coord_fixed()

# E: equal-width subplots. Trim the connector at the open-circle boundary.
sens <- read.delim(file.path(source_dir, "04_treatment_contrast_sensitivity_source.tsv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
celltype_order <- c(
  "Epidermal cells I", "Epidermal cells II", "Epidermal cells III", "Guard cells",
  "Mesophyll cells I", "Mesophyll cells II", "Procambial cells", "Vascular parenchyma cells",
  "Phloem parenchyma cells", "Defense detoxification cells",
  "Defense cell wall remodeling cells I", "Defense cell wall remodeling cells II",
  "Defense JA responsive cells", "Proliferating S phase cells", "Translation active cells"
)
short_celltype <- setNames(c(
  "Epidermal I", "Epidermal II", "Epidermal III", "Guard", "Mesophyll I",
  "Mesophyll II", "Procambial", "Vascular parenchyma", "Phloem parenchyma",
  "Defense detoxification", "Defense wall remodeling I", "Defense wall remodeling II",
  "Defense JA responsive", "Proliferating S phase", "Translation active"
), celltype_order)
sens$celltype_label <- factor(unname(short_celltype[sens$celltype]),
                              levels = rev(unname(short_celltype[celltype_order])))

p_retained <- ggplot(sens, aes(retained_contrast_percent, celltype_label)) +
  geom_segment(aes(x = 100, xend = retained_contrast_percent, yend = celltype_label),
               colour = "#D8D8D8", linewidth = 0.7) +
  geom_point(aes(fill = top200_overlap_fraction), shape = 21, size = 2.4,
             colour = "#333333", stroke = 0.25) +
  geom_vline(xintercept = 100, linetype = "22", colour = "#777777", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", retained_contrast_percent)), hjust = 1.35,
            family = "Arial", size = 2.05, colour = "#222222") +
  scale_fill_gradient(low = "#D9E5EC", high = "#8D3F58",
                      labels = scales::percent_format(accuracy = 1),
                      name = "Core genes among\ntop 200 shifts") +
  guides(fill = guide_colourbar(
    title.position = "top", title.hjust = 0.5,
    barwidth = unit(40, "mm"), barheight = unit(4.2, "mm")
  )) +
  scale_x_continuous(limits = c(max(0, floor(min(sens$retained_contrast_percent) - 3)), 101),
                     breaks = scales::pretty_breaks(4)) +
  labs(x = "CK-FCR contrast magnitude retained (%)", y = NULL,
       title = "Treatment-contrast sensitivity") +
  theme_pc(7.1) +
  theme(axis.text.y = element_text(size = 5.9), legend.position = "bottom",
        plot.margin = margin(1.5, 1.5, 1.5, 0.5, unit = "mm"))

connector_gap <- 0.00065
sens$connector_start <- sens$expected_signature_fraction +
  sign(sens$top200_overlap_fraction - sens$expected_signature_fraction) * connector_gap
p_overlap <- ggplot(sens, aes(top200_overlap_fraction, celltype_label)) +
  geom_segment(aes(x = connector_start, xend = top200_overlap_fraction,
                   yend = celltype_label), colour = "#BFC4C8", linewidth = 0.65) +
  geom_point(colour = "#8D3F58", size = 2.1) +
  geom_point(aes(x = expected_signature_fraction), shape = 21, fill = "white",
             colour = "#555555", stroke = 0.45, size = 1.8) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "Core genes among top 200 shifts", y = NULL,
       title = "Response-gene overlap") +
  theme_pc(7.1) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = margin(1.5, 1.5, 1.5, 0.5, unit = "mm"))

p_sens <- p_retained + p_overlap +
  plot_layout(widths = c(1, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    legend.direction = "horizontal",
    # Optical centring against the complete E panel, including the wide
    # cell-type labels on the left. The compensating right margin shifts the
    # visible colour bar left while preserving the two plot columns.
    legend.box.margin = margin(0, 26, 0, 0, unit = "mm")
  )

save_panel_preview(p_qc, "panelB_refined.png", 112, 75)
save_panel_preview(p_stress, "panelC_refined.png", 118, 77)
save_panel_preview(p_flow, "panelD_refined.png", 118, 92)
save_panel_preview(p_sens, "panelE_refined.png", 180, 94)

save_panel_formal(p_qc, "Figure_protoplasting_QC_correlation", 112, 75)
save_panel_formal(p_stress, "Figure_protoplasting_stress_overlap", 118, 77)
save_panel_formal(p_flow, "Figure_protoplasting_cluster_stability", 118, 92)
save_panel_formal(p_sens, "Figure_protoplasting_treatment_sensitivity", 180, 94)

# Assemble the preview with equal half-width columns. Panel E spans both
# columns, and its two equal-width subplots therefore track the C/D columns.
read_panel <- function(path) {
  image_read(path) |>
    image_background("white", flatten = TRUE) |>
    image_trim(fuzz = 2) |>
    image_border("white", "8x8")
}
paths <- c(panel_a_path,
           file.path(output_dir, "Figure_protoplasting_QC_correlation.png"),
           file.path(output_dir, "Figure_protoplasting_stress_overlap.png"),
           file.path(output_dir, "Figure_protoplasting_cluster_stability.png"),
           file.path(output_dir, "Figure_protoplasting_treatment_sensitivity.png"))
panels <- lapply(paths, read_panel)

canvas_w <- 3000
margin_x <- 80
margin_y <- 70
gap_x <- 55
gap_y <- 55
label_h <- 62
half_w <- (canvas_w - 2 * margin_x - gap_x) / 2
row1_h <- 1030
row2_h <- 1030
# Give panel E enough height to fill the full two-column width. This removes
# the artificial side padding created by aspect-ratio fitting and aligns its
# two internal plots with the C/D column geometry above.
row3_h <- 1548
y1 <- margin_y
y2 <- y1 + row1_h + gap_y
y3 <- y2 + row2_h + gap_y
canvas_h <- y3 + row3_h + margin_y
slots <- data.frame(
  label = LETTERS[1:5],
  x = c(margin_x, margin_x + half_w + gap_x,
        margin_x, margin_x + half_w + gap_x, margin_x),
  y = c(y1, y1, y2, y2, y3),
  w = c(half_w, half_w, half_w, half_w, canvas_w - 2 * margin_x),
  h = c(row1_h, row1_h, row2_h, row2_h, row3_h)
)
fit_box <- function(img, x, y, w, h) {
  info <- image_info(img); iw <- info$width[[1]]; ih <- info$height[[1]]
  content_h <- h - label_h
  scale <- min(w / iw, content_h / ih)
  rw <- iw * scale; rh <- ih * scale
  c(x = x + (w - rw) / 2, y = y + label_h + (content_h - rh) / 2, w = rw, h = rh)
}
px <- function(v) unit(v / canvas_w, "npc")
py <- function(v) unit(1 - v / canvas_h, "npc")
pw <- function(v) unit(v / canvas_w, "npc")
ph <- function(v) unit(v / canvas_h, "npc")
draw_page <- function() {
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
  for (i in seq_along(panels)) {
    s <- slots[i, ]; b <- fit_box(panels[[i]], s$x, s$y, s$w, s$h)
    grid.raster(as.raster(panels[[i]]), x = px(b[["x"]]), y = py(b[["y"]]),
                width = pw(b[["w"]]), height = ph(b[["h"]]),
                just = c("left", "top"), interpolate = TRUE)
    grid.text(s$label, x = px(s$x + 4), y = py(s$y + 3), just = c("left", "top"),
              gp = gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8, col = "#111111"))
  }
}

preview_path <- file.path(output_dir, "Supplementary_Figure_protoplasting_impact_refined_preview.png")
agg_png(preview_path, width = 180, height = 180 * canvas_h / canvas_w,
        units = "mm", res = 300, background = "white")
draw_page()
dev.off()

formal_stem <- file.path(output_dir, "Supplementary_Figure_protoplasting_impact")
formal_width_mm <- 180
formal_height_mm <- formal_width_mm * canvas_h / canvas_w
grDevices::cairo_pdf(paste0(formal_stem, ".pdf"),
                     width = formal_width_mm / 25.4,
                     height = formal_height_mm / 25.4,
                     family = "Arial", onefile = TRUE)
draw_page(); dev.off()
ragg::agg_png(paste0(formal_stem, ".png"),
              width = formal_width_mm, height = formal_height_mm,
              units = "mm", res = 600, background = "white")
draw_page(); dev.off()
ragg::agg_tiff(paste0(formal_stem, ".tiff"),
               width = formal_width_mm, height = formal_height_mm,
               units = "mm", res = 600, compression = "lzw",
               background = "white")
draw_page(); dev.off()
cat("PREVIEW", preview_path, "\n")
