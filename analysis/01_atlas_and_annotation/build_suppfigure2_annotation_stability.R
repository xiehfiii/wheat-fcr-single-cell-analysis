suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
})

set.seed(54)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else Sys.getenv("SCRNASEQ_PROJECT_DIR", unset = ".")
out_dir <- if (length(args) >= 2) args[[2]] else file.path(base_dir, "Supplementary_Figure2_annotation_stability")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

merged_file <- file.path(base_dir, "merged.rds")
top10_file <- file.path(base_dir, "merged_top10_markers.tsv")
evidence_file <- file.path(base_dir, "marker_specificity_revision_20260909", "08_final_marker_panel.tsv")

merged <- readRDS(merged_file)
meta <- merged[[]]
stopifnot(nrow(meta) == ncol(merged), all(as.character(meta$seurat_clusters) == as.character(meta$SCT_snn_res.0.6)))

cluster_ids <- as.character(0:14)
cluster_key <- meta %>%
  transmute(cluster = as.character(seurat_clusters), celltype = as.character(celltype_final)) %>%
  distinct() %>%
  arrange(as.integer(cluster))
stopifnot(nrow(cluster_key) == 15L, identical(cluster_key$cluster, cluster_ids))

celltype_levels <- cluster_key$celltype
names(celltype_levels) <- cluster_key$cluster

theme_pc <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.30, colour = "#282828"),
      axis.ticks = element_line(linewidth = 0.28, colour = "#282828"),
      axis.text = element_text(size = base_size - 0.8, colour = "#303030"),
      axis.title = element_text(size = base_size, colour = "#252525"),
      legend.title = element_text(size = base_size - 0.4),
      legend.text = element_text(size = base_size - 0.7),
      panel.grid = element_blank(),
      plot.margin = margin(3, 4, 3, 4)
    )
}

# Panel A: exact top-10 marker gene set from the supplied heatmap, summarized
# as cluster-level mean log-normalized expression for compact readability.
top10 <- read.delim(top10_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(cluster = as.character(cluster)) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  mutate(marker_rank = row_number()) %>%
  ungroup()

stopifnot(nrow(top10) == 150L, all(top10$cluster %in% cluster_ids))
present_genes <- intersect(top10$gene, rownames(merged[["RNA"]]))
if (length(present_genes) != nrow(top10)) {
  stop("Missing ", nrow(top10) - length(present_genes), " top-marker genes from the RNA assay")
}

rna_data <- GetAssayData(merged, assay = "RNA", layer = "data")[top10$gene, , drop = FALSE]
cluster_factor <- factor(as.character(meta$seurat_clusters), levels = cluster_ids)
avg_mat <- vapply(cluster_ids, function(cl) {
  Matrix::rowMeans(rna_data[, cluster_factor == cl, drop = FALSE])
}, numeric(nrow(rna_data)))
rownames(avg_mat) <- rownames(rna_data)
colnames(avg_mat) <- cluster_ids

row_sd <- apply(avg_mat, 1, sd)
row_z <- (avg_mat - rowMeans(avg_mat)) / ifelse(row_sd > 0, row_sd, 1)
row_z <- pmax(pmin(row_z, 2.5), -2.5)

marker_order <- top10 %>%
  arrange(as.integer(cluster), marker_rank) %>%
  mutate(row_index = rev(seq_len(n())))

heat_df <- as.data.frame(row_z) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "display_cluster", values_to = "z") %>%
  left_join(marker_order %>% select(gene, marker_cluster = cluster, marker_rank, row_index), by = "gene") %>%
  mutate(
    display_cluster = factor(display_cluster, levels = cluster_ids),
    row_index = as.numeric(row_index)
  )

block_centres <- marker_order %>%
  group_by(cluster) %>%
  summarise(row_index = mean(row_index), .groups = "drop") %>%
  mutate(label = paste0(cluster, "  ", celltype_levels[cluster]))
block_boundaries <- marker_order %>%
  group_by(cluster) %>%
  summarise(y = min(row_index) - 0.5, .groups = "drop") %>%
  filter(y > 0)

p_a <- ggplot(heat_df, aes(display_cluster, row_index, fill = z)) +
  geom_raster() +
  geom_hline(data = block_boundaries, aes(yintercept = y), colour = "white", linewidth = 0.36) +
  scale_fill_gradient2(
    low = "#3C7FA1", mid = "#F3F2ED", high = "#B65B4B",
    midpoint = 0, limits = c(-2.5, 2.5), name = "Row z-score"
  ) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_continuous(
    breaks = block_centres$row_index, labels = block_centres$label,
    expand = c(0, 0)
  ) +
  labs(x = "Cluster", y = "Top-10 marker block", tag = "A") +
  theme_minimal(base_size = 6.2, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.x.top = element_text(size = 5.5, colour = "#303030", margin = margin(b = 1)),
    axis.text.y = element_text(size = 5.0, colour = "#303030", lineheight = 0.88),
    axis.title.x = element_text(size = 6.0, margin = margin(t = 2)),
    axis.title.y = element_text(size = 6.0),
    legend.position = "right", legend.title = element_text(size = 5.5),
    legend.text = element_text(size = 5.1), legend.key.height = grid::unit(11, "mm"),
    plot.margin = margin(4, 4, 3, 4)
  )

# Panel B: selected annotation markers grouped by evidence origin.
evidence <- read.delim(evidence_file, stringsAsFactors = FALSE, check.names = FALSE)
evidence_expanded <- evidence %>%
  separate_rows(target_population, sep = ";\\s*") %>%
  mutate(
    evidence_class = case_when(
      evidence_tier %in% c(
        "Gold/reference marker from wheat",
        "Direct wheat reference signature",
        "Spatially supported wheat marker"
      ) ~ "Wheat reference",
      evidence_tier == "Cross-species canonical gold marker mapped to wheat" ~ "Cross-species RBH",
      evidence_tier == "Reference marker; not strict gold standard" ~ "Other reference",
      TRUE ~ "Dataset-specific"
    ),
    target_population = factor(target_population, levels = rev(celltype_levels)),
    evidence_class = factor(
      evidence_class,
      levels = c("Wheat reference", "Cross-species RBH", "Other reference", "Dataset-specific")
    )
  )

evidence_counts <- evidence_expanded %>%
  count(target_population, evidence_class, name = "n_markers") %>%
  complete(target_population, evidence_class, fill = list(n_markers = 0))

p_b <- ggplot(evidence_counts, aes(evidence_class, target_population, fill = n_markers)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = ifelse(n_markers > 0, n_markers, "")), size = 2.0, family = "Arial") +
  scale_fill_gradient(low = "#F0EEE8", high = "#356F77", name = "Selected\nmarkers") +
  labs(x = NULL, y = NULL, tag = "B") +
  theme_pc(6.2) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1, size = 5.4),
    axis.text.y = element_text(size = 5.0),
    legend.position = "right", legend.key.height = grid::unit(10, "mm")
  )

# Panels C-D: cluster stability around the selected resolution 0.6.
choose2 <- function(x) x * (x - 1) / 2
adjusted_rand <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  sum_nij <- sum(choose2(tab))
  sum_ai <- sum(choose2(rowSums(tab)))
  sum_bj <- sum(choose2(colSums(tab)))
  expected <- sum_ai * sum_bj / choose2(n)
  denom <- 0.5 * (sum_ai + sum_bj) - expected
  if (denom == 0) return(1)
  (sum_nij - expected) / denom
}

normalized_mi <- function(x, y) {
  tab <- table(x, y)
  pxy <- tab / sum(tab)
  px <- rowSums(pxy)
  py <- colSums(pxy)
  nz <- which(pxy > 0, arr.ind = TRUE)
  mi <- sum(vapply(seq_len(nrow(nz)), function(i) {
    r <- nz[i, 1]; c <- nz[i, 2]
    pxy[r, c] * log(pxy[r, c] / (px[r] * py[c]))
  }, numeric(1)))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 && hy == 0) return(1)
  mi / sqrt(hx * hy)
}

resolution_values <- seq(0.3, 0.9, by = 0.1)
resolution_cols <- paste0("SCT_snn_res.", resolution_values)
reference <- as.character(meta$SCT_snn_res.0.6)

stability_summary <- bind_rows(lapply(seq_along(resolution_values), function(i) {
  alt <- as.character(meta[[resolution_cols[i]]])
  data.frame(
    resolution = resolution_values[i],
    n_clusters = length(unique(alt)),
    ARI = adjusted_rand(reference, alt),
    NMI = normalized_mi(reference, alt)
  )
}))

stability_long <- stability_summary %>%
  pivot_longer(c(ARI, NMI), names_to = "metric", values_to = "value")

p_d <- ggplot(stability_long, aes(resolution, value, colour = metric, group = metric)) +
  annotate("rect", xmin = 0.575, xmax = 0.625, ymin = -Inf, ymax = Inf,
           fill = "#E9E5DC", alpha = 0.65) +
  geom_line(linewidth = 0.70) +
  geom_point(size = 2.0) +
  geom_text(
    data = stability_summary,
    aes(resolution, 0.52, label = paste0("k=", n_clusters)),
    inherit.aes = FALSE, family = "Arial", size = 1.85, colour = "#666666"
  ) +
  scale_colour_manual(values = c(ARI = "#B85F44", NMI = "#426F8B")) +
  scale_x_continuous(breaks = resolution_values) +
  scale_y_continuous(limits = c(0.48, 1.03), breaks = seq(0.5, 1, 0.1)) +
  labs(x = "Clustering resolution", y = "Agreement with resolution 0.6", colour = NULL, tag = "D") +
  theme_pc(6.4) +
  theme(legend.position = "top", legend.justification = "left")

f1_rows <- list()
for (res in c(0.4, 0.5, 0.7, 0.8)) {
  alt <- as.character(meta[[paste0("SCT_snn_res.", res)]])
  tab <- table(reference, alt)
  for (cl in cluster_ids) {
    ai <- sum(tab[cl, ])
    f1 <- vapply(colnames(tab), function(a) {
      nij <- tab[cl, a]
      bj <- sum(tab[, a])
      if ((ai + bj) == 0) return(0)
      2 * nij / (ai + bj)
    }, numeric(1))
    best <- which.max(f1)
    f1_rows[[length(f1_rows) + 1L]] <- data.frame(
      cluster = cl,
      celltype = celltype_levels[cl],
      resolution = res,
      best_match = names(f1)[best],
      best_f1 = f1[best]
    )
  }
}
f1_df <- bind_rows(f1_rows) %>%
  mutate(
    celltype = factor(celltype, levels = rev(celltype_levels)),
    resolution_label = factor(
      paste0("res. ", resolution),
      levels = paste0("res. ", c(0.4, 0.5, 0.7, 0.8))
    )
  )

p_c <- ggplot(f1_df, aes(resolution_label, celltype, fill = best_f1)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%.2f", best_f1)), family = "Arial", size = 1.72) +
  scale_fill_gradient(low = "#ECE9E1", high = "#33727A", limits = c(0.45, 1),
                      oob = scales::squish, name = "Best-match F1") +
  labs(x = "Alternative resolution", y = NULL, tag = "C") +
  theme_pc(6.2) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 5.5), axis.text.y = element_text(size = 5.0),
    axis.title.x = element_text(margin = margin(t = -4, unit = "pt")),
    legend.position = "right", legend.key.height = grid::unit(10, "mm")
  )

# Panel E: descriptive representation of the captured cell populations in the
# two libraries. This is deliberately descriptive because there is one library
# per condition and therefore no replicate-supported abundance inference.
palette_by_cluster <- c(
  "0" = "#EE6677", "1" = "#4477AA", "2" = "#228833",
  "3" = "#66CCEE", "4" = "#EE7733", "5" = "#AA3377",
  "6" = "#44AA99", "7" = "#332288", "8" = "#999933",
  "9" = "#882255", "10" = "#117777", "11" = "#CCBB44",
  "12" = "#DDCC77", "13" = "#666666", "14" = "#AA4499"
)
palette_by_celltype <- setNames(
  unname(palette_by_cluster[cluster_ids]),
  unname(celltype_levels[cluster_ids])
)

composition_meta <- meta %>%
  transmute(
    sample = factor(toupper(as.character(Samples)), levels = c("CK", "FCR")),
    cluster = factor(as.character(seurat_clusters), levels = cluster_ids)
  )
stopifnot(!anyNA(composition_meta$sample), !anyNA(composition_meta$cluster))

composition_counts <- table(composition_meta$sample, composition_meta$cluster)
composition_props <- sweep(composition_counts, 1, rowSums(composition_counts), "/")
composition_df <- bind_rows(lapply(seq_len(nrow(composition_props)), function(i) {
  prop <- as.numeric(composition_props[i, ])
  data.frame(
    sample = rownames(composition_props)[i],
    x = i,
    cluster = cluster_ids,
    celltype = unname(celltype_levels[cluster_ids]),
    n = as.integer(composition_counts[i, ]),
    proportion = prop,
    ymin = c(0, head(cumsum(prop), -1)),
    ymax = cumsum(prop),
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(celltype = factor(celltype, levels = unname(celltype_levels[cluster_ids])))

composition_ribbons <- bind_rows(lapply(cluster_ids, function(cl) {
  ck <- composition_df %>% filter(sample == "CK", cluster == cl)
  fcr <- composition_df %>% filter(sample == "FCR", cluster == cl)
  data.frame(
    cluster = cl,
    celltype = ck$celltype,
    x = c(1.25, 1.75, 1.75, 1.25),
    y = c(ck$ymin, fcr$ymin, fcr$ymax, ck$ymax),
    vertex = seq_len(4)
  )
}))

composition_boundaries <- data.frame(
  y_ck = c(0, cumsum(as.numeric(composition_props["CK", ]))),
  y_fcr = c(0, cumsum(as.numeric(composition_props["FCR", ])))
)

p_e <- ggplot() +
  geom_polygon(
    data = composition_ribbons,
    aes(x = x, y = y, group = cluster, fill = celltype),
    alpha = 0.16, colour = NA
  ) +
  geom_segment(
    data = composition_boundaries,
    aes(x = 1.25, xend = 1.75, y = y_ck, yend = y_fcr),
    colour = "#555555", linewidth = 0.25, alpha = 0.65
  ) +
  geom_rect(
    data = composition_df,
    aes(xmin = x - 0.25, xmax = x + 0.25, ymin = ymin, ymax = ymax, fill = celltype),
    colour = "#303030", linewidth = 0.22
  ) +
  scale_fill_manual(values = palette_by_celltype, drop = FALSE) +
  scale_x_continuous(breaks = 1:2, labels = c("CK", "FCR"), limits = c(0.65, 2.35), expand = c(0, 0)) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), "%"),
    limits = c(0, 1.01), expand = c(0, 0)
  ) +
  labs(
    x = NULL, y = "Cell fraction", fill = NULL, tag = "E",
    title = NULL,
    subtitle = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_pc(6.2) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(
    plot.title = element_text(size = 6.3, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 5.2, colour = "#666666", hjust = 0),
    axis.text.x = element_text(face = "bold", size = 5.6),
    axis.text.y = element_text(size = 5.0),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "vertical",
    legend.text = element_text(size = 5.0),
    legend.key.height = grid::unit(2.35, "mm"),
    legend.key.width = grid::unit(2.6, "mm"),
    legend.spacing.y = grid::unit(0.15, "mm"),
    legend.box.spacing = grid::unit(1.2, "mm"),
    plot.margin = margin(3, 2, 3, 2, unit = "pt")
  )

figure <- p_a / (p_b | p_c) / (p_d | p_e) +
  plot_layout(heights = c(1.25, 0.98, 0.82), widths = 1) &
  theme(
    plot.tag = element_text(family = "Arial", face = "bold", size = 8, colour = "#111111"),
    plot.tag.position = c(0.004, 0.996)
  )

preview_file <- file.path(out_dir, "Supplementary_Figure2_annotation_evidence_stability_preview.png")
agg_png(preview_file, width = 2161, height = 2717, res = 300, background = "white")
print(figure)
dev.off()

# Formal journal exports. The assembled figure retains A-D tags; standalone
# panels intentionally omit them so they can be reused in later layouts.
save_pdf_png_tiff <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  grDevices::cairo_pdf(
    file.path(out_dir, paste0(stem, ".pdf")),
    width = width_in, height = height_in, family = "Arial", onefile = TRUE
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    file.path(out_dir, paste0(stem, ".png")),
    width = width_in, height = height_in, units = "in", res = dpi,
    background = "white"
  )
  print(plot)
  dev.off()

  ragg::agg_tiff(
    file.path(out_dir, paste0(stem, ".tiff")),
    width = width_in, height = height_in, units = "in", res = dpi,
    background = "white", compression = "lzw"
  )
  print(plot)
  dev.off()
}

p_a_standalone <- p_a + labs(tag = NULL)
p_b_standalone <- p_b + labs(tag = NULL)
p_c_standalone <- p_c + labs(tag = NULL)
p_d_standalone <- p_d + labs(tag = NULL)
p_e_standalone <- p_e + labs(tag = NULL)

save_pdf_png_tiff(
  figure, "Supplementary_Figure2_annotation_evidence_stability",
  width_mm = 183, height_mm = 230
)
save_pdf_png_tiff(
  p_a_standalone, "Supplementary_Figure2A_top10_marker_heatmap",
  width_mm = 183, height_mm = 96
)
save_pdf_png_tiff(
  p_b_standalone, "Supplementary_Figure2B_annotation_evidence",
  width_mm = 90, height_mm = 105
)
save_pdf_png_tiff(
  p_c_standalone, "Supplementary_Figure2C_cluster_stability_F1",
  width_mm = 90, height_mm = 105
)
save_pdf_png_tiff(
  p_d_standalone, "Supplementary_Figure2D_resolution_sensitivity",
  width_mm = 90, height_mm = 78
)
save_pdf_png_tiff(
  p_e_standalone, "Supplementary_Figure2E_descriptive_cell_representation",
  width_mm = 93, height_mm = 92
)

write.csv(cluster_key, file.path(out_dir, "cluster_celltype_key.csv"), row.names = FALSE)
write.csv(top10, file.path(out_dir, "panelA_top10_marker_set.csv"), row.names = FALSE)
write.csv(heat_df, file.path(out_dir, "panelA_cluster_mean_expression_source.csv"), row.names = FALSE)
write.csv(evidence_counts, file.path(out_dir, "panelB_annotation_evidence_counts.csv"), row.names = FALSE)
write.csv(stability_summary, file.path(out_dir, "panelD_resolution_stability_summary.csv"), row.names = FALSE)
write.csv(f1_df, file.path(out_dir, "panelC_cluster_best_match_F1.csv"), row.names = FALSE)
write.csv(composition_df, file.path(out_dir, "panelE_descriptive_cell_representation.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

cat("PREVIEW=", preview_file, "\n", sep = "")
cat("CELLS=", ncol(merged), "\n", sep = "")
cat("TOP_MARKERS=", nrow(top10), "\n", sep = "")
cat("FORMAL_EXPORTS=18\n")
