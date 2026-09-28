args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript build_figure6_cell_communication_optimized.R HANDOFF_DIRECTORY OUTPUT_DIRECTORY")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(igraph)
  library(ggraph)
  library(openxlsx)
})

set.seed(54)
handoff_dir <- normalizePath(args[[1]], mustWork = TRUE)
data_dir <- file.path(handoff_dir, "data")
provenance_dir <- file.path(handoff_dir, "provenance")
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
tag_family <- "Times New Roman"
col_ck <- "#3F6F9F"
col_fcr <- "#C46A4A"
col_neutral <- "#A9AAA6"
col_teal <- "#147D83"
col_pale <- "#F4F4F1"
col_low <- "#4F79A7"
col_high <- "#C66B4E"

read_panel <- function(name) {
  read_tsv(file.path(data_dir, name), show_col_types = FALSE)
}

panel_a_data <- read_panel("Figure6A_pathway_source.tsv.gz")
panel_b_data <- read_panel("Figure6B_expression_source.tsv.gz")
network_edges <- read_panel("Figure6C_network_source.tsv.gz")
panel_de_data <- read_panel("Figure6DE_core_routes_source.tsv.gz")
core_routes <- read_panel("Figure6_core_route_audit.tsv.gz")
mapping <- read_tsv(
  file.path(provenance_dir, "celltype_annotation_mapping.tsv"),
  show_col_types = FALSE
)
cell_counts <- read.xlsx(
  file.path(data_dir, "wheat_CellChat_final_validated.xlsx"),
  sheet = "cell_counts", rowNames = TRUE
)

old_col <- intersect(c("old_celltype", "old_name", "original_celltype"), names(mapping))[[1]]
new_col <- intersect(c("new_celltype", "final_celltype", "submission_celltype"), names(mapping))[[1]]
name_map <- setNames(mapping[[new_col]], mapping[[old_col]])

old_order <- c(
  "Defense-detoxification cells", "Defense-thickening cells I",
  "Defense-thickening cells II", "Epidermis cells I", "Epidermis cells II",
  "Guard cells", "JA-mediated defense cells", "Mesophyll cells I",
  "Mesophyll cells II", "Phloem cells", "Proliferating-S cells",
  "Protein-synthesis cells", "Sclerenchyma cells", "Trichome cells",
  "Wounding-responsive cells"
)
stopifnot(all(old_order %in% names(name_map)))
cell_order <- unname(name_map[old_order])
cell_index <- setNames(seq_along(cell_order), cell_order)

rename_identity <- function(x) {
  out <- unname(name_map[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}

panel_b_data <- panel_b_data %>%
  mutate(celltype_legacy = celltype, celltype = rename_identity(celltype))
network_edges <- network_edges %>%
  mutate(
    from_legacy = from, to_legacy = to,
    from = rename_identity(from), to = rename_identity(to),
    sender = rename_identity(sender)
  )
panel_de_data <- panel_de_data %>%
  mutate(
    source_legacy = source, target_legacy = target,
    source = rename_identity(source), target = rename_identity(target)
  )
core_routes <- core_routes %>%
  mutate(
    source_legacy = source, target_legacy = target,
    source = rename_identity(source), target = rename_identity(target)
  )

rownames(cell_counts) <- rename_identity(rownames(cell_counts))
stopifnot(all(cell_order %in% rownames(cell_counts)))

cell_colours_old <- c(
  "Defense-detoxification cells" = "#C94C4C",
  "Defense-thickening cells I" = "#3B78A8",
  "Defense-thickening cells II" = "#58A65C",
  "Epidermis cells I" = "#8E5AA6",
  "Epidermis cells II" = "#E28A22",
  "Guard cells" = "#D96BA5",
  "JA-mediated defense cells" = "#A67BB5",
  "Mesophyll cells I" = "#A65E2E",
  "Mesophyll cells II" = "#4AA3D8",
  "Phloem cells" = "#243B7B",
  "Proliferating-S cells" = "#1B9E77",
  "Protein-synthesis cells" = "#9CCB6B",
  "Sclerenchyma cells" = "#D8B400",
  "Trichome cells" = "#E88C8C",
  "Wounding-responsive cells" = "#D12A78"
)
cell_colours <- setNames(unname(cell_colours_old[old_order]), cell_order)

theme_pub <- function(base_size = 7.1) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "#292929"),
      axis.ticks = element_line(linewidth = 0.28, colour = "#292929"),
      axis.text = element_text(colour = "#242424"),
      axis.title = element_text(colour = "#242424"),
      plot.title = element_text(face = "bold", size = 8.2, hjust = 0),
      strip.text = element_text(face = "bold", size = 7.1),
      strip.background = element_rect(fill = "#F0F0ED", colour = NA),
      legend.title = element_text(size = 6.3),
      legend.text = element_text(size = 5.7),
      legend.key.height = grid::unit(3.0, "mm"),
      panel.grid = element_blank(),
      plot.margin = margin(3, 4, 3, 3)
    )
}

pathway_levels <- c(
  "DCC-DCCR", "RALF-FER", "RALF-THE", "CIF-SGN", "DEP-DEPR",
  "EPF_EPFL-ER", "TaFIP-TaFIPR"
)
panel_a_data <- panel_a_data %>%
  mutate(
    pathway_name = factor(pathway_name, levels = rev(pathway_levels)),
    value_label = sprintf("%+.3f", normalized_delta_mean_probability)
  )

p_a <- ggplot(panel_a_data, aes(y = pathway_name, x = normalized_delta_mean_probability)) +
  geom_vline(xintercept = 0, linewidth = 0.32, colour = "#777777") +
  geom_segment(
    aes(x = 0, xend = normalized_delta_mean_probability, yend = pathway_name,
        colour = direction_class),
    linewidth = 1.05, lineend = "round"
  ) +
  geom_point(aes(colour = direction_class), size = 2.8) +
  geom_text(
    aes(x = label_x, label = value_label, hjust = label_hjust),
    size = 2.05, family = font_family, colour = "#202020"
  ) +
  scale_colour_manual(values = c(
    "FCR higher" = col_fcr, "CK higher" = col_ck,
    "No retained score" = col_neutral
  )) +
  scale_x_continuous(
    limits = c(-0.90, 0.80),
    breaks = c(-0.75, -0.25, 0, 0.25, 0.75),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "Pathway-level communication shift",
    x = "Normalized difference (FCR − CK)", y = NULL
  ) +
  theme_pub() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 6.0),
    axis.text.x = element_text(size = 5.8)
  )

dcc_gene_levels <- c("DCC1-A", "DCC1-B", "DCC1-D", "DCCR1")
panel_b_data <- panel_b_data %>%
  mutate(
    condition = factor(condition, levels = c("CK", "FCR")),
    gene = factor(gene, levels = dcc_gene_levels),
    celltype = factor(celltype, levels = rev(cell_order)),
    cell_id = factor(cell_index[as.character(celltype)], levels = rev(seq_along(cell_order)))
  )

p_b <- ggplot(
  panel_b_data,
  aes(x = gene, y = cell_id, size = fraction_expressing, colour = mean_LogNormalize)
) +
  geom_point(alpha = 0.96) +
  facet_grid(. ~ condition) +
  scale_y_discrete(drop = FALSE) +
  scale_size_continuous(
    range = c(0.25, 3.6), limits = c(0, 0.80),
    breaks = c(0.10, 0.30, 0.50, 0.70),
    labels = percent_format(accuracy = 1), name = "Cells expressing"
  ) +
  scale_colour_gradientn(
    colours = c("#F2F3F0", "#A9D4CF", col_teal, "#07555D"),
    values = rescale(c(0, 0.08, 0.35, 0.80)),
    limits = c(0, 0.80), oob = squish,
    name = "Mean expression\n(LogNormalize)"
  ) +
  labs(
    title = "Expression basis of the DCC1–DCCR1 module",
    x = NULL, y = "Cell-type index"
  ) +
  theme_pub() +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 5.8, face = "italic", angle = 35, hjust = 1),
    axis.text.y = element_text(size = 5.3),
    axis.title.y = element_text(size = 5.8, margin = margin(r = 2)),
    panel.spacing.x = grid::unit(2.5, "mm"),
    legend.position = "right",
    legend.box.spacing = grid::unit(1.0, "mm")
  ) +
  guides(
    size = guide_legend(order = 1, override.aes = list(colour = "#5E8F91")),
    colour = guide_colourbar(order = 2, barheight = grid::unit(16, "mm"),
                            barwidth = grid::unit(2.2, "mm"))
  )

angle_deg <- seq(90, 90 - 360 + 360 / length(cell_order), length.out = length(cell_order))
theta <- angle_deg * pi / 180
network_nodes <- tibble(
  name = cell_order,
  x = cos(theta), y = sin(theta),
  cell_id = seq_along(cell_order)
)

network_strength_max <- max(network_edges$strength, na.rm = TRUE)
node_count_max <- max(as.matrix(cell_counts[cell_order, c("ck", "fcr")]), na.rm = TRUE)

hex_luminance <- function(x) {
  rgb <- col2rgb(x) / 255
  0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
}
node_label_colours <- setNames(
  ifelse(hex_luminance(cell_colours) > 0.62, "#1B1B1B", "white"),
  names(cell_colours)
)

network_panel <- function(condition_key) {
  condition_column <- tolower(condition_key)
  nodes_use <- network_nodes %>%
    mutate(
      cells = as.numeric(cell_counts[name, condition_column]),
      label_colour = unname(node_label_colours[name])
    )
  edges_use <- network_edges %>%
    filter(condition == condition_key) %>%
    select(from, to, everything(), -condition)
  graph_use <- graph_from_data_frame(edges_use, directed = TRUE, vertices = nodes_use)

  ggraph(graph_use, layout = "manual", x = x, y = y) +
    geom_edge_arc(
      aes(filter = !self_route, width = strength, edge_colour = sender),
      strength = 0.18, alpha = 0.32,
      arrow = grid::arrow(length = grid::unit(0.9, "mm"), type = "closed"),
      start_cap = circle(2.0, "mm"), end_cap = circle(2.2, "mm"),
      lineend = "round", show.legend = FALSE
    ) +
    geom_edge_loop(
      aes(filter = self_route, width = strength, edge_colour = sender),
      alpha = 0.34,
      arrow = grid::arrow(length = grid::unit(0.9, "mm"), type = "closed"),
      start_cap = circle(2.0, "mm"), end_cap = circle(2.2, "mm"),
      lineend = "round", show.legend = FALSE
    ) +
    geom_node_point(
      aes(size = cells, fill = name), shape = 21,
      colour = "white", stroke = 0.45, show.legend = FALSE
    ) +
    geom_node_text(
      aes(label = cell_id, colour = label_colour),
      family = font_family, fontface = "bold", size = 2.0,
      show.legend = FALSE
    ) +
    scale_colour_identity() +
    scale_edge_width(range = c(0.18, 1.25), limits = c(0, network_strength_max)) +
    scale_edge_colour_manual(values = cell_colours) +
    scale_size_continuous(range = c(3.0, 6.0), limits = c(0, node_count_max)) +
    scale_fill_manual(values = cell_colours) +
    coord_fixed(xlim = c(-1.32, 1.32), ylim = c(-1.28, 1.28), clip = "off") +
    labs(title = condition_key) +
    theme_void(base_family = font_family) +
    theme(
      plot.title = element_text(face = "bold", size = 7.2, hjust = 0.5,
                                margin = margin(b = 1)),
      plot.margin = margin(1, 2, 1, 2)
    )
}

key_data <- tibble(
  celltype = cell_order,
  cell_id = seq_along(cell_order),
  column = rep(1:3, each = 5),
  row = rep(5:1, times = 3),
  label_colour = unname(node_label_colours[cell_order])
)

p_key <- ggplot(key_data, aes(x = column, y = row)) +
  geom_point(aes(fill = celltype), shape = 21, size = 3.0,
             colour = "white", stroke = 0.35, show.legend = FALSE) +
  geom_text(aes(label = cell_id, colour = label_colour), family = font_family,
            fontface = "bold", size = 1.80, show.legend = FALSE) +
  geom_text(aes(x = column + 0.07, label = celltype), hjust = 0,
            family = font_family, size = 1.80, colour = "#202020") +
  scale_fill_manual(values = cell_colours) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0.93, 3.98), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 5.5), expand = c(0, 0)) +
  labs(title = "Cell-type index (shared across B–E)") +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(face = "bold", size = 6.3, hjust = 0),
    plot.margin = margin(0, 3, 1, 3)
  )

p_c_inner <- (
  (network_panel("CK") | network_panel("FCR")) /
    p_key +
    plot_layout(heights = c(1.0, 0.46)) +
    plot_annotation(title = "Complete retained family-level communication networks",
                    theme = theme(plot.title = element_text(
                      family = font_family, face = "bold", size = 8.2, hjust = 0.5,
                      margin = margin(b = 1))))
)
p_c <- wrap_elements(full = p_c_inner)

panel_de_data <- panel_de_data %>%
  mutate(
    source_id = factor(cell_index[source], levels = rev(seq_along(cell_order))),
    target_id = factor(cell_index[target], levels = seq_along(cell_order))
  )
route_background <- expand_grid(
  source_id = factor(rev(seq_along(cell_order)), levels = rev(seq_along(cell_order))),
  target_id = factor(seq_along(cell_order), levels = seq_along(cell_order))
)

route_panel <- function(data, title_text, show_legend = TRUE) {
  ggplot() +
    geom_tile(
      data = route_background, aes(x = target_id, y = source_id),
      fill = col_pale, colour = "white", linewidth = 0.18
    ) +
    geom_point(
      data = data,
      aes(x = target_id, y = source_id, size = n_core_routes,
          colour = median_normalized_delta), alpha = 0.94
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    scale_size_continuous(
      range = c(1.5, 4.2), limits = c(1, 3), breaks = 1:3,
      name = "Core LR families"
    ) +
    scale_colour_gradient2(
      low = col_low, mid = "#F0F0EC", high = col_high,
      midpoint = 0, limits = c(-0.75, 0.75), oob = squish,
      name = "Median normalized\ndifference"
    ) +
    labs(title = title_text, x = "Receiver cell-type index", y = "Sender cell-type index") +
    theme_pub() +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(size = 5.2),
      axis.text.y = element_text(size = 5.2),
      axis.title = element_text(size = 6.2),
      legend.position = if (show_legend) "right" else "none",
      panel.grid = element_blank(),
      plot.margin = margin(3, 4, 2, 3)
    )
}

p_d <- route_panel(
  filter(panel_de_data, pathway_name == "DCC-DCCR"),
  "DCC1–DCCR1 high-confidence core routes (n = 48)", FALSE
)
p_e <- route_panel(
  filter(panel_de_data, pathway_name == "RALF-FER"),
  "RALF–FER high-confidence core routes (n = 62)", TRUE
)

figure6 <- (
  (p_a | p_b) / p_c / p_d / p_e +
    plot_layout(heights = c(0.88, 1.18, 0.84, 0.88)) +
    plot_annotation(tag_levels = "A")
) &
  theme(
    plot.tag = element_text(
      family = tag_family, face = "plain", size = 12,
      colour = "#111111"
    ),
    plot.tag.position = c(0.004, 0.996)
  )

width_mm <- 183
height_mm <- 230
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
stem <- file.path(outdir, "Figure6_cell_communication_optimized")

svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
print(figure6)
dev.off()

grDevices::cairo_pdf(
  paste0(stem, ".pdf"), width = width_in, height = height_in,
  family = font_family
)
print(figure6)
dev.off()

ragg::agg_tiff(
  paste0(stem, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw"
)
print(figure6)
dev.off()

ragg::agg_png(
  paste0(stem, ".png"), width = width_in, height = height_in,
  units = "in", res = 600
)
print(figure6)
dev.off()

write_tsv(panel_a_data, file.path(outdir, "Figure6A_pathway_source.tsv.gz"))
write_tsv(panel_b_data, file.path(outdir, "Figure6B_expression_source.tsv.gz"))
write_tsv(network_edges, file.path(outdir, "Figure6C_network_source.tsv.gz"))
write_tsv(panel_de_data, file.path(outdir, "Figure6DE_core_routes_source.tsv.gz"))
write_tsv(core_routes, file.path(outdir, "Figure6_core_route_audit.tsv.gz"))

summary_out <- tibble(
  metric = c(
    "DCC_core_routes", "RALF_FER_core_routes",
    "DCC_min_direction_stability", "RALF_FER_min_direction_stability",
    "CK_retained_network_edges", "FCR_retained_network_edges",
    "Panel_A_rows", "Panel_B_rows", "Panel_C_rows",
    "Panel_DE_rows", "Core_route_audit_rows"
  ),
  value = c(
    sum(core_routes$pathway_name == "DCC-DCCR"),
    sum(core_routes$pathway_name == "RALF-FER"),
    min(core_routes$direction_stability[core_routes$pathway_name == "DCC-DCCR"]),
    min(core_routes$direction_stability[core_routes$pathway_name == "RALF-FER"]),
    sum(network_edges$condition == "CK"),
    sum(network_edges$condition == "FCR"),
    nrow(panel_a_data), nrow(panel_b_data), nrow(network_edges),
    nrow(panel_de_data), nrow(core_routes)
  )
)
write_tsv(summary_out, file.path(outdir, "Figure6_QA_summary.tsv"))

legend_text <- paste0(
  "Figure 6 | Cell-type-resolved remodeling of peptide–receptor communication in wheat leaf sheath following Fusarium crown rot challenge. ",
  "(A) Pathway-level normalized differences in mean family communication probability between FCR and CK. Positive values indicate higher normalized communication in FCR, whereas negative values indicate higher communication in CK. ",
  "(B) Expression basis of the DCC1–DCCR1 module, showing mean LogNormalize expression (colour) and the fraction of expressing cells (dot size) for DCC1-A, DCC1-B, DCC1-D and DCCR1 in each annotated cell population and condition. ",
  "(C) Complete retained family-level communication networks in CK and FCR. Both conditions use identical circular node coordinates and a shared edge-width scale; node size represents cell number, node colour denotes cell identity, edge colour denotes the sender, and arrows indicate communication direction. All 88 CK edges and 55 FCR edges retained by the validated workflow are shown. ",
  "(D,E) High-confidence core routes for DCC1–DCCR1 (D; n = 48) and the candidate RALF–FER module (E; n = 62). Dot size indicates the number of core ligand–receptor families and colour indicates the median normalized difference, using a shared zero-centred scale. ",
  "Numbers 1–15 identify the same cell-type order across panels B–E. These comparisons are descriptive and do not constitute treatment-effect hypothesis tests: cells were not treated as biological replicates, and cell-resampling stability is not a biological confidence interval. ",
  "The DCC1–DCCR1 assignment is supported by direct wheat mature-peptide–receptor evidence, without implying separate binding validation for every A/B/D precursor gene. The RALF–FER module is a wheat candidate inferred by stringent orthology mapping from experimental evidence in other species and has not been directly validated in wheat."
)
writeLines(legend_text, file.path(outdir, "Figure6_Legend_EN.md"), useBytes = TRUE)

change_log <- c(
  "# Figure 6 change log",
  "",
  "## Visual and naming changes",
  "- Replaced legacy cell-type labels with the final submission labels from provenance/celltype_annotation_mapping.tsv.",
  "- Introduced a shared 1–15 cell-type index to keep every full label legible at 183 mm width.",
  "- Used identical circular node coordinates and shared node/edge scales for CK and FCR.",
  "- Reduced network-edge opacity, refined curvature and arrowheads, and moved the full cell-identity key below panel C.",
  "- Applied a restrained colour-blind-aware blue–neutral–warm scale shared by panels D and E.",
  "- Standardized Arial internal typography and 12 pt Times New Roman panel letters.",
  "",
  "## Analysis changes",
  "- None. No CellChat result was re-screened, no threshold or statistic was changed, and no route or network edge was added or removed.",
  "- Panel A, B, C and D/E values were read directly from the frozen panel source tables in the handoff package.",
  "- Random seed was fixed at 54; node placement is deterministic."
)
writeLines(change_log, file.path(outdir, "Figure6_CHANGELOG.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(outdir, "R_sessionInfo.txt"))

cat("FIGURE6_OPTIMIZED_BUILD_OK\n")
