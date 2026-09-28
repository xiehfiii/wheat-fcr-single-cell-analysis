suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(magick)
  library(grid)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: Rscript script.R volcano.png overlap.csv scenic_summary.csv knockout.csv rss.pdf tf_umap.png output.png")
}

volcano_path <- args[[1]]
overlap_path <- args[[2]]
scenic_path <- args[[3]]
ko_path <- args[[4]]
rss_path <- args[[5]]
tf_umap_path <- args[[6]]
out_path <- args[[7]]

gene_erf <- "TraesCS4A02G001300"
gene_mgbp <- "TraesCS3D02G094200"
col_erf <- "#6256C7"
col_mgbp <- "#D55E5E"
col_warm <- "#D36B50"
col_blue <- "#4D78A8"
col_dark <- "#1F2933"
col_mid <- "#72808E"
col_light <- "#D8DEE4"

theme_set(
  theme_classic(base_family = "Arial", base_size = 7) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = col_dark),
      axis.ticks = element_line(linewidth = 0.35, colour = col_dark),
      plot.title = element_text(size = 7.2, face = "bold", colour = col_dark,
                                hjust = 0.5, margin = margin(b = 3)),
      axis.title = element_text(size = 6.5, colour = col_dark),
      axis.text = element_text(size = 5.7, colour = "#3E4A55")
    )
)

overlap <- read.csv(overlap_path, check.names = FALSE, stringsAsFactors = FALSE)
overlap <- overlap[!duplicated(overlap$Gene), , drop = FALSE]
overlap <- overlap[order(overlap$Induced_LogFC, decreasing = TRUE), , drop = FALSE]
overlap$Rank <- seq_len(nrow(overlap))
stopifnot(nrow(overlap) == 320)

scenic <- read.csv(scenic_path, check.names = FALSE, stringsAsFactors = FALSE)
erf_row <- scenic[scenic$TF_ID == gene_erf, , drop = FALSE]
stopifnot(nrow(erf_row) == 1)
targets <- unique(strsplit(erf_row$All_Targets[[1]], ";", fixed = TRUE)[[1]])
targets <- targets[nzchar(targets)]
stopifnot(length(targets) == 44, gene_mgbp %in% targets)

ko <- read.csv(ko_path, check.names = FALSE, stringsAsFactors = FALSE)
ko_gene_col <- names(ko)[2]
ko_sig <- ko[is.finite(ko$Z) & ko$Z > 2 & is.finite(ko$p.adj) & ko$p.adj < 0.05, , drop = FALSE]
ko_genes <- unique(as.character(ko_sig[[ko_gene_col]]))
shared <- sort(intersect(targets, ko_genes))
stopifnot(length(ko_genes) == 137, length(shared) == 5, gene_mgbp %in% shared)
targets_downstream <- setdiff(targets, gene_erf)
ko_downstream <- setdiff(ko_genes, gene_erf)
shared_downstream <- sort(intersect(targets_downstream, ko_downstream))
stopifnot(length(targets_downstream) == 43, length(ko_downstream) == 136,
          length(shared_downstream) == 4, gene_mgbp %in% shared_downstream)

# A: preserve the original volcano and add restrained direct labels only.
volcano_img <- image_read(volcano_path)[1] |> image_background("white", flatten = TRUE)
vi <- image_info(volcano_img)
vw <- vi$width[[1]]
vh <- vi$height[[1]]

pA <- ggplot() +
  annotation_raster(as.raster(volcano_img), xmin = 0, xmax = vw, ymin = 0, ymax = vh) +
  coord_fixed(ratio = 1, xlim = c(0, vw), ylim = c(0, vh), expand = FALSE, clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(2, 2, 1, 2))

# B: the intersection alone establishes the FCR-induced epidermal marker set.
circle_df <- function(cx, cy, rx, ry, set_name) {
  a <- seq(0, 2 * pi, length.out = 361)
  data.frame(x = cx + rx * cos(a), y = cy + ry * sin(a), set = set_name)
}
circles <- rbind(
  circle_df(0.40, 0.48, 0.34, 0.36, "Marker"),
  circle_df(0.70, 0.48, 0.34, 0.36, "Up")
)
pB <- ggplot(circles, aes(x, y, group = set, fill = set, colour = set)) +
  geom_polygon(alpha = 0.25, linewidth = 0.65) +
  annotate("text", x = 0.06, y = 0.94, label = "Epidermal cells I markers\n(n = 1,808)",
           family = "Arial", size = 2.55, fontface = "bold", lineheight = 0.94,
           colour = "#334155") +
  annotate("text", x = 1.04, y = 0.94, label = "FCR-upregulated genes\n(n = 1,556)",
           family = "Arial", size = 2.55, fontface = "bold", lineheight = 0.94,
           colour = "#334155") +
  annotate("label", x = 0.55, y = 0.51, label = "320",
           family = "Arial", size = 4.0, fontface = "bold", linewidth = 0,
           label.padding = unit(0.18, "lines"), fill = "white", colour = col_dark) +
  annotate("text", x = 0.55, y = 0.37, label = "FCR-induced marker\ngenes",
           family = "Arial", size = 1.95, fontface = "bold", lineheight = 0.90,
           colour = "#52606D") +
  scale_fill_manual(values = c(Marker = "#4C78A8", Up = "#E07A5F")) +
  scale_colour_manual(values = c(Marker = "#365F8C", Up = "#B85C46")) +
  coord_equal(xlim = c(-0.24, 1.34), ylim = c(0.01, 1.03), clip = "off") +
  labs(title = NULL) +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 7.2, face = "bold", hjust = 0.5,
                                  colour = col_dark, margin = margin(b = 2)),
        legend.position = "none", plot.margin = margin(4, 4, 2, 4))

# C: rank the 320-gene set, then show the specific SCENIC candidate edge.
overlap$Target <- ifelse(overlap$Gene == gene_mgbp, "TaMGBP1",
                         ifelse(overlap$Gene == gene_erf, "TaERF87", "Other"))
labC <- overlap[overlap$Target != "Other", , drop = FALSE]
labC$display <- paste0("italic(", labC$Target, ")~'rank ", labC$Rank, "'")
pC1 <- ggplot(overlap, aes(Rank, Induced_LogFC)) +
  geom_line(linewidth = 0.36, colour = "#B9C1C9") +
  geom_point(aes(colour = Target), size = 0.82, alpha = 0.82) +
  geom_point(data = labC, aes(colour = Target), size = 1.45, shape = 21,
             fill = "white", stroke = 0.48) +
  geom_text_repel(data = labC, aes(label = display, colour = Target), parse = TRUE,
                  family = "Arial", size = 2.15, seed = 54,
                  box.padding = 0.25, point.padding = 0.22, min.segment.length = 0,
                  segment.size = 0.25, segment.colour = "#6B7580", max.overlaps = Inf) +
  scale_colour_manual(values = c(Other = "#B9C1C9", TaMGBP1 = col_mgbp, TaERF87 = col_erf)) +
  scale_x_continuous(limits = c(-15, max(overlap$Rank) + 8), expand = c(0, 0)) +
  labs(x = "Rank among 320 genes", y = expression(log[2] * " fold change")) +
  theme(legend.position = "none", plot.margin = margin(5, 6, 2, 2),
        axis.title = element_text(size = 5.8), axis.text = element_text(size = 5.2))

pC2 <- ggplot() +
  annotate("label", x = 0.40, y = 0.76, label = "italic(TaERF87)~'rank 13'", parse = TRUE,
           family = "Arial", size = 2.7, lineheight = 0.95,
           fill = "white", colour = col_erf, linewidth = 0.65,
           label.padding = unit(0.20, "lines")) +
  annotate("segment", x = 0.40, xend = 0.40, y = 0.63, yend = 0.38,
           colour = "#7A8793", linewidth = 0.65,
           arrow = arrow(length = unit(1.7, "mm"), type = "closed")) +
  annotate("text", x = 0.56, y = 0.51, hjust = 0, label = "SCENIC\ncandidate\nregulation",
           family = "Arial", fontface = "bold", size = 1.70, lineheight = 0.90,
           colour = "#56616D") +
  annotate("label", x = 0.40, y = 0.24, label = "italic(TaMGBP1)~'rank 2'", parse = TRUE,
           family = "Arial", size = 2.7, lineheight = 0.95,
           fill = "white", colour = col_mgbp, linewidth = 0.65,
           label.padding = unit(0.20, "lines")) +
  coord_cartesian(xlim = c(0.10, 0.98), ylim = c(0.08, 0.94), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(8, 4, 4, 4))

pC <- (pC1 | pC2) + plot_layout(widths = c(1.25, 0.75))

# D: compact rendering of the original scTenifoldKnk perturbation dot plot.
# The perturbed TF itself is excluded from the downstream response display.
ko_plot <- ko[ko[[ko_gene_col]] != gene_erf & is.finite(ko$Z) & is.finite(ko$p.value), , drop = FALSE]
ko_plot$neglog10p <- -log10(pmax(ko_plot$p.value, 1e-30))
ko_plot$status <- ifelse(ko_plot$Z > 2 & ko_plot$p.adj < 0.05, "Significant", "Other")
ko_plot$shared <- ko_plot[[ko_gene_col]] %in% shared_downstream
ko_mgbp <- ko_plot[ko_plot[[ko_gene_col]] == gene_mgbp, , drop = FALSE]

pK <- ggplot(ko_plot, aes(Z, neglog10p)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.32,
             colour = "#C76A62") +
  geom_vline(xintercept = 2, linetype = "dashed", linewidth = 0.32,
             colour = "#5578A7") +
  geom_point(aes(colour = status), size = 0.42, alpha = 0.48) +
  geom_point(data = ko_mgbp, shape = 21, size = 1.40,
             stroke = 0.45, fill = "white", colour = col_mgbp) +
  annotate("segment", x = ko_mgbp$Z[[1]] + 0.10, y = ko_mgbp$neglog10p[[1]] + 0.15,
           xend = 2.94, yend = 5.78, linewidth = 0.25, colour = col_mgbp) +
  annotate("text", x = 3.22, y = 5.42, hjust = 0.5,
           label = "italic(TaMGBP1)", parse = TRUE,
           family = "Arial", size = 1.45, colour = col_mgbp) +
  scale_colour_manual(values = c(Other = "#B9C1C9", Significant = "#263746")) +
  scale_x_continuous(limits = c(-4.2, 4.2), breaks = c(-4, -2, 0, 2, 4),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 30.5), breaks = c(0, 10, 20, 30),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Z-score", y = expression(-log[10] * "(p value)")) +
  theme(legend.position = "none",
        axis.title = element_text(size = 5.5),
        axis.text = element_text(size = 5.0),
        plot.margin = margin(5, 4, 3, 4))

# E: use an UpSet-style display so this validation step cannot be confused
# with the marker/upregulated Venn diagram in B.
upset_counts <- data.frame(
  group = factor(c("SCENIC only", "Shared", "Knockout only"),
                 levels = c("SCENIC only", "Shared", "Knockout only")),
  n = c(39, 4, 132),
  type = c("SCENIC", "Shared", "KO")
)
pDbar <- ggplot(upset_counts, aes(group, n, fill = type)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.35) +
  geom_text(aes(label = n), vjust = -0.35, family = "Arial", fontface = "bold",
            size = 2.45, colour = col_dark) +
  scale_fill_manual(values = c(SCENIC = "#6F86B5", Shared = col_mgbp, KO = "#D58B61")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = "Genes") +
  theme_classic(base_family = "Arial", base_size = 6.2) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 5.5), axis.text.y = element_text(size = 5.2),
        legend.position = "none", plot.margin = margin(2, 2, 0, 2))

matrix_df <- data.frame(
  group = rep(levels(upset_counts$group), each = 2),
  set = rep(c("SCENIC", "scTenifoldKnk"), times = 3),
  present = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE)
)
pDmat <- ggplot(matrix_df, aes(group, set)) +
  geom_segment(data = data.frame(x = "Shared", xend = "Shared",
                                 y = "SCENIC", yend = "scTenifoldKnk"),
               aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
               linewidth = 0.55, colour = col_dark) +
  geom_point(aes(fill = present), shape = 21, size = 2.35, stroke = 0.4,
             colour = "#5C6874") +
  scale_fill_manual(values = c(`TRUE` = col_dark, `FALSE` = "#E4E8EC")) +
  scale_x_discrete(limits = c("SCENIC only", "Shared", "Knockout only"),
                   labels = c("SCENIC\nonly", "Shared", "Knockout\nonly")) +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = "Arial") +
  theme(axis.text.x = element_text(size = 5.1, colour = col_dark, margin = margin(t = 1)),
        axis.text.y = element_text(size = 5.1, colour = col_dark, hjust = 1,
                                   margin = margin(r = 2)),
        legend.position = "none", plot.margin = margin(0, 2, 2, 2))

pDupset <- pDbar / pDmat + plot_layout(heights = c(0.72, 0.28))
shared_card <- data.frame(
  label = c("TraesCS1B02G276500", "TaMGBP1", "TraesCS5B02G449500", "TraesCS6A02G018200"),
  y = c(0.76, 0.59, 0.42, 0.25),
  type = c("Other", "TaMGBP1", "Other", "Other")
)
pDcard <- ggplot(shared_card, aes(x = 0.5, y = y)) +
  geom_point(aes(x = 0.08, colour = type), size = 1.7) +
  geom_text(data = subset(shared_card, type == "Other"),
            aes(x = 0.15, label = label, colour = type), hjust = 0,
            family = "Arial", size = 1.76) +
  annotate("text", x = 0.15, y = 0.59, hjust = 0, label = "italic(TaMGBP1)", parse = TRUE,
           family = "Arial", size = 1.80, colour = col_mgbp) +
  scale_colour_manual(values = c(Other = "#65717D", TaMGBP1 = col_mgbp)) +
  coord_cartesian(xlim = c(0, 1.18), ylim = c(0.16, 0.84), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(legend.position = "none", plot.margin = margin(4, 5, 4, 5))

pD <- (pDupset | pDcard) + plot_layout(widths = c(1.10, 0.90))

# F: keep the original UMAP values and remove the embedded subplot titles.
tf_img <- image_read(tf_umap_path)[1] |> image_background("white", flatten = TRUE) |> image_trim(fuzz = 2)
ti <- image_info(tf_img)
tw <- ti$width[[1]]
th <- ti$height[[1]]
pE <- ggplot() +
  annotation_raster(as.raster(tf_img), xmin = 0, xmax = tw, ymin = 0, ymax = th) +
  coord_fixed(ratio = 1, xlim = c(0, tw), ylim = c(0, th * 0.94), expand = FALSE) +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(1, 2, 1, 2))

# G: end with TF filtering and the exact RSS heatmap.
funnel <- data.frame(
  y = c(3, 2, 1),
  label = c("FCR-responsive\nmarker genes", "PlantTFDB\ntranscription factors",
            "SCENIC\nregulons"),
  n = c(320, 9, 7),
  fill = c("#F3F6F8", "#EAF0F4", "#DFE9F0")
)
pD1 <- ggplot(funnel) +
  geom_label(aes(x = 2, y = y, label = paste0(label, "\n", n), fill = fill),
             family = "Arial", fontface = "bold", size = 2.25, lineheight = 0.9,
             colour = col_dark, linewidth = 0.35,
             label.padding = unit(0.18, "lines"), label.r = unit(0.12, "lines")) +
  scale_fill_identity() +
  annotate("segment", x = 2, xend = 2, y = 2.66, yend = 2.36,
           arrow = arrow(length = unit(1.1, "mm"), type = "closed"), linewidth = 0.42,
           colour = col_mid) +
  annotate("segment", x = 2, xend = 2, y = 1.66, yend = 1.36,
           arrow = arrow(length = unit(1.1, "mm"), type = "closed"), linewidth = 0.42,
           colour = col_mid) +
  coord_cartesian(xlim = c(0, 4), ylim = c(0.45, 3.55), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(6, 2, 5, 2))

rss_img <- image_read_pdf(rss_path, density = 300)[1] |>
  image_background("white", flatten = TRUE) |>
  image_trim(fuzz = 2)
pF2 <- wrap_elements(full = rasterGrob(as.raster(rss_img), interpolate = TRUE))
pF <- (pD1 | pF2) + plot_layout(widths = c(0.52, 1.48))

A_panel <- wrap_elements(full = rasterGrob(as.raster(volcano_img), interpolate = TRUE))
B_panel <- pB
C_panel <- wrap_elements(full = pC)
K_panel <- pK
D_panel <- wrap_elements(full = pD)
E_panel <- pE
F_panel <- wrap_elements(full = pF)

row1 <- A_panel + B_panel + plot_layout(widths = c(0.82, 1.18))
row2 <- C_panel + K_panel + D_panel + plot_layout(widths = c(1.00, 0.70, 1.30))

final_plot <- row1 / row2 / E_panel / F_panel +
  plot_layout(heights = c(1.00, 0.92, 1.20, 0.62)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.margin = margin(4, 4, 4, 4)
    )
  ) &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 8,
                                colour = "black"),
        plot.tag.position = c(0.005, 0.995))

out_base <- sub("\\.[^.]+$", "", out_path)

ragg::agg_png(paste0(out_base, ".png"), width = 183, height = 277,
              units = "mm", res = 600, background = "white")
print(final_plot)
dev.off()

grDevices::cairo_pdf(paste0(out_base, ".pdf"), width = 183 / 25.4,
                     height = 277 / 25.4, family = "Arial")
print(final_plot)
dev.off()

ragg::agg_tiff(paste0(out_base, ".tiff"), width = 183, height = 277,
               units = "mm", res = 600, compression = "lzw",
               background = "white")
print(final_plot)
dev.off()

if (identical(Sys.getenv("EXPORT_PANELS"), "1")) {
  save_panel_triplet <- function(plot_obj, filename, width_mm, height_mm) {
    ragg::agg_png(paste0(filename, ".png"), width = width_mm, height = height_mm,
                  units = "mm", res = 600, background = "white")
    print(plot_obj)
    dev.off()

    grDevices::cairo_pdf(paste0(filename, ".pdf"), width = width_mm / 25.4,
                         height = height_mm / 25.4, family = "Arial")
    print(plot_obj)
    dev.off()

    ragg::agg_tiff(paste0(filename, ".tiff"), width = width_mm, height = height_mm,
                   units = "mm", res = 600, compression = "lzw",
                   background = "white")
    print(plot_obj)
    dev.off()
  }

  panel_specs <- list(
    list(A_panel, "Figure2_panelA_volcano", 89, 89),
    list(B_panel, "Figure2_panelB_FCR_induced_marker_genes", 89, 75),
    list(C_panel, "Figure2_panelC_ranked_genes_SCENIC", 68, 75),
    list(K_panel, "Figure2_panelD_virtual_knockout", 68, 75),
    list(D_panel, "Figure2_panelE_SCENIC_knockout_intersection", 82, 75),
    list(E_panel, "Figure2_panelF_TaERF87_expression_regulon_UMAP", 183, 90),
    list(F_panel, "Figure2_panelG_TF_filtering_RSS_heatmap", 183, 58)
  )

  for (spec in panel_specs) {
    save_panel_triplet(spec[[1]], spec[[2]], spec[[3]], spec[[4]])
  }
}

cat("QA_COUNTS\n")
cat("overlap", nrow(overlap), "\n")
cat("TaMGBP1_rank", overlap$Rank[overlap$Gene == gene_mgbp], "\n")
cat("TaERF87_rank", overlap$Rank[overlap$Gene == gene_erf], "\n")
cat("SCENIC_targets", length(targets), "\n")
cat("KO_significant", length(ko_genes), "\n")
cat("shared", length(shared), paste(shared, collapse = ";"), "\n")
cat("downstream_shared", length(shared_downstream), paste(shared_downstream, collapse = ";"), "\n")
