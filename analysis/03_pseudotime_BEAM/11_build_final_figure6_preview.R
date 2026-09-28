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
if (length(args) != 4L) {
  stop("Usage: Rscript 11_build_final_figure6_preview.R SELECTED_CDS_RDS HEATMAP_COMPONENT_RDS GO_DISPLAY_CSV OUTPUT_DIRECTORY")
}

set.seed(54)
cds_file <- normalizePath(args[[1]], mustWork = TRUE)
heatmap_file <- normalizePath(args[[2]], mustWork = TRUE)
go_file <- normalizePath(args[[3]], mustWork = TRUE)
outdir <- args[[4]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

font_body <- "Arial"
font_tag <- "Times New Roman"
col_ck <- "#4C78A8"
col_fcr <- "#C96F55"
col_state4 <- "#B8664F"
col_state5 <- "#4F7698"
col_pre <- "#C9C9C7"
state_colours <- c("#456A95", "#63A58C", "#D2A24C", "#B46A87", "#7B6FA8")

cds <- readRDS(cds_file)
comp <- readRDS(heatmap_file)
go <- read.csv(go_file, stringsAsFactors = FALSE, check.names = FALSE)
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
  scale_colour_gradientn(colours = c("#EDF1F4", "#9DBFD1", "#4D86AB", "#183D63"), name = "Pseudotime") +
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
  ) + guides(fill = "none", colour = guide_legend(override.aes = list(fill = NA, linewidth = 0.9, linetype = 1)))

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

mat <- comp$matrix
tree_order <- rownames(mat)[comp$tree$order]
module_key <- data.frame(
  gene_id = names(comp$modules),
  module = paste0("M", as.integer(comp$modules)),
  stringsAsFactors = FALSE
)
module_labels <- c(M1 = "M1", M2 = "M2")
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
  segment = factor(comp$annotation_col[[1]],
                   levels = c("State 4 terminal", "Pre-branch", "State 5 terminal"))
)
segment_centres <- col_anno %>% group_by(segment) %>% summarise(x = mean(matrix_index), .groups = "drop")

p_colanno <- ggplot(col_anno, aes(matrix_index, 1, fill = segment)) +
  geom_raster() +
  scale_fill_manual(values = c("State 4 terminal" = col_state4, "Pre-branch" = col_pre,
                               "State 5 terminal" = col_state5), guide = "none") +
  scale_x_continuous(position = "top", breaks = segment_centres$x,
                     labels = as.character(segment_centres$segment), expand = c(0, 0)) +
  coord_cartesian(expand = FALSE) +
  labs(tag = "D", title = "Branch point 1: terminal-state divergence") +
  theme_void(base_family = font_body) +
  theme(plot.title = element_text(size = 6.8, face = "bold", hjust = 0.5, margin = margin(b = 2)),
        axis.text.x.top = element_text(size = 5.4, colour = "#222222", margin = margin(b = 1)),
        plot.margin = margin(2, 2, 0, 1))

p_heat <- ggplot(heat, aes(matrix_index, gene_id, fill = z)) +
  geom_raster(interpolate = FALSE) +
  facet_grid(module_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = "#315E7A", mid = "#F4F2EC", high = "#B65C4B",
                       midpoint = 0, limits = c(-3, 3), oob = squish, name = "Row z-score") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 6.4, base_family = font_body) +
  theme(
    panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    strip.placement = "outside", strip.background = element_blank(),
    strip.text.y.left = element_text(size = 5.4, angle = 90, face = "bold"),
    panel.spacing.y = grid::unit(0.7, "mm"),
    legend.position = "right", legend.title = element_text(size = 5.4),
    legend.text = element_text(size = 5.0), legend.key.height = grid::unit(8, "mm"),
    plot.margin = margin(0, 2, 2, 1)
  )

p_heat_full <- p_colanno / p_heat + plot_layout(heights = c(0.085, 0.915))

ratio_to_num <- function(x) vapply(strsplit(x, "/", fixed = TRUE), function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1))
go <- go %>%
  mutate(
    gene_ratio = ratio_to_num(GeneRatio),
    neglog10_fdr = -log10(pmax(p.adjust, .Machine$double.xmin)),
    description_wrapped = vapply(Description, function(x) paste(strwrap(x, width = 31), collapse = "\n"), character(1)),
    module_label = factor(module_labels[module], levels = module_labels)
  ) %>%
  arrange(module_label, desc(neglog10_fdr)) %>%
  mutate(term_order = factor(paste(module, ID), levels = rev(paste(module, ID))))

p_go <- ggplot(go, aes(gene_ratio, term_order, size = Count, colour = neglog10_fdr)) +
  geom_point(alpha = 0.92) +
  facet_grid(module_label ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = setNames(go$description_wrapped, paste(go$module, go$ID))) +
  scale_size_continuous(range = c(1.7, 4.2), name = "Gene count") +
  scale_colour_gradient(low = "#B9C7CE", high = "#304F68", name = expression(-log[10]~FDR)) +
  labs(x = "Gene ratio", y = NULL, title = "GO biological processes") +
  theme_classic(base_size = 6.4, base_family = font_body) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.30), axis.ticks.x = element_line(linewidth = 0.28),
    axis.text.x = element_text(size = 5.1), axis.text.y = element_text(size = 5.0, lineheight = 0.94),
    axis.title.x = element_text(size = 5.8), plot.title = element_text(size = 6.8, face = "bold", hjust = 0.5),
    strip.background = element_blank(), strip.text = element_text(size = 5.3, face = "bold"),
    panel.spacing.y = grid::unit(1.2, "mm"),
    legend.position = "right", legend.title = element_text(size = 5.2),
    legend.text = element_text(size = 4.9),
    plot.margin = margin(5, 2, 2, 4)
  )

top <- p_state | p_sample | p_time
middle_gene <- (p_gene | p_gene_facet) + plot_layout(widths = c(0.38, 0.62))
middle <- (p_density | middle_gene) + plot_layout(widths = c(0.30, 0.70))
bottom <- (p_heat_full | p_go) + plot_layout(widths = c(0.55, 0.45))
figure <- top / middle / bottom + plot_layout(heights = c(0.78, 0.75, 1.22)) &
  theme(
    plot.tag = element_text(family = font_tag, face = "plain", size = 12, colour = "#111111"),
    plot.tag.position = c(0.005, 0.995)
  )

ragg::agg_png(file.path(outdir, "Figure6_pseudotime_final_preview.png"),
              width = 183 / 25.4, height = 205 / 25.4, units = "in", res = 300, background = "white")
print(figure)
dev.off()

write.table(source_cells, gzfile(file.path(outdir, "Figure6_source_cells.tsv.gz")),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(heat %>% select(gene_id, module, matrix_index, z),
          file.path(outdir, "Figure6D_heatmap_source.csv"), row.names = FALSE)
write.csv(go, file.path(outdir, "Figure6D_GO_display_source.csv"), row.names = FALSE)
write.csv(medians, file.path(outdir, "Figure6B_pseudotime_medians.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo_final_figure.txt"))
