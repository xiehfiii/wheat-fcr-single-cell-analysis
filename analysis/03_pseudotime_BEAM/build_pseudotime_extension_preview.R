suppressPackageStartupMessages({
  library(monocle)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(Matrix)
  library(ragg)
})

set.seed(54)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else Sys.getenv("PSEUDOTIME_PROJECT_DIR", unset = ".")
out_dir <- if (length(args) >= 2) args[[2]] else file.path(base_dir, "extension_preview")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cds_path <- file.path(base_dir, "selected_preview", "Epidermal_cells_I_monocle2_selected.rds")
state_path <- file.path(base_dir, "selected_preview", "state_sample_composition.csv")
sensitivity_path <- file.path(base_dir, "sensitivity", "trajectory_sensitivity_metrics.csv")
scenic_path <- if (length(args) >= 3) args[[3]] else file.path(base_dir, "SCENIC_Summary.csv")

cds <- readRDS(cds_path)
meta <- as.data.frame(pData(cds))
meta$cell_id <- rownames(meta)
meta$Samples <- toupper(as.character(meta$Samples))

gene_erf <- "TraesCS4A02G001300"
gene_mgbp <- "TraesCS3D02G094200"
stopifnot(gene_erf %in% rownames(cds), gene_mgbp %in% rownames(cds))

theme_pc <- function(base_size = 8.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#252525"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#252525"),
      axis.text = element_text(colour = "#333333"),
      axis.title = element_text(colour = "#252525"),
      plot.title = element_text(face = "bold", size = base_size + 0.8, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.3, colour = "#555555", hjust = 0),
      strip.background = element_rect(fill = "#F2F2F0", colour = NA),
      strip.text = element_text(face = "bold", colour = "#303030", size = base_size - 0.2),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(5, 7, 5, 7)
    )
}

pal_sample <- c(CK = "#2E6F95", FCR = "#D9822B")
pal_branch <- c(
  "State 4 fate" = "#B75D3E",
  "State 5 fate" = "#4F7894",
  "State 2 fate" = "#B75D3E",
  "States 4/5 fate" = "#4F7894"
)

# A: state composition -------------------------------------------------------
state_df <- read.csv(state_path, stringsAsFactors = FALSE) %>%
  mutate(
    Samples = toupper(Samples),
    State = factor(State, levels = sort(unique(as.integer(State)))),
    label = ifelse(fraction >= 0.055, paste0(n, "\n", sprintf("%.1f%%", 100 * fraction)), "")
  )

p_state <- ggplot(state_df, aes(State, fraction, fill = Samples)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.35) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            size = 2.35, family = "Arial", colour = "white", lineheight = 0.9) +
  scale_fill_manual(values = pal_sample, breaks = c("CK", "FCR")) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "State composition by sample",
    subtitle = "Counts and within-State percentages",
    x = "Transcriptional State", y = "Cell fraction", fill = NULL
  ) +
  theme_pc() +
  theme(legend.position = "top", legend.justification = "left")

# B: sensitivity summary -----------------------------------------------------
sens <- read.csv(sensitivity_path, stringsAsFactors = FALSE) %>%
  filter(run %in% c("unsupervised_HVG1000", "unsupervised_HVG2000", "unsupervised_HVG2000_dim10")) %>%
  mutate(
    setting = recode(run,
      unsupervised_HVG1000 = "HVG 1,000\n15 dimensions",
      unsupervised_HVG2000 = "HVG 2,000\n15 dimensions",
      unsupervised_HVG2000_dim10 = "HVG 2,000\n10 dimensions"
    ),
    setting = factor(setting, levels = c(
      "HVG 1,000\n15 dimensions",
      "HVG 2,000\n15 dimensions",
      "HVG 2,000\n10 dimensions"
    ))
  )

sens_long <- sens %>%
  select(setting, states, branch_points) %>%
  pivot_longer(c(states, branch_points), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, states = "States", branch_points = "Branch points"))

p_sens <- ggplot(sens_long, aes(setting, value, colour = metric, group = metric)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_text(aes(label = value), vjust = -0.75, size = 2.6, family = "Arial", show.legend = FALSE) +
  geom_text(
    data = sens,
    aes(setting, 0.65, label = paste0("AUC ", sprintf("%.3f", fcr_vs_ck_pseudotime_auc))),
    inherit.aes = FALSE, size = 2.35, family = "Arial", colour = "#555555"
  ) +
  scale_colour_manual(values = c("States" = "#596F8A", "Branch points" = "#C26B45")) +
  scale_y_continuous(breaks = 0:7, limits = c(0.2, 7.8), expand = c(0, 0)) +
  labs(
    title = "Trajectory sensitivity",
    subtitle = "CK-to-FCR ordering is stable; fine topology varies",
    x = NULL, y = "Number", colour = NULL
  ) +
  theme_pc() +
  theme(
    legend.position = "top", legend.justification = "left",
    axis.text.x = element_text(size = 7.2, lineheight = 0.9)
  )

# C: target and regulon dynamics --------------------------------------------
scenic <- read.csv(scenic_path, stringsAsFactors = FALSE, check.names = FALSE)
target_string <- scenic$All_Targets[scenic$TF_ID == gene_erf][1]
regulon_targets <- unique(strsplit(target_string, ";", fixed = TRUE)[[1]])
regulon_targets <- setdiff(intersect(regulon_targets, rownames(cds)), gene_erf)
if (length(regulon_targets) < 5L) stop("Too few TaERF87 regulon targets were found in the trajectory object")

sf <- as.numeric(sizeFactors(cds))
lognorm <- function(gene_ids, object = cds) {
  x <- as.matrix(exprs(object)[gene_ids, , drop = FALSE])
  log1p(sweep(x, 2, as.numeric(sizeFactors(object)), "/"))
}

reg_mat <- lognorm(regulon_targets)
gene_sd <- apply(reg_mat, 1, sd)
reg_mat <- reg_mat[is.finite(gene_sd) & gene_sd > 0, , drop = FALSE]
reg_z <- t(scale(t(reg_mat)))
reg_score <- colMeans(reg_z, na.rm = TRUE)
mgbp_expr <- as.numeric(lognorm(gene_mgbp))

trend_df <- bind_rows(
  data.frame(cell_id = colnames(cds), Pseudotime = meta$Pseudotime,
             value = reg_score, feature = "TaERF87 regulon score"),
  data.frame(cell_id = colnames(cds), Pseudotime = meta$Pseudotime,
             value = mgbp_expr, feature = "TaMGBP1 expression")
) %>%
  mutate(feature = factor(feature, levels = c("TaERF87 regulon score", "TaMGBP1 expression")))

p_trend <- ggplot(trend_df, aes(Pseudotime, value)) +
  geom_point(size = 0.28, alpha = 0.16, colour = "#7F878D") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6), se = FALSE,
              linewidth = 1.05, colour = "#B75D3E") +
  facet_wrap(~feature, scales = "free_y", nrow = 1,
             labeller = as_labeller(c(
               "TaERF87 regulon score" = "italic('TaERF87')~regulon~score",
               "TaMGBP1 expression" = "italic('TaMGBP1')~expression"
             ), label_parsed)) +
  labs(
    title = "Regulatory module dynamics along pseudotime",
    subtitle = paste0("Descriptive GAM; regulon score summarizes ", nrow(reg_mat), " SCENIC target genes"),
    x = "Pseudotime", y = NULL
  ) +
  theme_pc() +
  theme(strip.text = element_text(face = "plain"))

# D: branch-specific fate dynamics ------------------------------------------
extract_branch <- function(bp) {
  bc <- buildBranchCellDataSet(
    cds,
    branch_point = bp,
    progenitor_method = "duplicate",
    branch_labels = c("Branch_A", "Branch_B")
  )
  bm <- as.data.frame(pData(bc))
  vals <- lognorm(c(gene_erf, gene_mgbp), object = bc)
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
  geom_point(size = 0.22, alpha = 0.08) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = FALSE, linewidth = 0.9) +
  facet_grid(branch_point ~ gene, scales = "free_y",
             labeller = labeller(gene = as_labeller(
               c(TaERF87 = "italic('TaERF87')", TaMGBP1 = "italic('TaMGBP1')"), label_parsed
             ))) +
  scale_colour_manual(values = pal_branch) +
  labs(
    title = "Fate-dependent expression dynamics",
    subtitle = "Curves are descriptive; progenitor cells are duplicated by Monocle2 for branch fitting",
    x = "Branched pseudotime", y = "log1p normalized expression", colour = NULL
  ) +
  theme_pc() +
  theme(
    legend.position = "top", legend.justification = "left",
    strip.text.x = element_text(face = "plain")
  )

# Assemble -------------------------------------------------------------------
fig <- ((p_state | p_sens) + plot_layout(widths = c(1, 1.12))) /
  p_trend /
  p_branch +
  plot_layout(heights = c(0.9, 0.82, 1.65)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12))

preview_path <- file.path(out_dir, "Pseudotime_extension_analyses_preview.png")
agg_png(preview_path, width = 2400, height = 2900, res = 300, background = "white")
print(fig)
dev.off()

write.csv(state_df, file.path(out_dir, "panelA_state_sample_composition.csv"), row.names = FALSE)
write.csv(sens, file.path(out_dir, "panelB_trajectory_sensitivity.csv"), row.names = FALSE)
write.csv(trend_df, file.path(out_dir, "panelC_regulon_target_pseudotime.csv"), row.names = FALSE)
write.csv(branch_df, file.path(out_dir, "panelD_branch_gene_dynamics.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

cat("PREVIEW=", preview_path, "\n", sep = "")
cat("CELLS=", ncol(cds), "\n", sep = "")
cat("REGULON_TARGETS_USED=", nrow(reg_mat), "\n", sep = "")
cat("STATE_COUNTS=", paste(tapply(rep(1, nrow(meta)), meta$State, sum), collapse = ","), "\n", sep = "")
