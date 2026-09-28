############################################################
## Wheat CellChat: statistically bounded, no-replicate version
## Standalone validated run. This script starts in a separate Rscript process,
## reads the Seurat object from disk, and never clears or modifies the user's
## interactive R workspace.
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
})

merged <- readRDS("/home/xiehangfei/scRNAseq/merged.rds")

lr_pairs <- data.frame(
  ligand = c(
    "TraesCS1D02G057600",
    "TraesCS1B02G309600", "TraesCS3A02G346000", "TraesCS3D02G415500",
    "TraesCS4B02G277700", "TraesCS6B02G347500",
    "TraesCS5A02G093000", "TraesCS5B02G099100", "TraesCS5D02G105300",
    "TraesCS5A02G093000", "TraesCS5B02G099100", "TraesCS5D02G105300",
    "TraesCS5A02G093000", "TraesCS5B02G099100", "TraesCS5D02G105300",
    "TraesCS2A02G276500", "TraesCS2A02G276600",
    "TraesCS2B02G294200", "TraesCS2B02G294300",
    "TraesCS2D02G275500", "TraesCS2D02G275600",
    "TraesCS4A02G024100",
    "TraesCS4A02G024100", "TraesCS4A02G024200", "TraesCS4A02G024300",
    "TraesCS4B02G279200", "TraesCS4D02G277800", "TraesCS4D02G277900"
  ),
  receptor = c(
    "TraesCS2B02G251700",
    rep("TraesCS7D02G166100", 5),
    rep("TraesCS1D02G228900", 3),
    rep("TraesCS4A02G133800", 3),
    rep("TraesCS7B02G241900", 3),
    rep("TraesCS1A02G342600", 6),
    "TraesCS3B02G310000",
    rep("TraesCS3D02G276200", 6)
  ),
  pathway_name = c(
    "CIF-SGN",
    rep("EPF_EPFL-ER", 5),
    rep("RALF-FER", 6),
    rep("RALF-THE", 3),
    rep("DCC-DCCR", 6),
    "DEP-DEPR",
    rep("TaFIP-TaFIPR", 6)
  ),
  stringsAsFactors = FALSE
)

lr_pairs$interaction_name <- paste(
  lr_pairs$ligand, lr_pairs$receptor, sep = "_"
)
lr_pairs$agonist <- ""
lr_pairs$antagonist <- ""
lr_pairs$co_A_receptor <- ""
lr_pairs$co_I_receptor <- ""
lr_pairs$annotation <- "Secreted Signaling"
lr_pairs$interaction_name_2 <- paste(
  lr_pairs$ligand, lr_pairs$receptor, sep = " - "
)
lr_pairs <- unique(lr_pairs)

stopifnot(inherits(merged, "Seurat"))
stopifnot(is.data.frame(lr_pairs))
stopifnot(all(c("ligand", "receptor") %in% colnames(lr_pairs)))

master_seed       <- getOption("wheat_cellchat.seed", 54L)
sample_col        <- "Samples"
celltype_col      <- "celltype"
assay_use         <- "RNA"
condition_ck      <- "ck"
condition_fcr     <- "fcr"
min_cells_group   <- 10L
min_expression_pct <- 0.10
trim_use          <- 0.10
nboot_cellchat    <- 5000L
nboot_stability   <- 500L
stability_cutoff  <- 0.90
effect_cutoff     <- 0.20
outdir <- "/home/xiehangfei/scRNAseq/cell_comm_final_validated"

if (length(nboot_stability) != 1L || !is.finite(nboot_stability) ||
    nboot_stability < 1L || nboot_stability != as.integer(nboot_stability)) {
  stop("nboot_stability must be one positive integer.")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
options(stringsAsFactors = FALSE)

get_data_layer <- function(seu, assay) {
  tryCatch(
    GetAssayData(seu, assay = assay, layer = "data"),
    error = function(e) GetAssayData(seu, assay = assay, slot = "data")
  )
}

normalise_condition <- function(x) tolower(trimws(as.character(x)))

meta_all <- merged@meta.data
stopifnot(all(c(sample_col, celltype_col) %in% colnames(meta_all)))

meta_all$.condition <- normalise_condition(meta_all[[sample_col]])
meta_all$.celltype  <- as.character(meta_all[[celltype_col]])

keep_condition <- meta_all$.condition %in% c(condition_ck, condition_fcr)
keep_cell <- rownames(meta_all)[keep_condition & !is.na(meta_all$.celltype)]
seu_use <- merged[, keep_cell]
seu_use$.condition <- meta_all[keep_cell, ".condition"]
seu_use$.celltype  <- meta_all[keep_cell, ".celltype"]

count_table <- table(seu_use$.celltype, seu_use$.condition)
if (!all(c(condition_ck, condition_fcr) %in% colnames(count_table))) {
  stop("Both ck and fcr must be present in the sample column.")
}

common_celltypes <- rownames(count_table)[
  count_table[, condition_ck] >= min_cells_group &
    count_table[, condition_fcr] >= min_cells_group
]
if (length(common_celltypes) < 2L) {
  stop("Fewer than two cell types have at least min_cells_group cells in both conditions.")
}

seu_use <- subset(seu_use, subset = .celltype %in% common_celltypes)
cell_order <- sort(common_celltypes)

############################################################
## Prepare one common non-negative RNA LogNormalize layer
## on the local object only; the global `merged` is unchanged.
############################################################

assay_names <- names(seu_use@assays)

if (length(assay_use) != 1L || is.na(assay_use) || !is.character(assay_use)) {
  stop("assay_use must be one non-missing character value.")
}

if (!(assay_use %in% assay_names)) {
  stop(
    paste0(
      "The Seurat object does not contain the requested assay: ",
      assay_use,
      ". Available assays: ",
      paste(assay_names, collapse = ", ")
    )
  )
}

DefaultAssay(seu_use) <- assay_use
message("Selected assay: ", DefaultAssay(seu_use))

# A merged Seurat v5 object can retain sample-specific counts/data layers.
# Join them before normalization so CellChat reads one common matrix.
if (inherits(seu_use[[assay_use]], "Assay5")) {
  layer_names_before <- SeuratObject::Layers(seu_use[[assay_use]])
  message("RNA layers before JoinLayers: ", paste(layer_names_before, collapse = ", "))

  if (length(layer_names_before) > 1L) {
    seu_use <- SeuratObject::JoinLayers(seu_use, assay = assay_use)
  }

  layer_names_after <- SeuratObject::Layers(seu_use[[assay_use]])
  message("RNA layers after JoinLayers: ", paste(layer_names_after, collapse = ", "))
}

seu_use <- NormalizeData(
  seu_use,
  assay = assay_use,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

rna_data_check <- get_data_layer(seu_use, assay_use)
if (nrow(rna_data_check) == 0L || ncol(rna_data_check) == 0L) {
  stop("The normalized RNA data layer is empty.")
}
rna_nonzero_values <- if (inherits(rna_data_check, "sparseMatrix")) {
  rna_data_check@x
} else {
  as.numeric(rna_data_check)
}
if (any(!is.finite(rna_nonzero_values))) {
  stop("The normalized RNA data layer contains non-finite values.")
}
if (any(rna_nonzero_values < 0)) {
  stop("The normalized RNA data layer contains negative values.")
}

message("CellChat input: RNA assay, common LogNormalize data layer")
message("RNA matrix dimensions: ", nrow(rna_data_check), " genes x ", ncol(rna_data_check), " cells")

############################################################
## Prepare the hypothesis-driven LR database
############################################################

lr_db <- lr_pairs %>%
  mutate(
    ligand = as.character(ligand),
    receptor = as.character(receptor),
    pathway_name = if ("pathway_name" %in% names(.)) {
      as.character(pathway_name)
    } else if ("pathway" %in% names(.)) {
      as.character(pathway)
    } else {
      "Uncategorized"
    },
    interaction_name = if ("interaction_name" %in% names(.)) {
      make.unique(as.character(interaction_name))
    } else {
      make.unique(paste(ligand, receptor, sep = "_"))
    }
  ) %>%
  filter(!is.na(ligand), !is.na(receptor), ligand != "", receptor != "") %>%
  distinct(interaction_name, .keep_all = TRUE)

dcc1_genes <- c(
  "TraesCS2A02G276600", "TraesCS2B02G294300", "TraesCS2D02G275600"
)
dcc2_genes <- c(
  "TraesCS2A02G276500", "TraesCS2B02G294200", "TraesCS2D02G275500"
)

ligand_family <- function(x) {
  case_when(
    x %in% dcc1_genes ~ "DCC1_family",
    x %in% dcc2_genes ~ "DCC2_family",
    TRUE ~ x
  )
}

lr_db <- lr_db %>%
  mutate(
    ligand_family = ligand_family(ligand),
    receptor_family = receptor,
    signal_family = paste(pathway_name, ligand_family, receptor_family, sep = "|")
  )

all_lr_genes <- unique(c(lr_db$ligand, lr_db$receptor))
expr_all <- get_data_layer(seu_use, assay_use)
missing_genes <- setdiff(all_lr_genes, rownames(expr_all))
if (length(missing_genes) > 0L) {
  warning("LR genes absent from the assay and removed: ", paste(missing_genes, collapse = ", "))
}

lr_db <- lr_db %>%
  filter(ligand %in% rownames(expr_all), receptor %in% rownames(expr_all))
if (nrow(lr_db) == 0L) {
  stop("None of the ligand-receptor pairs is present in the selected assay.")
}
all_lr_genes <- unique(c(lr_db$ligand, lr_db$receptor))

if (any(expr_all[all_lr_genes, , drop = FALSE] < 0)) {
  stop("The selected data layer contains negative values. Use a non-negative normalized data layer, not scale.data or integrated residuals.")
}

build_custom_cellchat_db <- function(lr_table, genes) {
  interaction <- lr_table %>%
    transmute(
      interaction_name = as.character(interaction_name),
      pathway_name = as.character(pathway_name),
      ligand = as.character(ligand),
      receptor = as.character(receptor),
      agonist = "",
      antagonist = "",
      co_A_receptor = "",
      co_I_receptor = "",
      evidence = "curated_direct_or_high_quality_mapping",
      annotation = "Secreted Signaling",
      interaction_name_2 = paste(ligand, receptor, sep = " - ")
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)
  rownames(interaction) <- make.unique(interaction$interaction_name)

  complex <- as.data.frame(
    matrix(character(0), nrow = 0, ncol = 4,
           dimnames = list(NULL, paste0("subunit_", 1:4))),
    stringsAsFactors = FALSE
  )
  cofactor <- as.data.frame(
    matrix(character(0), nrow = 0, ncol = 16,
           dimnames = list(NULL, paste0("cofactor", 1:16))),
    stringsAsFactors = FALSE
  )
  gene_info <- data.frame(
    Symbol = genes,
    Name = genes,
    EntrezGene.ID = NA_integer_,
    Ensembl.Gene.ID = "",
    MGI.ID = "",
    Gene.group.name = "custom_wheat_LR_gene",
    row.names = genes,
    stringsAsFactors = FALSE
  )

  list(
    interaction = interaction,
    complex = complex,
    cofactor = cofactor,
    geneInfo = gene_info
  )
}

CellChatDB_wheat <- build_custom_cellchat_db(lr_db, all_lr_genes)

############################################################
## Run CellChat separately within each observed sample
## These p/q values are within-sample screening quantities only.
############################################################

group_expression <- function(data_mat, labels, genes, groups) {
  avg <- matrix(0, nrow = length(groups), ncol = length(genes),
                dimnames = list(groups, genes))
  pct <- avg
  for (g in groups) {
    cells_g <- names(labels)[labels == g]
    x <- data_mat[genes, cells_g, drop = FALSE]
    avg[g, ] <- Matrix::rowMeans(x)
    pct[g, ] <- Matrix::rowMeans(x > 0)
  }
  list(avg = avg, pct = pct)
}

run_one_condition <- function(seu, condition_name) {
  cells_condition <- colnames(seu)[seu$.condition == condition_name]
  seu_sub <- seu[, cells_condition]
  data_input <- get_data_layer(seu_sub, assay_use)[all_lr_genes, , drop = FALSE]

  labels <- as.character(seu_sub$.celltype)
  names(labels) <- colnames(seu_sub)
  meta_input <- data.frame(
    labels = labels,
    samples = factor(rep(condition_name, length(labels))),
    row.names = names(labels),
    stringsAsFactors = FALSE
  )

  cellchat <- createCellChat(
    object = data_input,
    meta = meta_input,
    group.by = "labels"
  )
  cellchat@DB <- CellChatDB_wheat
  cellchat <- subsetData(cellchat)

  pair_use <- cellchat@DB$interaction[
    cellchat@DB$interaction$ligand %in% rownames(data_input) &
      cellchat@DB$interaction$receptor %in% rownames(data_input),
    , drop = FALSE
  ]
  if ("interaction_name" %in% colnames(pair_use)) {
    rownames(pair_use) <- make.unique(as.character(pair_use$interaction_name))
  }
  cellchat@LR$LRsig <- pair_use

  cellchat <- computeCommunProb(
    cellchat,
    type = "truncatedMean",
    trim = trim_use,
    raw.use = TRUE,
    population.size = FALSE,
    distance.use = FALSE,
    nboot = nboot_cellchat,
    seed.use = master_seed
  )
  cellchat <- filterCommunication(cellchat, min.cells = min_cells_group)

  prob_arr <- cellchat@net$prob
  pval_arr <- cellchat@net$pval
  sources <- dimnames(prob_arr)[[1]]
  targets <- dimnames(prob_arr)[[2]]
  interactions <- dimnames(prob_arr)[[3]]

  out <- expand.grid(
    source = sources,
    target = targets,
    interaction_name = interactions,
    stringsAsFactors = FALSE
  )
  out$prob <- as.vector(prob_arr)
  out$p_cellchat <- as.vector(pval_arr)

  out <- out %>%
    left_join(
      lr_db %>% select(interaction_name, ligand, receptor, pathway_name,
                       ligand_family, receptor_family, signal_family),
      by = "interaction_name"
    )

  expr_stats <- group_expression(data_input, labels, all_lr_genes, sources)
  source_i <- match(out$source, rownames(expr_stats$avg))
  target_i <- match(out$target, rownames(expr_stats$avg))
  ligand_i <- match(out$ligand, colnames(expr_stats$avg))
  receptor_i <- match(out$receptor, colnames(expr_stats$avg))

  out$ligand_mean <- expr_stats$avg[cbind(source_i, ligand_i)]
  out$receptor_mean <- expr_stats$avg[cbind(target_i, receptor_i)]
  out$ligand_pct <- expr_stats$pct[cbind(source_i, ligand_i)]
  out$receptor_pct <- expr_stats$pct[cbind(target_i, receptor_i)]

  # Finite-permutation correction: p can never be exactly zero.
  # This remains a within-sample screening p value, not a CK-FCR p value.
  exceedances <- round(pmin(1, pmax(0, out$p_cellchat)) * nboot_cellchat)
  out$p_empirical <- (exceedances + 1) / (nboot_cellchat + 1)
  out$q_within_sample <- p.adjust(out$p_empirical, method = "BH")

  out <- out %>%
    mutate(
      condition = condition_name,
      passes_expression = ligand_pct >= min_expression_pct &
        receptor_pct >= min_expression_pct,
      supported_within_sample = prob > 0 & passes_expression &
        q_within_sample <= 0.05
    ) %>%
    select(
      condition, source, target, interaction_name, pathway_name,
      ligand, receptor, ligand_family, receptor_family, signal_family,
      prob, p_cellchat, p_empirical, q_within_sample,
      ligand_pct, receptor_pct, ligand_mean, receptor_mean,
      passes_expression, supported_within_sample
    )

  list(cellchat = cellchat, table = out)
}

result_ck <- run_one_condition(seu_use, condition_ck)
result_fcr <- run_one_condition(seu_use, condition_fcr)

############################################################
## Descriptive CK-FCR comparison using RAW probabilities
## No differential p value is calculated because n=1 per condition.
############################################################

join_keys <- c(
  "source", "target", "interaction_name", "pathway_name",
  "ligand", "receptor", "ligand_family", "receptor_family", "signal_family"
)

comparison_raw <- full_join(
  result_ck$table %>%
    select(all_of(join_keys), prob_ck = prob,
           ligand_pct_ck = ligand_pct, receptor_pct_ck = receptor_pct,
           ligand_mean_ck = ligand_mean, receptor_mean_ck = receptor_mean,
           within_support_ck = supported_within_sample),
  result_fcr$table %>%
    select(all_of(join_keys), prob_fcr = prob,
           ligand_pct_fcr = ligand_pct, receptor_pct_fcr = receptor_pct,
           ligand_mean_fcr = ligand_mean, receptor_mean_fcr = receptor_mean,
           within_support_fcr = supported_within_sample),
  by = join_keys
) %>%
  mutate(
    across(c(prob_ck, prob_fcr, starts_with("ligand_pct_"),
             starts_with("receptor_pct_"), starts_with("ligand_mean_"),
             starts_with("receptor_mean_")), ~replace_na(.x, 0)),
    across(starts_with("within_support_"), ~replace_na(.x, FALSE)),
    delta_prob_raw = prob_fcr - prob_ck,
    normalized_delta = (prob_fcr - prob_ck) /
      (prob_fcr + prob_ck + .Machine$double.eps),
    log2_ratio_when_both_positive = if_else(
      prob_ck > 0 & prob_fcr > 0,
      log2(prob_fcr / prob_ck),
      NA_real_
    )
  ) %>%
  arrange(desc(delta_prob_raw))

# Collapse A/B/D DCC homoeologs so they are not counted as independent pathways.
family_comparison <- comparison_raw %>%
  group_by(source, target, pathway_name, ligand_family,
           receptor_family, signal_family) %>%
  summarise(
    n_gene_pairs = n(),
    cellchat_prob_ck_family_mean = mean(prob_ck),
    cellchat_prob_fcr_family_mean = mean(prob_fcr),
    ligand_pct_ck_family_max = max(ligand_pct_ck),
    ligand_pct_fcr_family_max = max(ligand_pct_fcr),
    receptor_pct_ck_family_max = max(receptor_pct_ck),
    receptor_pct_fcr_family_max = max(receptor_pct_fcr),
    ligand_mean_ck_family_mean = mean(ligand_mean_ck),
    ligand_mean_fcr_family_mean = mean(ligand_mean_fcr),
    receptor_mean_ck_family_mean = mean(receptor_mean_ck),
    receptor_mean_fcr_family_mean = mean(receptor_mean_fcr),
    any_within_support_ck = any(within_support_ck),
    any_within_support_fcr = any(within_support_fcr),
    .groups = "drop"
  ) %>%
  mutate(
    expression_supported_any_condition =
      pmax(ligand_pct_ck_family_max, ligand_pct_fcr_family_max) >= min_expression_pct &
      pmax(receptor_pct_ck_family_max, receptor_pct_fcr_family_max) >= min_expression_pct,
    delta_cellchat_prob = cellchat_prob_fcr_family_mean - cellchat_prob_ck_family_mean,
    normalized_delta_cellchat =
      (cellchat_prob_fcr_family_mean - cellchat_prob_ck_family_mean) /
      (cellchat_prob_fcr_family_mean + cellchat_prob_ck_family_mean +
         .Machine$double.eps),
    log2_ratio_when_both_positive = if_else(
      cellchat_prob_ck_family_mean > 0 & cellchat_prob_fcr_family_mean > 0,
      log2(cellchat_prob_fcr_family_mean / cellchat_prob_ck_family_mean),
      NA_real_
    )
  ) %>%
  arrange(desc(delta_cellchat_prob))

############################################################
## Repeated stratified cell resampling
## This measures cell-sampling stability, NOT biological uncertainty.
############################################################

expr_all <- get_data_layer(seu_use, assay_use)[all_lr_genes, , drop = FALSE]
meta_use <- seu_use@meta.data[colnames(expr_all), , drop = FALSE]

cells_by_condition_group <- lapply(c(condition_ck, condition_fcr), function(cond) {
  setNames(lapply(cell_order, function(ct) {
    rownames(meta_use)[meta_use$.condition == cond & meta_use$.celltype == ct]
  }), cell_order)
})
names(cells_by_condition_group) <- c(condition_ck, condition_fcr)

n_equal <- setNames(vapply(cell_order, function(ct) {
  min(length(cells_by_condition_group[[condition_ck]][[ct]]),
      length(cells_by_condition_group[[condition_fcr]][[ct]]))
}, integer(1)), cell_order)

if (any(n_equal < min_cells_group)) {
  stop("At least one retained cell type has too few cells for stratified resampling.")
}

bootstrap_group_summary <- function(cond) {
  ans_mean <- matrix(0, nrow = length(cell_order), ncol = length(all_lr_genes),
                     dimnames = list(cell_order, all_lr_genes))
  ans_fraction <- ans_mean
  for (ct in cell_order) {
    pool <- cells_by_condition_group[[cond]][[ct]]
    sampled <- sample(pool, size = n_equal[[ct]], replace = TRUE)
    sampled_expression <- expr_all[, sampled, drop = FALSE]
    ans_mean[ct, ] <- Matrix::rowMeans(sampled_expression)
    ans_fraction[ct, ] <- Matrix::rowMeans(sampled_expression > 0)
  }
  list(mean = ans_mean, fraction = ans_fraction)
}

route_design <- expand.grid(
  source = cell_order,
  target = cell_order,
  interaction_name = lr_db$interaction_name,
  stringsAsFactors = FALSE
) %>%
  left_join(
    lr_db %>% select(interaction_name, pathway_name, ligand, receptor,
                     ligand_family, receptor_family, signal_family),
    by = "interaction_name"
  ) %>%
  mutate(
    family_route = paste(source, target, signal_family, sep = "||")
  )

family_levels <- unique(route_design$family_route)
family_index <- match(route_design$family_route, family_levels)
family_size <- tabulate(family_index, nbins = length(family_levels))

if (anyNA(route_design[, c("source", "target", "ligand", "receptor",
                           "pathway_name", "signal_family", "family_route")])) {
  stop("route_design contains missing annotations; check interaction_name uniqueness.")
}
if (anyNA(family_index) || any(family_size < 1L)) {
  stop("Family-route indexing failed.")
}

family_meta <- route_design %>%
  group_by(family_route) %>%
  slice(1L) %>%
  ungroup() %>%
  select(family_route, source, target, pathway_name,
         ligand_family, receptor_family, signal_family) %>%
  slice(match(family_levels, family_route))

if (!identical(family_meta$family_route, family_levels)) {
  stop("Family-route metadata and averaging matrix are not aligned.")
}

# Sparse averaging operator: each row is one family route and each column is
# one gene-level LR route. A/B/D homoeolog scores are averaged, never summed.
family_averager <- Matrix::sparseMatrix(
  i = family_index,
  j = seq_along(family_index),
  x = 1 / family_size[family_index],
  dims = c(length(family_levels), nrow(route_design))
)

src_i <- match(route_design$source, cell_order)
tgt_i <- match(route_design$target, cell_order)
lig_i <- match(route_design$ligand, all_lr_genes)
rec_i <- match(route_design$receptor, all_lr_genes)
if (anyNA(c(src_i, tgt_i, lig_i, rec_i))) {
  stop("A route-design index is missing; the resampling analysis cannot proceed.")
}

boot_ck <- matrix(NA_real_, nrow = length(family_levels), ncol = nboot_stability)
boot_fcr <- boot_ck
boot_expression_support_ck <- matrix(
  FALSE, nrow = length(family_levels), ncol = nboot_stability,
  dimnames = dimnames(boot_ck)
)
boot_expression_support_fcr <- boot_expression_support_ck

set.seed(master_seed)
for (b in seq_len(nboot_stability)) {
  summary_ck <- bootstrap_group_summary(condition_ck)
  summary_fcr <- bootstrap_group_summary(condition_fcr)

  gene_score_ck <-
    summary_ck$mean[cbind(src_i, lig_i)] *
    summary_ck$mean[cbind(tgt_i, rec_i)]
  gene_score_fcr <-
    summary_fcr$mean[cbind(src_i, lig_i)] *
    summary_fcr$mean[cbind(tgt_i, rec_i)]

  gene_support_ck <-
    summary_ck$fraction[cbind(src_i, lig_i)] >= min_expression_pct &
    summary_ck$fraction[cbind(tgt_i, rec_i)] >= min_expression_pct
  gene_support_fcr <-
    summary_fcr$fraction[cbind(src_i, lig_i)] >= min_expression_pct &
    summary_fcr$fraction[cbind(tgt_i, rec_i)] >= min_expression_pct

  if (any(!is.finite(gene_score_ck)) || any(!is.finite(gene_score_fcr))) {
    stop("Non-finite expression product generated in resample ", b, ".")
  }

  boot_ck[, b] <- as.vector(family_averager %*% gene_score_ck)
  boot_fcr[, b] <- as.vector(family_averager %*% gene_score_fcr)
  boot_expression_support_ck[, b] <-
    as.vector(family_averager %*% as.numeric(gene_support_ck)) > 0
  boot_expression_support_fcr[, b] <-
    as.vector(family_averager %*% as.numeric(gene_support_fcr)) > 0
}

if (any(!is.finite(boot_ck)) || any(!is.finite(boot_fcr))) {
  stop("The completed cell-resampling matrices contain non-finite values.")
}

boot_denominator <- boot_fcr + boot_ck
boot_normalized_delta <- matrix(
  0,
  nrow = nrow(boot_ck), ncol = ncol(boot_ck),
  dimnames = dimnames(boot_ck)
)
positive_denominator <- boot_denominator > 0
boot_normalized_delta[positive_denominator] <-
  (boot_fcr[positive_denominator] - boot_ck[positive_denominator]) /
  boot_denominator[positive_denominator]

qrow <- function(x, p) apply(x, 1, quantile, probs = p, names = FALSE)

fraction_fcr_greater <- rowMeans(boot_fcr > boot_ck)
fraction_ck_greater <- rowMeans(boot_ck > boot_fcr)
fraction_equal <- rowMeans(boot_ck == boot_fcr)
expression_support_frequency_ck <- rowMeans(boot_expression_support_ck)
expression_support_frequency_fcr <- rowMeans(boot_expression_support_fcr)

stability_table <- family_meta %>%
  mutate(
    expression_product_ck_median = apply(boot_ck, 1, median),
    expression_product_ck_q025 = qrow(boot_ck, 0.025),
    expression_product_ck_q975 = qrow(boot_ck, 0.975),
    expression_product_fcr_median = apply(boot_fcr, 1, median),
    expression_product_fcr_q025 = qrow(boot_fcr, 0.025),
    expression_product_fcr_q975 = qrow(boot_fcr, 0.975),
    normalized_delta_median = apply(boot_normalized_delta, 1, median),
    normalized_delta_q025 = qrow(boot_normalized_delta, 0.025),
    normalized_delta_q975 = qrow(boot_normalized_delta, 0.975),
    fraction_fcr_greater = fraction_fcr_greater,
    fraction_ck_greater = fraction_ck_greater,
    fraction_equal = fraction_equal,
    direction_stability = pmax(fraction_fcr_greater, fraction_ck_greater),
    expression_support_frequency_ck = expression_support_frequency_ck,
    expression_support_frequency_fcr = expression_support_frequency_fcr,
    direction = case_when(
      fraction_equal == 1 ~ "no_observed_difference",
      fraction_fcr_greater >= fraction_ck_greater ~ "FCR_higher",
      TRUE ~ "CK_higher"
    )
  ) %>%
  left_join(
    family_comparison %>%
      select(source, target, pathway_name, ligand_family, receptor_family,
             signal_family, expression_supported_any_condition,
             cellchat_prob_ck_family_mean, cellchat_prob_fcr_family_mean,
             normalized_delta_cellchat, any_within_support_ck,
             any_within_support_fcr),
    by = c("source", "target", "pathway_name", "ligand_family",
           "receptor_family", "signal_family")
  ) %>%
  mutate(
    directional_expression_support_stability = if_else(
      direction == "FCR_higher",
      expression_support_frequency_fcr,
      expression_support_frequency_ck
    ),
    robust_descriptive_candidate = expression_supported_any_condition &
      pmax(expression_product_ck_median,
           expression_product_fcr_median) > 0 &
      direction_stability >= stability_cutoff &
      directional_expression_support_stability >= stability_cutoff &
      abs(normalized_delta_median) >= effect_cutoff
  ) %>%
  arrange(desc(direction_stability), desc(abs(normalized_delta_median)))

# Equal-weight pathway summaries. These are means across family-level routes,
# not sums, so pathways with more homoeolog pairs are not automatically larger.
pathway_descriptive <- family_comparison %>%
  group_by(pathway_name) %>%
  summarise(
    n_family_routes = n(),
    mean_family_prob_ck = mean(cellchat_prob_ck_family_mean),
    mean_family_prob_fcr = mean(cellchat_prob_fcr_family_mean),
    median_family_prob_ck = median(cellchat_prob_ck_family_mean),
    median_family_prob_fcr = median(cellchat_prob_fcr_family_mean),
    fraction_routes_fcr_higher = mean(cellchat_prob_fcr_family_mean >
                                        cellchat_prob_ck_family_mean),
    .groups = "drop"
  ) %>%
  mutate(
    normalized_delta_mean_probability =
      (mean_family_prob_fcr - mean_family_prob_ck) /
      (mean_family_prob_fcr + mean_family_prob_ck + .Machine$double.eps)
  ) %>%
  arrange(desc(normalized_delta_mean_probability))

############################################################
## Hard validation gate: export is forbidden if resampling failed
############################################################

stability_numeric_columns <- c(
  "expression_product_ck_median",
  "expression_product_ck_q025",
  "expression_product_ck_q975",
  "expression_product_fcr_median",
  "expression_product_fcr_q025",
  "expression_product_fcr_q975",
  "normalized_delta_median",
  "normalized_delta_q025",
  "normalized_delta_q975",
  "fraction_fcr_greater",
  "fraction_ck_greater",
  "fraction_equal",
  "direction_stability",
  "expression_support_frequency_ck",
  "expression_support_frequency_fcr",
  "directional_expression_support_stability"
)

if (!all(stability_numeric_columns %in% names(stability_table))) {
  stop("The stability table is missing required numeric columns.")
}
if (nrow(stability_table) != length(family_levels)) {
  stop("The stability table does not contain exactly one row per family route.")
}
finite_by_column <- vapply(
  stability_table[stability_numeric_columns],
  function(x) all(is.finite(x)),
  logical(1)
)
if (!all(finite_by_column)) {
  stop(
    "Non-finite stability output detected in: ",
    paste(names(finite_by_column)[!finite_by_column], collapse = ", ")
  )
}
if (any(stability_table$normalized_delta_q025 >
        stability_table$normalized_delta_median) ||
    any(stability_table$normalized_delta_median >
        stability_table$normalized_delta_q975)) {
  stop("At least one resampling interval does not contain its median.")
}
if (any(stability_table$direction_stability < 0 |
        stability_table$direction_stability > 1)) {
  stop("direction_stability is outside [0, 1].")
}
if (any(stability_table$directional_expression_support_stability < 0 |
        stability_table$directional_expression_support_stability > 1)) {
  stop("directional_expression_support_stability is outside [0, 1].")
}

validation_summary <- data.frame(
  check = c(
    "family_routes",
    "finite_stability_rows",
    "expression_supported_routes",
    "robust_descriptive_routes",
    "DCC1_expression_supported_routes",
    "DCC1_robust_descriptive_routes"
  ),
  value = c(
    nrow(stability_table),
    sum(complete.cases(stability_table[stability_numeric_columns])),
    sum(stability_table$expression_supported_any_condition),
    sum(stability_table$robust_descriptive_candidate),
    sum(stability_table$ligand_family == "DCC1_family" &
        stability_table$expression_supported_any_condition),
    sum(stability_table$ligand_family == "DCC1_family" &
        stability_table$robust_descriptive_candidate)
  ),
  stringsAsFactors = FALSE
)

############################################################
## Export: terminology deliberately avoids differential significance
############################################################

analysis_notes <- data.frame(
  item = c(
    "Independent biological units",
    "Allowed inference",
    "CellChat p/q meaning",
    "CK-FCR comparison",
    "Stability range meaning",
    "Stability score",
    "Homoeolog handling",
    "Random seed"
  ),
  value = c(
    "One CK sample and one FCR sample; cells are subsamples, not biological replicates.",
    "Descriptive association within the two observed samples only.",
    "Within-sample cell-label permutation screening; not a CK-versus-FCR treatment test.",
    "Raw continuous scores and effect sizes; no differential p value is reported.",
    "2.5%-97.5% cell-resampling range; not a biological confidence interval.",
    "A descriptive robustness flag requires full-data expression support, direction stability >= 0.90, directional expression-threshold support in >= 0.90 of cell resamples, and |normalized difference| >= 0.20; it is not a significance test.",
    "DCC1 A/B/D genes are averaged as one DCC1 family; DCC2 A/B/D genes are averaged separately as one DCC2 family. DCC1 and DCC2 are never added together.",
    "54 is the master seed; repeated resamples differ reproducibly within that RNG stream."
  ),
  stringsAsFactors = FALSE
)

analysis_parameters <- data.frame(
  parameter = c(
    "assay", "normalization", "normalization_scale_factor", "master_seed",
    "minimum_cells_per_celltype_per_condition",
    "minimum_expression_fraction", "truncated_mean_trim",
    "CellChat_permutations", "cell_resampling_iterations",
    "direction_stability_cutoff", "normalized_effect_cutoff"
  ),
  value = as.character(c(
    assay_use, "LogNormalize", 10000, master_seed, min_cells_group,
    min_expression_pct, trim_use,
    nboot_cellchat, nboot_stability, stability_cutoff, effect_cutoff
  )),
  stringsAsFactors = FALSE
)

software_versions <- data.frame(
  software = c("R", "Seurat", "CellChat", "Matrix", "dplyr", "ggplot2", "openxlsx"),
  version = c(
    R.version.string,
    as.character(packageVersion("Seurat")),
    as.character(packageVersion("CellChat")),
    as.character(packageVersion("Matrix")),
    as.character(packageVersion("dplyr")),
    as.character(packageVersion("ggplot2")),
    as.character(packageVersion("openxlsx"))
  ),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()
addWorksheet(wb, "analysis_notes")
writeData(wb, "analysis_notes", analysis_notes)
addWorksheet(wb, "parameters")
writeData(wb, "parameters", analysis_parameters)
addWorksheet(wb, "software_versions")
writeData(wb, "software_versions", software_versions)
addWorksheet(wb, "validation_summary")
writeData(wb, "validation_summary", validation_summary)
addWorksheet(wb, "cell_counts")
writeData(wb, "cell_counts", as.data.frame.matrix(table(meta_use$.celltype, meta_use$.condition)), rowNames = TRUE)
addWorksheet(wb, "input_LR_database")
writeData(wb, "input_LR_database", lr_db)
addWorksheet(wb, "CK_within_sample")
writeData(wb, "CK_within_sample", result_ck$table)
addWorksheet(wb, "FCR_within_sample")
writeData(wb, "FCR_within_sample", result_fcr$table)
addWorksheet(wb, "raw_descriptive_comparison")
writeData(wb, "raw_descriptive_comparison", comparison_raw)
addWorksheet(wb, "family_descriptive")
writeData(wb, "family_descriptive", family_comparison)
addWorksheet(wb, "pathway_descriptive")
writeData(wb, "pathway_descriptive", pathway_descriptive)
addWorksheet(wb, "cell_sampling_stability")
writeData(wb, "cell_sampling_stability", stability_table)
addWorksheet(wb, "DCC_stability")
writeData(wb, "DCC_stability", stability_table %>% filter(grepl("DCC", signal_family)))

saveWorkbook(
  wb,
  file.path(outdir, "wheat_CellChat_final_validated.xlsx"),
  overwrite = TRUE
)

############################################################
## Publication-oriented descriptive figures
## All CK-FCR displays are descriptive because n=1 per condition.
############################################################

theme_comm <- theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

circle_summary <- family_comparison %>%
  group_by(source, target) %>%
  summarise(
    retained_count_ck = sum(any_within_support_ck),
    retained_count_fcr = sum(any_within_support_fcr),
    # Equal-weight mean across all family hypotheses, with unsupported routes
    # contributing zero. This preserves a common denominator while preventing
    # visually dense edges driven only by routes that failed support screening.
    mean_family_strength_ck = mean(if_else(
      any_within_support_ck, cellchat_prob_ck_family_mean, 0
    )),
    mean_family_strength_fcr = mean(if_else(
      any_within_support_fcr, cellchat_prob_fcr_family_mean, 0
    )),
    .groups = "drop"
  )

cell_short <- c(
  "Defense-detoxification cells" = "Def-detox",
  "Defense-thickening cells I" = "Def-thick I",
  "Defense-thickening cells II" = "Def-thick II",
  "Epidermis cells I" = "Epidermis I",
  "Epidermis cells II" = "Epidermis II",
  "Guard cells" = "Guard",
  "JA-mediated defense cells" = "JA-defense",
  "Mesophyll cells I" = "Mesophyll I",
  "Mesophyll cells II" = "Mesophyll II",
  "Phloem cells" = "Phloem",
  "Proliferating-S cells" = "Prolif-S",
  "Protein-synthesis cells" = "Protein-synth",
  "Sclerenchyma cells" = "Sclerenchyma",
  "Trichome cells" = "Trichome",
  "Wounding-responsive cells" = "Wounding"
)

short_cell_name <- function(x) {
  out <- unname(cell_short[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

route_signal_name <- function(pathway, ligand_family) {
  ifelse(
    ligand_family == "DCC1_family", "DCC1-DCCR",
    ifelse(ligand_family == "DCC2_family", "DCC2-DCCR", pathway)
  )
}

to_sr_matrix <- function(data, value_col) {
  value_name <- rlang::as_name(rlang::ensym(value_col))
  mat <- matrix(
    0, nrow = length(cell_order), ncol = length(cell_order),
    dimnames = list(cell_order, cell_order)
  )
  ri <- match(data$source, cell_order)
  ci <- match(data$target, cell_order)
  mat[cbind(ri, ci)] <- data[[value_name]]
  mat
}

count_ck_mat <- to_sr_matrix(circle_summary, retained_count_ck)
count_fcr_mat <- to_sr_matrix(circle_summary, retained_count_fcr)
strength_ck_mat <- to_sr_matrix(circle_summary, mean_family_strength_ck)
strength_fcr_mat <- to_sr_matrix(circle_summary, mean_family_strength_fcr)

group_size_ck <- as.numeric(table(factor(
  meta_use$.celltype[meta_use$.condition == condition_ck], levels = cell_order
)))
group_size_fcr <- as.numeric(table(factor(
  meta_use$.celltype[meta_use$.condition == condition_fcr], levels = cell_order
)))

plot_circle_file <- function(mat, vertex_weight, filename, title_text,
                             common_edge_max) {
  pdf(file.path(outdir, filename), width = 12.5, height = 12.5,
      useDingbats = FALSE)
  if (any(mat > 0)) {
    # Base-R circle implementation avoids the CellChat/igraph `is.R` failure
    # seen with older CellChat and current igraph releases.
    n_node <- nrow(mat)
    theta <- seq(pi / 2, pi / 2 - 2 * pi, length.out = n_node + 1)[-1]
    node_x <- cos(theta)
    node_y <- sin(theta)
    node_col <- grDevices::hcl.colors(n_node, palette = "Dark 3")
    node_scale <- sqrt(vertex_weight / max(vertex_weight))

    old_par <- par(mar = c(1.5, 1.5, 3.2, 1.5), xpd = NA)
    plot.new()
    plot.window(xlim = c(-1.48, 1.48), ylim = c(-1.48, 1.48), asp = 1)

    edge_index <- which(mat > 0, arr.ind = TRUE)
    for (k in seq_len(nrow(edge_index))) {
      i <- edge_index[k, 1]
      j <- edge_index[k, 2]
      w <- mat[i, j]
      edge_lwd <- 0.5 + 6 * sqrt(w / common_edge_max)
      edge_col <- grDevices::adjustcolor(node_col[i], alpha.f = 0.38)

      if (i != j) {
        arrows(
          node_x[i], node_y[i], node_x[j], node_y[j],
          length = 0.055, angle = 22,
          lwd = edge_lwd, col = edge_col
        )
      } else {
        symbols(
          node_x[i] * 1.055, node_y[i] * 1.055,
          circles = 0.105,
          inches = FALSE, add = TRUE,
          fg = edge_col, lwd = edge_lwd
        )
      }
    }

    points(
      node_x, node_y,
      pch = 21, bg = node_col, col = "white", lwd = 1.2,
      cex = 2.5 + 4.5 * node_scale
    )
    for (i in seq_len(n_node)) {
      text(
        1.28 * node_x[i], 1.28 * node_y[i],
        labels = short_cell_name(cell_order[i]),
        cex = 0.78,
        srt = if (node_x[i] < 0) theta[i] * 180 / pi + 180 else
          theta[i] * 180 / pi,
        adj = if (node_x[i] < 0) 1 else 0
      )
    }
    title(main = title_text, font.main = 2, cex.main = 1.15)
    par(old_par)
  } else {
    plot.new()
    text(0.5, 0.5, paste0(title_text, "\nNo retained routes"), cex = 1.1)
  }
  dev.off()
}

count_edge_max <- max(c(count_ck_mat, count_fcr_mat), na.rm = TRUE)
strength_edge_max <- max(c(strength_ck_mat, strength_fcr_mat), na.rm = TRUE)
if (!is.finite(count_edge_max) || count_edge_max <= 0) count_edge_max <- 1
if (!is.finite(strength_edge_max) || strength_edge_max <= 0) strength_edge_max <- 1

plot_circle_file(
  count_ck_mat, group_size_ck, "CK_retained_LR_count_circle.pdf",
  "CK: retained family-level LR route count", count_edge_max
)
plot_circle_file(
  count_fcr_mat, group_size_fcr, "FCR_retained_LR_count_circle.pdf",
  "FCR: retained family-level LR route count", count_edge_max
)
plot_circle_file(
  strength_ck_mat, group_size_ck, "CK_retained_strength_circle.pdf",
  "CK: retained family-level communication strength", strength_edge_max
)
plot_circle_file(
  strength_fcr_mat, group_size_fcr, "FCR_retained_strength_circle.pdf",
  "FCR: retained family-level communication strength", strength_edge_max
)

# CK and FCR heatmaps on one common scale. The route-family mean prevents DCC
# A/B/D homoeologs and pathways with more mapped pairs from dominating by sum.
sr_strength_long <- bind_rows(
  circle_summary %>% transmute(source, target, condition = "CK",
                               mean_family_probability = mean_family_strength_ck),
  circle_summary %>% transmute(source, target, condition = "FCR",
                               mean_family_probability = mean_family_strength_fcr)
) %>%
  mutate(
    source = factor(source, levels = rev(cell_order)),
    target = factor(target, levels = cell_order)
  )

pdf(file.path(outdir, "CK_FCR_family_strength_heatmaps.pdf"),
    width = 15, height = 7, useDingbats = FALSE)
print(
  ggplot(sr_strength_long, aes(target, source, fill = mean_family_probability)) +
    geom_tile(color = "grey92", linewidth = 0.25) +
    facet_wrap(~condition, nrow = 1) +
    scale_fill_gradient(low = "white", high = "#B2182B",
                        name = "Mean family\nprobability",
                        labels = scales::label_scientific(digits = 1)) +
    labs(
      x = "Receiver", y = "Sender",
      title = "Family-level communication strength",
      subtitle = "Common color scale"
    ) +
    theme_comm + theme(legend.position = "right")
)
dev.off()

sr_delta <- circle_summary %>%
  mutate(
    normalized_delta = (mean_family_strength_fcr - mean_family_strength_ck) /
      (mean_family_strength_fcr + mean_family_strength_ck + .Machine$double.eps),
    source = factor(source, levels = rev(cell_order)),
    target = factor(target, levels = cell_order)
  )

pdf(file.path(outdir, "CK_FCR_descriptive_delta_heatmap.pdf"),
    width = 10, height = 9, useDingbats = FALSE)
print(
  ggplot(sr_delta, aes(target, source, fill = normalized_delta)) +
    geom_tile(color = "grey92", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1),
      name = "Normalized\ndifference"
    ) +
    labs(
      x = "Receiver", y = "Sender",
      title = "Descriptive change in mean family-level communication",
      subtitle = "Normalized contrast: FCR - CK"
    ) +
    theme_comm
)
dev.off()

# Always display the strongest expression-supported routes. The robustness flag
# is encoded instead of deleting all rows when no route crosses the cutoff.
top_plot <- stability_table %>%
  filter(robust_descriptive_candidate,
         is.finite(normalized_delta_median),
         is.finite(normalized_delta_q025),
         is.finite(normalized_delta_q975)) %>%
  mutate(rank_score = abs(normalized_delta_median) * direction_stability) %>%
  slice_max(order_by = rank_score, n = 30, with_ties = FALSE) %>%
  mutate(
    route_label = paste(
      short_cell_name(source), "->", short_cell_name(target),
      route_signal_name(pathway_name, ligand_family), sep = " | "
    ),
    route_label = reorder(route_label, normalized_delta_median),
    direction_label = if_else(
      normalized_delta_median >= 0, "FCR higher", "CK higher"
    )
  )

pdf(file.path(outdir, "top_stable_descriptive_routes.pdf"),
    width = 12, height = 10, useDingbats = FALSE)
if (nrow(top_plot) > 0L) {
  print(
    ggplot(top_plot, aes(normalized_delta_median, route_label,
                         color = direction_label)) +
      geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
      geom_segment(
        aes(x = normalized_delta_q025, xend = normalized_delta_q975,
            y = route_label, yend = route_label),
        linewidth = 0.55
      ) +
      geom_point(size = 2.3) +
      scale_color_manual(values = c(
        "CK higher" = "#4477AA", "FCR higher" = "#CC6677"
      ), name = NULL) +
      labs(
        x = "Normalized descriptive difference (FCR - CK)",
        y = NULL,
        title = "Most stable family-level LR contrasts",
        subtitle = "Median and 2.5%-97.5% cell-resampling range"
      ) +
      theme_comm + theme(
        axis.text.x = element_text(angle = 0),
        legend.position = "bottom",
        plot.margin = margin(8, 12, 8, 8)
      )
  )
} else {
  plot.new()
  text(0.5, 0.5, "No expression-supported routes were available", cex = 1.1)
}
dev.off()

# Direct CK/FCR route display based on CellChat probabilities, independent of
# whether a route meets the cell-resampling robustness cutoff.
top_retained <- family_comparison %>%
  filter(expression_supported_any_condition) %>%
  mutate(max_probability = pmax(cellchat_prob_ck_family_mean,
                                cellchat_prob_fcr_family_mean)) %>%
  slice_max(max_probability, n = 30, with_ties = FALSE) %>%
  mutate(route_label = paste(
    short_cell_name(source), "->", short_cell_name(target),
    route_signal_name(pathway_name, ligand_family), sep = " | "
  )) %>%
  select(route_label, pathway_name,
         CK = cellchat_prob_ck_family_mean,
         FCR = cellchat_prob_fcr_family_mean) %>%
  pivot_longer(c(CK, FCR), names_to = "condition", values_to = "probability") %>%
  mutate(route_label = reorder(route_label, probability, FUN = max))

pdf(file.path(outdir, "top_retained_LR_routes.pdf"),
    width = 12, height = 10, useDingbats = FALSE)
print(
  ggplot(top_retained, aes(probability, route_label, color = condition)) +
    geom_line(aes(group = route_label), color = "grey75", linewidth = 0.45) +
    geom_point(size = 2.3) +
    scale_color_manual(values = c(CK = "#4477AA", FCR = "#CC6677")) +
    labs(
      x = "Family-level CellChat probability", y = NULL,
      title = "Top expression-supported LR routes",
      subtitle = "Expression-supported family-level routes"
    ) +
    theme_comm + theme(axis.text.x = element_text(angle = 0))
)
dev.off()

dcc_heatmap <- stability_table %>%
  filter(ligand_family == "DCC1_family",
         robust_descriptive_candidate,
         is.finite(normalized_delta_median)) %>%
  mutate(
    source = factor(source, levels = rev(cell_order)),
    target = factor(target, levels = cell_order)
  )

dcc_background <- expand.grid(
  source = factor(rev(cell_order), levels = rev(cell_order)),
  target = factor(cell_order, levels = cell_order),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

pdf(file.path(outdir, "DCC_family_stability_heatmap.pdf"),
    width = 15, height = 8, useDingbats = FALSE)
if (nrow(dcc_heatmap) > 0L) {
  print(
    ggplot() +
      geom_tile(
        data = dcc_background,
        aes(target, source),
        fill = "grey96", color = "white", linewidth = 0.25
      ) +
      geom_tile(
        data = dcc_heatmap,
        aes(target, source, fill = normalized_delta_median),
        color = "white", linewidth = 0.25
      ) +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, limits = c(-1, 1),
        name = "Normalized\ndifference"
      ) +
      labs(
        x = "Receiver", y = "Sender",
        title = "Robust DCC1-DCCR cell-sampling pattern",
        subtitle = "Colored routes pass direction, effect-size and expression-support stability criteria"
      ) +
      theme_comm
  )
} else {
  plot.new()
  text(0.5, 0.5, "No expression-supported DCC-family routes", cex = 1.1)
}
dev.off()

dcc_genes <- c(dcc1_genes, dcc2_genes)
dcc_expression <- bind_rows(
  result_ck$table %>%
    transmute(condition = "CK", celltype = source, gene = ligand,
              expression_fraction = ligand_pct,
              mean_log_normalized_expression = ligand_mean),
  result_fcr$table %>%
    transmute(condition = "FCR", celltype = source, gene = ligand,
              expression_fraction = ligand_pct,
              mean_log_normalized_expression = ligand_mean)
) %>%
  filter(gene %in% dcc_genes) %>%
  distinct() %>%
  mutate(
    family = if_else(gene %in% dcc1_genes, "DCC1", "DCC2"),
    celltype = factor(celltype, levels = rev(cell_order))
  )

pdf(file.path(outdir, "DCC_homoeolog_expression_dotplot.pdf"),
    width = 13, height = 9, useDingbats = FALSE)
print(
  ggplot(dcc_expression,
         aes(gene, celltype, size = expression_fraction,
             color = mean_log_normalized_expression)) +
    geom_point() +
    facet_grid(family ~ condition, scales = "free_x", space = "free_x") +
    scale_size(range = c(0, 6), limits = c(0, 1),
               name = "Fraction\nexpressing") +
    scale_color_viridis_c(option = "C", name = "Mean log-normalized\nexpression") +
    labs(
      x = "DCC homoeolog", y = "Cell type",
      title = "DCC1 and DCC2 homoeolog expression",
      subtitle = "Gene-level expression shown separately; no homoeolog summation"
    ) +
    theme_comm
)
dev.off()

pathway_plot <- pathway_descriptive %>%
  select(pathway_name, CK = mean_family_prob_ck,
         FCR = mean_family_prob_fcr) %>%
  pivot_longer(c(CK, FCR), names_to = "condition", values_to = "mean_probability")

pdf(file.path(outdir, "pathway_strength_CK_FCR.pdf"),
    width = 9, height = 6, useDingbats = FALSE)
print(
  ggplot(pathway_plot,
         aes(reorder(pathway_name, mean_probability, FUN = max),
             mean_probability, color = condition, group = pathway_name)) +
    geom_line(color = "grey75", linewidth = 0.5) +
    geom_point(size = 2.8) +
    coord_flip() +
    scale_color_manual(values = c(CK = "#4477AA", FCR = "#CC6677")) +
    labs(
      x = NULL, y = "Mean family-level probability",
      title = "Pathway-level communication summary",
      subtitle = "Equal-weight mean across family-level routes"
    ) +
    theme_comm + theme(axis.text.x = element_text(angle = 0))
)
dev.off()

saveRDS(result_ck$cellchat, file.path(outdir, "cellchat_ck_within_sample.rds"))
saveRDS(result_fcr$cellchat, file.path(outdir, "cellchat_fcr_within_sample.rds"))

message("Completed. Output directory: ", normalizePath(outdir))
message("Important: no CK-versus-FCR inferential p value is reported because n=1 per condition.")
message("Generated 11 PDF figures, one XLSX workbook, and two CellChat RDS objects.")
