suppressPackageStartupMessages({
  library(monocle)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 06_audit_branch_points.R SELECTED_CDS_RDS OUTPUT_DIRECTORY")
}

set.seed(54)
input_file <- normalizePath(args[[1]], mustWork = TRUE)
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cds <- readRDS(input_file)
target_gene <- "TraesCS4A02G001300"
branch_points <- cds@auxOrderingData[["DDRTree"]]$branch_points
if (length(branch_points) < 1L) stop("No DDRTree branch point was found")

audit_one_branch <- function(bp) {
  branch_cds <- buildBranchCellDataSet(
    cds,
    branch_point = bp,
    progenitor_method = "duplicate",
    branch_labels = c("Branch_A", "Branch_B")
  )
  meta <- as.data.frame(pData(branch_cds))
  meta$cell_id <- rownames(meta)
  meta$Samples <- tolower(as.character(meta$Samples))
  meta$Branch <- as.character(meta$Branch)
  meta$State <- as.character(meta$State)

  # Duplicated progenitor cells have duplicated row names; preserve both records
  # for branch-model bookkeeping, but also report unique biological cells.
  meta$original_cell_id <- sub("\\.1$", "", meta$cell_id)
  meta$is_duplicated_record <- duplicated(meta$original_cell_id) |
    duplicated(meta$original_cell_id, fromLast = TRUE)

  expr_target <- as.numeric(exprs(branch_cds)[target_gene, ]) /
    as.numeric(sizeFactors(branch_cds))
  meta$target_log1p_norm <- log1p(expr_target)

  branch_summary <- meta %>%
    group_by(Branch) %>%
    summarise(
      branch_records = n(),
      unique_cells = n_distinct(original_cell_id),
      CK_records = sum(Samples == "ck"),
      FCR_records = sum(Samples == "fcr"),
      FCR_fraction = FCR_records / branch_records,
      median_pseudotime = median(Pseudotime),
      target_fraction_expressing = mean(target_log1p_norm > 0),
      target_mean_log1p_norm = mean(target_log1p_norm),
      .groups = "drop"
    ) %>%
    mutate(branch_point = bp, .before = 1)

  state_summary <- meta %>%
    count(Branch, State, Samples, name = "records") %>%
    group_by(Branch) %>%
    mutate(branch_records = sum(records), fraction_of_branch = records / branch_records) %>%
    ungroup() %>%
    mutate(branch_point = bp, .before = 1)

  sample_summary <- meta %>%
    count(Branch, Samples, name = "records") %>%
    group_by(Branch) %>%
    mutate(branch_records = sum(records), fraction_of_branch = records / branch_records) %>%
    ungroup() %>%
    mutate(branch_point = bp, .before = 1)

  duplicated_summary <- meta %>%
    summarise(
      branch_point = bp,
      branch_records = n(),
      unique_cells = n_distinct(original_cell_id),
      duplicated_progenitor_records = sum(is_duplicated_record),
      duplicated_progenitor_unique_cells = n_distinct(original_cell_id[is_duplicated_record])
    )

  write.csv(meta, file.path(outdir, sprintf("branch_point_%d_cell_records.csv", bp)), row.names = FALSE)
  saveRDS(branch_cds, file.path(outdir, sprintf("branch_point_%d_branch_cds.rds", bp)), compress = FALSE)

  list(
    branch = branch_summary,
    state = state_summary,
    sample = sample_summary,
    duplicate = duplicated_summary
  )
}

audits <- lapply(seq_along(branch_points), audit_one_branch)
write.csv(bind_rows(lapply(audits, `[[`, "branch")),
          file.path(outdir, "branch_point_summary.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(audits, `[[`, "state")),
          file.path(outdir, "branch_point_state_composition.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(audits, `[[`, "sample")),
          file.path(outdir, "branch_point_sample_composition.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(audits, `[[`, "duplicate")),
          file.path(outdir, "branch_point_duplicate_audit.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))

print(bind_rows(lapply(audits, `[[`, "branch")))
