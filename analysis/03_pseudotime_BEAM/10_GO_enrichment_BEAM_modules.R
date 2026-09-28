suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Ta.eg.db)
  library(AnnotationDbi)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript 10_GO_enrichment_BEAM_modules.R MODULE_CSV BACKGROUND_GENE_CSV BRANCH_POINT OUTPUT_DIRECTORY")
}

set.seed(54)
module_file <- normalizePath(args[[1]], mustWork = TRUE)
background_file <- normalizePath(args[[2]], mustWork = TRUE)
branch_point <- as.integer(args[[3]])
outdir <- args[[4]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

modules <- read.csv(module_file, stringsAsFactors = FALSE)
background_df <- read.csv(background_file, stringsAsFactors = FALSE)
background <- unique(background_df$gene_id)

valid_keys <- keys(org.Ta.eg.db, keytype = "GID")
background_mapped <- intersect(background, valid_keys)

run_one <- function(module_name) {
  input <- unique(modules$gene_id[modules$module == module_name])
  mapped <- intersect(input, background_mapped)
  if (length(mapped) < 10L) return(NULL)

  ego <- enrichGO(
    gene = mapped,
    universe = background_mapped,
    OrgDb = org.Ta.eg.db,
    keyType = "GID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE
  )
  full <- as.data.frame(ego)
  if (nrow(full) > 0L) {
    full$module <- module_name
    full$branch_point <- branch_point
  }

  simplified <- tryCatch(
    as.data.frame(simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)),
    error = function(e) full
  )
  if (nrow(simplified) > 0L) {
    simplified$module <- module_name
    simplified$branch_point <- branch_point
  }

  mapping <- data.frame(
    branch_point = branch_point,
    module = module_name,
    input_genes = length(input),
    mapped_genes = length(mapped),
    mapping_fraction = length(mapped) / length(input),
    background_genes = length(background),
    background_mapped_genes = length(background_mapped),
    significant_BP_terms_BH_0.05 = sum(full$p.adjust < 0.05, na.rm = TRUE)
  )
  list(full = full, simplified = simplified, mapping = mapping)
}

res <- lapply(sort(unique(modules$module)), run_one)
res <- Filter(Negate(is.null), res)
full_all <- bind_rows(lapply(res, `[[`, "full"))
simplified_all <- bind_rows(lapply(res, `[[`, "simplified"))
mapping_all <- bind_rows(lapply(res, `[[`, "mapping"))

display <- simplified_all %>%
  filter(p.adjust < 0.05, Count >= 10) %>%
  group_by(module) %>%
  arrange(p.adjust, desc(Count), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(selection_rule = "Top 5 nonredundant BP terms per module by BH-adjusted P among terms supported by >=10 module genes; semantic similarity cutoff 0.70")

write.csv(full_all, file.path(outdir, sprintf("BEAM_branch_point_%d_GO_BP_full.csv", branch_point)), row.names = FALSE)
write.csv(simplified_all, file.path(outdir, sprintf("BEAM_branch_point_%d_GO_BP_simplified.csv", branch_point)), row.names = FALSE)
write.csv(display, file.path(outdir, sprintf("BEAM_branch_point_%d_GO_BP_display_terms.csv", branch_point)), row.names = FALSE)
write.csv(mapping_all, file.path(outdir, sprintf("BEAM_branch_point_%d_GO_mapping_QA.csv", branch_point)), row.names = FALSE)
write.csv(data.frame(gene_id = background),
          file.path(outdir, sprintf("BEAM_branch_point_%d_GO_background_all.csv", branch_point)), row.names = FALSE)
write.csv(data.frame(gene_id = background_mapped),
          file.path(outdir, sprintf("BEAM_branch_point_%d_GO_background_mapped.csv", branch_point)), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, sprintf("sessionInfo_GO_branch_point_%d.txt", branch_point)))

print(mapping_all)
print(display %>% select(module, ID, Description, GeneRatio, BgRatio, p.adjust, Count))
