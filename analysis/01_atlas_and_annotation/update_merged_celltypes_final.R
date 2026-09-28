suppressPackageStartupMessages(library(Seurat))

set.seed(54)
base_dir <- "/home/xiehangfei/scRNAseq"
src <- file.path(base_dir, "merged.rds")
tmp <- file.path(base_dir, "merged.celltype_revision.tmp.rds")
bak <- file.path(base_dir, "merged_before_celltype_revision_20260908.rds")
outdir <- file.path(base_dir, "NO1_atlas_20260908")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

old_to_final <- c(
  "Defense-detoxification cells" = "Defense detoxification cells",
  "Defense-thickening cells I" = "Defense cell wall remodeling cells I",
  "Defense-thickening cells II" = "Defense cell wall remodeling cells II",
  "Epidermis cells I" = "Epidermal cells I",
  "Epidermis cells II" = "Epidermal cells II",
  "Guard cells" = "Guard cells",
  "JA-mediated defense cells" = "Defense JA responsive cells",
  "Mesophyll cells I" = "Mesophyll cells I",
  "Mesophyll cells II" = "Mesophyll cells II",
  "Phloem cells" = "Procambial cells",
  "Proliferating-S cells" = "Proliferating S phase cells",
  "Protein-synthesis cells" = "Translation active cells",
  "Sclerenchyma cells" = "Vascular parenchyma cells",
  "Trichome cells" = "Phloem parenchyma cells",
  "Wounding-responsive cells" = "Epidermal cells III"
)

final_levels <- unname(old_to_final)
functional_levels <- c(
  "Defense detoxification cells",
  "Defense cell wall remodeling cells I",
  "Defense cell wall remodeling cells II",
  "Defense JA responsive cells",
  "Proliferating S phase cells",
  "Translation active cells"
)

obj <- readRDS(src)
stopifnot(inherits(obj, "Seurat"))
stopifnot("celltype" %in% colnames(obj@meta.data))
old <- as.character(obj$celltype)
unknown <- setdiff(unique(old), names(old_to_final))
if (length(unknown)) stop("Unmapped labels: ", paste(unknown, collapse = ", "))

before <- data.frame(
  old_celltype = names(table(old)),
  n_cells = as.integer(table(old)),
  stringsAsFactors = FALSE
)
before$final_celltype <- unname(old_to_final[before$old_celltype])

obj$celltype_original <- old
new_labels <- unname(old_to_final[old])
obj$celltype <- factor(new_labels, levels = final_levels)
obj$celltype_final <- factor(new_labels, levels = final_levels)
obj$annotation_basis <- ifelse(new_labels %in% functional_levels,
                               "transcriptional state", "anatomical identity")
Idents(obj) <- "celltype"

after <- as.data.frame(table(celltype = obj$celltype), stringsAsFactors = FALSE)
colnames(after)[2] <- "n_cells"
after <- after[after$n_cells > 0, ]
stopifnot(nrow(before) == 15L, nrow(after) == 15L)
stopifnot(sum(before$n_cells) == ncol(obj), sum(after$n_cells) == ncol(obj))
stopifnot(all(before$n_cells == after$n_cells[match(before$final_celltype, after$celltype)]))

write.table(before, file.path(outdir, "celltype_annotation_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(after, file.path(outdir, "celltype_counts_final.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (file.exists(tmp)) file.remove(tmp)
saveRDS(obj, tmp, compress = "gzip")
rm(obj)
gc()

check <- readRDS(tmp)
stopifnot(inherits(check, "Seurat"))
stopifnot(ncol(check) == sum(after$n_cells))
stopifnot(all(c("celltype_original", "celltype", "celltype_final", "annotation_basis") %in%
              colnames(check@meta.data)))
check_counts <- table(check$celltype)
stopifnot(identical(names(check_counts), final_levels))
stopifnot(all(as.integer(check_counts) == after$n_cells))
stopifnot(all(as.character(Idents(check)) == as.character(check$celltype)))

validation <- data.frame(
  check = c("Seurat object reload", "cell count preserved", "15 labels present",
            "per-label counts preserved", "active identity updated",
            "RNA assay retained", "UMAP retained"),
  passed = c(TRUE, ncol(check) == sum(after$n_cells), length(check_counts) == 15L,
             all(as.integer(check_counts) == after$n_cells),
             all(as.character(Idents(check)) == as.character(check$celltype)),
             "RNA" %in% names(check@assays), "umap" %in% names(check@reductions)),
  stringsAsFactors = FALSE
)
write.table(validation, file.path(outdir, "merged_rds_validation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(validation$passed)) stop("Validation failed; original file was not replaced.")
rm(check)
gc()

if (file.exists(bak)) stop("Server backup path already exists; refusing to overwrite it: ", bak)
if (!file.rename(src, bak)) stop("Could not rename original to backup.")
if (!file.rename(tmp, src)) {
  file.rename(bak, src)
  stop("Could not promote validated file; original restored.")
}

cat("Validated replacement completed.\n")
cat("Current object: ", src, "\n", sep = "")
cat("Server backup: ", bak, "\n", sep = "")
print(before, row.names = FALSE)
print(validation, row.names = FALSE)
