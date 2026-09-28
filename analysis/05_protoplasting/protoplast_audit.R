set.seed(54)

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

options(stringsAsFactors = FALSE)

root <- "/home/xiehangfei/scRNAseq"
set_dir <- file.path(root, "protoplast_optimization_20260905")
outdir <- file.path(root, "protoplast_audit_suppfig_20260907")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "source_data"), showWarnings = FALSE)
dir.create(file.path(outdir, "source_data", "input_gene_sets"), showWarnings = FALSE)

obj_path <- file.path(root, "merged.rds")
marker_path <- file.path(root, "merged_all_markers.tsv")
deg_dir <- file.path(root, "All_Clusters_DEG_Results")
go_dir <- file.path(root, "All_Clusters_GO_Results")
scenic_path <- file.path(root, "SCENIC_Summary.csv")

stopifnot(file.exists(obj_path), file.exists(marker_path), file.exists(scenic_path))
seu <- readRDS(obj_path)
stopifnot(all(c("celltype", "Samples") %in% colnames(seu@meta.data)))
stopifnot("umap" %in% Reductions(seu))
seu$condition <- factor(tolower(as.character(seu$Samples)), levels = c("ck", "fcr"), labels = c("CK", "FCR"))

celltype_order <- names(sort(table(seu$celltype), decreasing = TRUE))
seu$celltype <- factor(as.character(seu$celltype), levels = celltype_order)

rice_map <- fread(file.path(set_dir, "high_confidence_rice_to_wheat_protoplast_RBH.tsv"))
ath_map <- fread(file.path(set_dir, "arabidopsis_to_wheat_protoplast_strict_RBH.tsv"))
rice_all <- unique(rice_map$wheat_gene)
rice_ultra <- unique(rice_map[ultra_core == 1, wheat_gene])
ath_all <- unique(ath_map$s)
conserved <- intersect(rice_all, ath_all)

rna_genes <- rownames(seu[["RNA"]])
sets_raw <- list(
  Conserved = conserved,
  Rice_direct = rice_all,
  Arabidopsis_direct = ath_all,
  Ultra_core = rice_ultra
)
sets <- lapply(sets_raw, intersect, y = rna_genes)
if (length(sets$Conserved) < 20) stop("Conserved signature has fewer than 20 detected genes.")

gene_set_summary <- data.frame(
  signature = names(sets),
  mapped_wheat_genes = lengths(sets_raw),
  detected_in_RNA = lengths(sets),
  used_for_UCell = c(TRUE, FALSE, FALSE, TRUE),
  role = c("Primary cross-study score", "Overlap audit", "Overlap audit", "Sensitivity score")
)
write.table(gene_set_summary, file.path(outdir, "source_data", "01_gene_set_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

gene_membership <- bind_rows(lapply(names(sets_raw), function(nm) {
  data.frame(signature = nm, wheat_gene = sets_raw[[nm]], detected_in_RNA = sets_raw[[nm]] %in% rna_genes)
}))
write.table(gene_membership, file.path(outdir, "source_data", "02_gene_set_membership.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

literature <- tribble(
  ~study, ~species, ~tissue, ~experimental_definition, ~reported_size, ~use_in_this_audit, ~doi, ~url,
  "Zhu et al. 2025", "Oryza sativa", "Root", "Bulk intact tissue versus protoplast; induced >=2-fold", "Source supplementary table", "Strict rice-to-wheat RBH; direct and conserved layers", "10.1038/s41586-025-08941-z", "https://doi.org/10.1038/s41586-025-08941-z",
  "Liew et al. 2024", "Arabidopsis thaliana", "Germinating embryo", "Bulk intact embryo versus protoplast; FDR <0.01 and |log2FC|>1.5", "782 induced; 420 reduced", "Strict Arabidopsis-to-wheat RBH; induced genes only", "10.1038/s41477-024-01771-3", "https://doi.org/10.1038/s41477-024-01771-3",
  "Frank et al. 2023", "Lotus japonicus", "Root", "Bulk whole root versus protoplast", "655 induced", "Contextual catalogue only; not merged into score", "10.1038/s41467-023-42911-1", "https://doi.org/10.1038/s41467-023-42911-1",
  "Xu et al. 2021", "Oryza sativa", "Root and shoot", "Published protoplasting-response exclusion set", "61-gene subset reused by later rice studies", "Contextual catalogue; covered by newer rice direct comparison", "10.1016/j.devcel.2020.12.015", "https://doi.org/10.1016/j.devcel.2020.12.015",
  "Zong et al. 2022", "Oryza sativa", "Inflorescence", "Sensitivity analysis after excluding 61 published protoplasting genes", "61", "Methodological precedent only", "10.1111/nph.18008", "https://doi.org/10.1111/nph.18008",
  "Shahan et al. 2022", "Arabidopsis thaliana", "Root and shoot", "Compared multiple published protoplasting signatures", "346 and larger tissue-specific lists", "Contextual catalogue only because signatures were tissue/context dependent", "10.1111/tpj.15659", "https://doi.org/10.1111/tpj.15659"
)
write.table(literature, file.path(outdir, "source_data", "03_literature_gene_set_catalogue.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Rank-based per-cell scores. Counts are used only for ranks; no expression values are overwritten.
seu <- AddModuleScore_UCell(
  seu,
  features = list(Protoplast_conserved = sets$Conserved, Protoplast_ultra = sets$Ultra_core),
  assay = "RNA", slot = "counts", maxRank = 3000, ncores = 8,
  name = "_UCell"
)
score_cols <- grep("^Protoplast_(conserved|ultra).*UCell$", colnames(seu@meta.data), value = TRUE)
if (length(score_cols) != 2) stop("Unexpected UCell score column names: ", paste(score_cols, collapse = ", "))
primary_col <- score_cols[grepl("conserved", score_cols)]
ultra_col <- score_cols[grepl("ultra", score_cols)]

um <- Embeddings(seu, "umap")
cell_scores <- data.frame(
  cell = colnames(seu),
  UMAP_1 = um[, 1], UMAP_2 = um[, 2],
  cluster = as.character(seu$celltype),
  condition = as.character(seu$condition),
  conserved_score = seu@meta.data[[primary_col]],
  ultra_score = seu@meta.data[[ultra_col]],
  stringsAsFactors = FALSE
)
cell_score_con <- gzfile(file.path(outdir, "source_data", "04_cell_level_scores.tsv.gz"), open = "wt")
write.table(cell_scores, cell_score_con, sep = "\t", row.names = FALSE, quote = FALSE)
close(cell_score_con)

score_summary <- cell_scores %>%
  group_by(cluster, condition) %>%
  summarise(n_cells = n(), median = median(conserved_score), q1 = quantile(conserved_score, .25),
            q3 = quantile(conserved_score, .75), mean = mean(conserved_score), sd = sd(conserved_score), .groups = "drop")
write.table(score_summary, file.path(outdir, "source_data", "05_score_by_cluster_condition.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Cell-stratified bootstrap is a stability analysis, not biological-replicate inference.
set.seed(54)
B <- 500L
clusters <- levels(seu$celltype)
boot_delta <- lapply(clusters, function(ct) {
  ck <- cell_scores$conserved_score[cell_scores$cluster == ct & cell_scores$condition == "CK"]
  fr <- cell_scores$conserved_score[cell_scores$cluster == ct & cell_scores$condition == "FCR"]
  obs <- median(fr) - median(ck)
  bs <- replicate(B, median(sample(fr, length(fr), replace = TRUE)) - median(sample(ck, length(ck), replace = TRUE)))
  data.frame(cluster = ct, n_CK = length(ck), n_FCR = length(fr), delta_median = obs,
             cell_resample_low = quantile(bs, .025), cell_resample_high = quantile(bs, .975),
             sign_stability = max(mean(bs > 0), mean(bs < 0)))
}) %>% bind_rows()
write.table(boot_delta, file.path(outdir, "source_data", "06_cell_stratified_bootstrap_deltas.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

sens_delta <- cell_scores %>%
  group_by(cluster, condition) %>%
  summarise(Conserved = median(conserved_score), Ultra_core = median(ultra_score), .groups = "drop") %>%
  pivot_longer(c(Conserved, Ultra_core), names_to = "signature", values_to = "median") %>%
  pivot_wider(names_from = condition, values_from = median) %>%
  mutate(delta = FCR - CK)
write.table(sens_delta, file.path(outdir, "source_data", "07_signature_sensitivity_deltas.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

universe_genes <- unique(rna_genes)
membership_sets <- list(Conserved = sets$Conserved, Rice_direct = sets$Rice_direct,
                        Arabidopsis_direct = sets$Arabidopsis_direct, Ultra_core = sets$Ultra_core)
overlap_stats <- function(genes, signature, universe = universe_genes) {
  genes <- unique(intersect(genes, universe)); sig <- unique(intersect(membership_sets[[signature]], universe))
  k <- length(intersect(genes, sig)); n <- length(genes); K <- length(sig); N <- length(universe)
  data.frame(signature = signature, list_size = n, overlap_n = k,
             overlap_fraction = ifelse(n > 0, k/n, 0),
             expected_fraction = K/N,
             fold_over_expected = ifelse(n > 0 && k > 0, (k/n)/(K/N), 0),
             enrichment_p = ifelse(n > 0, phyper(k - 1, K, N - K, n, lower.tail = FALSE), NA_real_))
}

markers <- fread(marker_path) %>%
  mutate(cluster = recode(cluster, `Defense-secondary metabolism cells` = "Sclerenchyma cells")) %>%
  filter(avg_log2FC > 0, cluster %in% as.character(unique(seu$celltype)))
marker_audit <- lapply(unique(markers$cluster), function(ct) {
  z <- markers %>% filter(cluster == ct) %>% arrange(desc(avg_log2FC))
  bind_rows(lapply(c(20L, 50L, nrow(z)), function(topn) {
    genes <- head(z$gene, topn)
    bind_rows(lapply(names(membership_sets), function(s) overlap_stats(genes, s))) %>%
      mutate(cluster = ct, scope = ifelse(topn >= nrow(z), "all_positive", paste0("top", topn)))
  }))
}) %>% bind_rows() %>% group_by(signature, scope) %>% mutate(enrichment_fdr = p.adjust(enrichment_p, "BH")) %>% ungroup()
write.table(marker_audit, file.path(outdir, "source_data", "08_cluster_marker_overlap_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

deg_files <- list.files(deg_dir, pattern = "_Significant_DEGs\\.csv$", full.names = TRUE)
deg_audit <- lapply(deg_files, function(f) {
  nm <- sub("_Significant_DEGs\\.csv$", "", basename(f))
  ct <- gsub("_", " ", nm, fixed = TRUE)
  if (ct == "Defense-secondary metabolism cells") ct <- "Sclerenchyma cells"
  if (!ct %in% as.character(unique(seu$celltype))) return(NULL)
  z <- fread(f) %>% arrange(desc(abs(logFC)))
  top50 <- bind_rows(lapply(names(membership_sets), function(s) overlap_stats(head(z$Symbol, 50L), s))) %>%
    mutate(cluster = ct, scope = "top50_abs_logFC")
  all_deg <- bind_rows(lapply(names(membership_sets), function(s) overlap_stats(z$Symbol, s))) %>%
    mutate(cluster = ct, scope = "all_reported")
  bind_rows(top50, all_deg)
}) %>% bind_rows() %>% group_by(signature, scope) %>% mutate(enrichment_fdr = p.adjust(enrichment_p, "BH")) %>% ungroup()
write.table(deg_audit, file.path(outdir, "source_data", "09_DEG_overlap_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cluster_map <- data.frame(seurat_cluster = as.character(seu$seurat_clusters), celltype = as.character(seu$celltype)) %>% distinct()
write.table(cluster_map, file.path(outdir, "source_data", "10_SCENIC_cluster_mapping.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
scenic <- fread(scenic_path)
scenic_audit <- lapply(seq_len(nrow(scenic)), function(i) {
  targets <- unlist(strsplit(scenic$All_Targets[i], ";", fixed = TRUE))
  bind_rows(lapply(names(membership_sets), function(s) {
    overlap_stats(targets, s) %>% mutate(TF_in_signature = scenic$TF_ID[i] %in% membership_sets[[s]])
  })) %>%
    mutate(Regulon = scenic$Regulon[i], TF_ID = scenic$TF_ID[i], Best_Cluster = as.character(scenic$Best_Cluster[i]),
           Max_RSS = scenic$Max_RSS[i], Target_Count_reported = scenic$Target_Count[i])
}) %>% bind_rows() %>% left_join(cluster_map, by = c("Best_Cluster" = "seurat_cluster")) %>%
  group_by(signature) %>% mutate(enrichment_fdr = p.adjust(enrichment_p, "BH")) %>% ungroup()
scenic_con <- gzfile(file.path(outdir, "source_data", "11_SCENIC_regulon_overlap_audit.tsv.gz"), open = "wt")
write.table(scenic_audit, scenic_con, sep = "\t", row.names = FALSE, quote = FALSE)
close(scenic_con)

go_files <- list.files(go_dir, pattern = "_GO_All_Results\\.csv$", full.names = TRUE)
go_audit <- lapply(go_files, function(f) {
  z <- fread(f)
  stem <- sub("_GO_All_Results\\.csv$", "", basename(f))
  direction <- ifelse(grepl("_up$", stem), "up", "down")
  ct <- sub("_(up|down)$", "", stem) %>% gsub("_", " ", .)
  if (!ct %in% as.character(unique(seu$celltype)) || nrow(z) == 0) return(NULL)
  z <- z %>% arrange(p.adjust) %>% mutate(rank_within_file = row_number())
  bind_rows(lapply(seq_len(nrow(z)), function(i) {
    genes <- unlist(strsplit(z$geneID[i], "/", fixed = TRUE))
    overlap_stats(genes, "Conserved") %>%
      mutate(cluster = ct, direction = direction, GO_ID = z$ID[i], Description = z$Description[i],
             ontology = z$Ontology[i], GO_p_adjust = z$p.adjust[i], GO_rank = z$rank_within_file[i])
  }))
}) %>% bind_rows() %>% mutate(overlap_enrichment_fdr = p.adjust(enrichment_p, "BH"))
write.table(go_audit, file.path(outdir, "source_data", "12_GO_pathway_overlap_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Compact headline table.
headline <- data.frame(
  metric = c("Cells", "Cell types", "Primary conserved genes detected", "Rice direct mapped genes detected",
             "Arabidopsis direct mapped genes detected", "Ultra-core genes detected",
             "Top-50 marker overlap median (%)", "Reported DEG overlap median (%)",
             "Top SCENIC regulons with TF in primary signature", "GO terms with >=25% primary-signature genes"),
  value = c(ncol(seu), length(unique(seu$celltype)), length(sets$Conserved), length(sets$Rice_direct),
            length(sets$Arabidopsis_direct), length(sets$Ultra_core),
            100 * median(marker_audit$overlap_fraction[marker_audit$signature == "Conserved" & marker_audit$scope == "top50"], na.rm = TRUE),
            100 * median(deg_audit$overlap_fraction[deg_audit$signature == "Conserved" & deg_audit$scope == "all_reported"], na.rm = TRUE),
            sum(scenic_audit$signature == "Conserved" & scenic_audit$TF_in_signature & scenic_audit$Max_RSS >= quantile(scenic$Max_RSS, .90, na.rm=TRUE)),
            sum(go_audit$overlap_fraction >= .25 & go_audit$overlap_n >= 2, na.rm = TRUE))
)
write.table(headline, file.path(outdir, "source_data", "13_headline_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ---------- Figure ----------
base_theme <- theme_classic(base_family = "Times New Roman", base_size = 9) +
  theme(axis.text = element_text(colour = "black"), axis.title = element_text(colour = "black"),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"), plot.margin = margin(5, 7, 5, 5))
pal_cond <- c(CK = "#3B6FB6", FCR = "#D97732")

lims <- quantile(cell_scores$conserved_score, c(.01, .99), na.rm = TRUE)
pA <- ggplot(cell_scores, aes(UMAP_1, UMAP_2, colour = conserved_score)) +
  geom_point(size = .18, alpha = .75) +
  scale_colour_gradientn(colours = c("#F4F1EA", "#5AB4AC", "#01665E", "#542788"), limits = lims, oob = squish,
                         name = "Conserved\nprotoplast score") +
  coord_equal() + labs(x = "UMAP 1", y = "UMAP 2") + base_theme +
  guides(colour = guide_colourbar(direction = "horizontal", barwidth = unit(34, "mm"), barheight = unit(2.8, "mm"))) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        legend.position = "bottom", legend.title = element_text(size = 8))

sum_plot <- score_summary %>% mutate(cluster = factor(cluster, levels = rev(celltype_order)))
pB <- ggplot(sum_plot, aes(median, cluster, colour = condition)) +
  geom_errorbar(aes(xmin = q1, xmax = q3), orientation = "y", width = 0, linewidth = .7, position = position_dodge(width = .45)) +
  geom_point(size = 1.9, position = position_dodge(width = .45)) +
  scale_colour_manual(values = pal_cond) +
  labs(x = "Median score (point; IQR line)", y = NULL, colour = NULL) + base_theme +
  theme(legend.position = "top", axis.text.y = element_text(size = 7.2))

delta_plot <- sens_delta %>% mutate(cluster = factor(cluster, levels = rev(celltype_order)))
mx <- max(abs(delta_plot$delta), na.rm = TRUE)
pC <- ggplot(delta_plot, aes(signature, cluster, fill = delta)) +
  geom_tile(colour = "white", linewidth = .4) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C44E52", midpoint = 0,
                       limits = c(-mx, mx), name = expression(Delta*" median")) +
  labs(x = NULL, y = NULL) + base_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), axis.text.y = element_text(size = 7.2))

ma <- marker_audit %>% filter(signature == "Conserved", scope == "top50") %>%
  mutate(cluster = factor(cluster, levels = rev(celltype_order)))
pD <- ggplot(ma, aes(100 * overlap_fraction, cluster)) +
  geom_vline(xintercept = 100 * unique(ma$expected_fraction)[1], linetype = 2, colour = "grey45", linewidth = .45) +
  geom_segment(aes(x = 0, xend = 100 * overlap_fraction, yend = cluster), colour = "#74A9CF", linewidth = .7) +
  geom_point(size = 2, colour = "#045A8D") +
  labs(x = "Conserved genes among top 50 markers (%)", y = NULL) + base_theme +
  theme(axis.text.y = element_text(size = 7.2))

da <- deg_audit %>% filter(signature == "Conserved", scope == "all_reported") %>%
  mutate(cluster = factor(cluster, levels = rev(celltype_order)))
pE <- ggplot(da, aes(100 * overlap_fraction, cluster)) +
  geom_vline(xintercept = 100 * unique(da$expected_fraction)[1], linetype = 2, colour = "grey45", linewidth = .45) +
  geom_segment(aes(x = 0, xend = 100 * overlap_fraction, yend = cluster), colour = "#F4A582", linewidth = .7) +
  geom_point(aes(fill = enrichment_fdr < .05), shape = 21, size = 2.2, colour = "#8C2D04") +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "white"), guide = "none") +
  labs(x = "Conserved genes among reported DEGs (%)", y = NULL) + base_theme +
  theme(axis.text.y = element_text(size = 7.2), legend.position = "right")

sc_top <- scenic_audit %>% filter(signature == "Conserved") %>%
  group_by(celltype) %>% slice_max(Max_RSS, n = 5, with_ties = FALSE) %>% ungroup() %>%
  arrange(desc(fold_over_expected), desc(overlap_n)) %>% slice_head(n = 10) %>%
  mutate(label = sub("\\(\\+\\)$", "", Regulon),
         label = factor(label, levels = rev(unique(label))))
pF <- ggplot(sc_top, aes(100 * overlap_fraction, label, size = pmax(overlap_n, 1), colour = fold_over_expected)) +
  geom_point(alpha = .9) +
  scale_colour_gradient(low = "#80CDC1", high = "#542788", name = "Fold over\nexpected") +
  scale_size_continuous(range = c(1.5, 4.5), guide = "none") +
  labs(x = "Primary-signature targets in top regulon (%)", y = NULL) + base_theme +
  theme(axis.text.y = element_text(size = 6.4))

go_top <- go_audit %>% filter(GO_rank <= 15, overlap_n >= 1) %>%
  arrange(desc(overlap_fraction), GO_p_adjust) %>% distinct(Description, .keep_all = TRUE) %>% slice_head(n = 10) %>%
  mutate(label = str_wrap(Description, 42),
         label = factor(label, levels = rev(unique(label))))
if (nrow(go_top) > 0) {
  pG <- ggplot(go_top, aes(100 * overlap_fraction, label, size = overlap_n, colour = -log10(pmax(GO_p_adjust, 1e-300)))) +
    geom_point(alpha = .9) +
    scale_colour_gradient(low = "#7FCDBB", high = "#2C3E75", name = expression(-log[10]*" GO FDR")) +
    scale_size_continuous(range = c(1.7, 5), guide = "none") +
    labs(x = "Primary-signature genes in enriched GO term (%)", y = NULL) + base_theme +
    theme(axis.text.y = element_text(size = 6.2))
} else {
  pG <- ggplot() + annotate("text", x = 0, y = 0, label = "No top GO terms contained primary-signature genes", family = "Times New Roman", size = 3) + theme_void()
}

figure <- ((pA | pB) + plot_layout(widths = c(.82, 1.18))) /
  ((pC | pD) + plot_layout(widths = c(.82, 1.18))) /
  ((pE | pF) + plot_layout(widths = c(1, 1))) /
  pG +
  plot_layout(heights = c(1.05, 1.18, 1.25, 1.2)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(family = "Times New Roman", face = "bold", size = 12)))

ggsave(file.path(outdir, "Supplementary_Figure_protoplast_audit.pdf"), figure,
       width = 183, height = 255, units = "mm", device = cairo_pdf)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(outdir, "Supplementary_Figure_protoplast_audit.svg"), figure,
         width = 183, height = 255, units = "mm", device = svglite::svglite)
}
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(file.path(outdir, "Supplementary_Figure_protoplast_audit.tiff"), figure,
         width = 183, height = 255, units = "mm", dpi = 600, compression = "lzw", device = ragg::agg_tiff)
  ggsave(file.path(outdir, "Supplementary_Figure_protoplast_audit_preview.png"), figure,
         width = 183, height = 255, units = "mm", dpi = 200, device = ragg::agg_png)
}

caption <- paste0(
  "Supplementary Fig. Sx | Audit of protoplasting-associated transcriptional signals in the wheat single-cell atlas. ",
  "(A) UMAP coloured by a rank-based UCell score computed from ", length(sets$Conserved),
  " detected wheat genes supported by independent rice and Arabidopsis intact-tissue-versus-protoplast experiments after strict reciprocal protein mapping. ",
  "(B) Median score and interquartile range within each annotated cell population in CK and FCR. ",
  "(C) FCR-minus-CK differences in median score for the conserved signature and the 10-gene ultra-core rice signature. ",
  "(D) Fraction of the 50 strongest positive markers per population belonging to the conserved signature; the dashed line denotes the transcriptome-wide expected fraction. ",
  "(E) Fraction of reported CK-FCR DEGs belonging to the conserved signature; colour represents BH-adjusted hypergeometric overlap evidence and the dashed line denotes expectation. These overlaps are a sensitivity audit of the supplied DEG lists and do not provide biological-replicate inference. ",
  "(F) Conserved-signature target representation in the most cell-population-specific SCENIC regulons (five highest RSS regulons per population; the 18 largest enrichments shown). ",
  "(G) Conserved-signature representation in leading enriched GO terms. Cell-level score intervals and 500 cell-stratified bootstrap iterations quantify within-dataset stability only; they are not confidence intervals for population-level treatment effects."
)
writeLines(caption, file.path(outdir, "Supplementary_Figure_protoplast_audit_legend.txt"))

methods <- paste0(
  "Protoplasting-associated genes were curated from direct bulk comparisons of intact tissues and isolated protoplasts in rice (Zhu et al., 2025) and Arabidopsis (Liew et al., 2024). ",
  "Only induced genes were used. Published protein sequences were mapped to IWGSC wheat proteins by reciprocal-best-hit BLASTP; the rice mapping required E <= 1e-50, amino-acid identity >=70%, query and subject coverage >=70%, and retained A/B/D co-orthologues within 97% of the best bit score. ",
  "The Arabidopsis mapping required E <= 1e-20, identity >=40%, query and subject coverage >=60%, reciprocal recovery of the source gene, and the same 97% bit-score rule. ",
  "The primary signature comprised wheat genes supported by both mappings. Per-cell scores were calculated from RNA counts using the rank-based UCell method (maxRank=3000). ",
  "Scores were summarized within each annotated population and condition. FCR-CK differences were reported as descriptive effect sizes; 500 bootstrap iterations resampled cells separately within each population and condition to assess sensitivity to cell composition and outlying cells, without treating cells as biological replicates. ",
  "Overlap of the primary signature with cluster markers, supplied CK-FCR DEG lists, SCENIC regulon targets, and enriched GO-term gene sets was assessed against the detected RNA-gene universe using one-sided hypergeometric tests followed by Benjamini-Hochberg correction."
)
writeLines(methods, file.path(outdir, "Methods_protoplast_audit.txt"))

saveRDS(list(gene_sets = sets_raw, detected_gene_sets = sets, score_columns = score_cols,
             score_summary = score_summary, bootstrap_delta = boot_delta,
             marker_audit = marker_audit, deg_audit = deg_audit,
             scenic_audit = scenic_audit, go_audit = go_audit),
        file.path(outdir, "protoplast_audit_analysis_objects.rds"), compress = "xz")

writeLines(c(
  "Reproducible protoplasting-signal audit; random seed 54.",
  "The original merged.rds was read only and was not modified or overwritten.",
  "Primary scoring uses the cross-study conserved signature; broad species-specific sets are retained for overlap sensitivity analyses.",
  "CK-FCR differences are descriptive because cells are not independent biological replicates.",
  "The bootstrap quantifies stability to cell resampling within the observed dataset only."
), file.path(outdir, "README.txt"))

input_files <- c(
  "rice_protoplast_gene_sets.tsv",
  "high_confidence_rice_to_wheat_protoplast_RBH.tsv",
  "arabidopsis_embryo_protoplast_genes.tsv",
  "arabidopsis_to_wheat_protoplast_strict_RBH.tsv"
)
file.copy(file.path(set_dir, input_files), file.path(outdir, "source_data", "input_gene_sets", input_files), overwrite = TRUE)

file.copy(file.path(root, "protoplast_audit.R"), file.path(outdir, "protoplast_audit.R"), overwrite = TRUE)
message("Completed: ", outdir)
