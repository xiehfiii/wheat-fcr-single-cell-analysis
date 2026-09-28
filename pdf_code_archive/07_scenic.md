# SCENIC gene-regulatory network

PDF transcription only. The PDF may wrap code or omit glyphs;
use the curated scripts for executable analysis.

## Source PDF page 91

```text
if(my_tf %in% gene_cluster_df$GeneID){
print(paste(" ", my_tf, "   Cluster:",
gene_cluster_df[gene_cluster_df$GeneID == my_tf, "Cluster"]))
} else {
print(paste(" ", my_tf, "    qval                   "))
}
write.csv(gene_cluster_df, "Epidermis_BEAM.csv", row.names = FALSE)



#######################################################################
#########
# 1.     Fate 1 (     )
#        hm
left_means <- rowMeans(hm$heatmap_matrix[, 1:round(ncol(hm$heatmap_matrix)/3)])
target_gene_id <- names(which.max(left_means))

# 2.              State
#  Monocle
p_check <- plot_genes_jitter(mycds[target_gene_id, ],
grouping = "State",
color_by = "State")
print(p_check)

# 3.
avg_expr_by_state <- tapply(exprs(mycds)[target_gene_id, ], pData(mycds)$State, mean)
best_state <- names(which.max(avg_expr_by_state))

cat("\n================ ================\n")
cat("       Fate 1       :", target_gene_id, "\n")
cat("           State          :\n")
print(avg_expr_by_state)
cat("\n     Fate 1 (     )        State", best_state, "\n")
cat("     Fate 2 (     )                   State\n")
cat("======================================\n")
```

## Source PDF page 93

```text
mamba create -n scenic python=3.9 "numpy<1.24" pandas pyarrow python-flatbuffers cython cmake -c conda-forge -c bioconda --channel-priority flexible -y && mamba run -n scenic pip install pyscenic loompy scanpy ctxcore tqdm -i https://pypi.tuna.tsinghua.edu.cn/simple --default-timeout=1000

mamba activate scenic
pyscenic -h

wget https://planttfdb.gao-lab.org/download/motif/Tae_TF_binding_motifs_information.txt
mv Tae_TF_binding_motifs_information.txt motif.tsv
wget https://planttfdb.gao-lab.org/download/motif/Tae_TF_binding_motifs.meme.gz
gunzip Tae_TF_binding_motifs.meme.gz
wget https://planttfdb.gao-lab.org/download/TF_list/Tae_TF_list.txt.gz
gunzip Tae_TF_list.txt.gz
mv Tae_TF_list.txt wheat_tfs.txt

vim extract_promoter.sh

#!/bin/bash

# =================    =================
GTF="Triticum_aestivum.gtf"
GENOME="Triticum_aestivum.fasta"
SIZES="genome.sizes"
# ===========================================

#
if [ ! -f "$GTF" ]; then echo "   GTF: $GTF"; exit 1; fi
if [ ! -f "$GENOME" ]; then echo "   Fasta: $GENOME"; exit 1; fi

# -------------------------------------------------
# 1.              (        )
# -------------------------------------------------
if [ ! -s "$SIZES" ]; then
echo "1.             ..."
if [ -f "${GENOME}.fai" ]; then
cut -f1,2 "${GENOME}.fai" > "$SIZES"
else
awk 'BEGIN {OFS="\t"} /^>/ {$1=$1; print substr($1,2), 0; next} {len[substr($1,2)]+=length($0)} END {for(c in len) print c, len
```

## Source PDF page 94

```text
fi
else
echo "1.        $SIZES        "
fi

# -------------------------------------------------
# 2.         (          )
# -------------------------------------------------
echo "2.    GTF          ..."

#
# 1.      $3 == "transcript" (      )
# 2.   gsub(/^gene:/, "", val)  gene:          ID
awk 'BEGIN {OFS="\t"} $3 == "transcript" {
gene_id = "";
for(i=9;i<=NF;i++) {
if ($i == "gene_id") {
val = $(i+1);
gsub(/[";]/, "", val); #
gsub(/^gene:/, "", val); #  gene:    (  gene:Traes... -> Traes...)
gene_id = val;
break;
}
}

#    BED    : chr, start-1, end, name, score, strand
if (gene_id != "") {
print $1, $4 - 1, $5, gene_id, ".", $7;
}
}' "$GTF" > genes_raw.bed

#                            ID
sort -k4,4 -u genes_raw.bed > genes.bed
rm genes_raw.bed

echo "           : $(wc -l < genes.bed)"
if [ ! -s genes.bed ]; then echo " : genes.bed           awk"; exit 1; fi

# -------------------------------------------------
# 3.   Promoter    (    3000bp)
# -------------------------------------------------
echo "3.             (-l 3000)..."
bedtools flank -i genes.bed -g "$SIZES" -l 3000 -r 0 -s > promoters.bed
```

## Source PDF page 95

```text

# -------------------------------------------------
# 4.
# -------------------------------------------------
echo "4.      Fasta   ..."
bedtools getfasta -fi "$GENOME" -bed promoters.bed -s -name -fo wheat.promoter3K.V2.fa

echo "              : wheat.promoter3K.V2.fa"

bash extract_promoter.sh
raw_data <- read_tsv("motif.tsv")

motif <- raw_data %>%
select(Matrix_id, Gene_id, Datasource) %>%
rename(
motif = Matrix_id,
TF = Gene_id,
source = Datasource
)

head(motif)

motif2TF <- motif %>% dplyr::transmute(`#motif_id`=motif, motif_name=motif, motif_description=TF,
source_name=source, source_version=1.1, gene_name=TF,
motif_similarity_qvalue=0.000000, similar_motif_id="None",
similar_motif_description="None", orthologous_identity=1.000000,
orthologous_gene_name="None", orthologous_species="None",
description="gene is directly annotated")
write.table(motif2TF, file = "wheat_motif2TF.tbl", sep = "\t", row.names = F, quote = F)
```

## Source PDF page 96

```text
lines <- readLines("Tae_TF_binding_motifs.meme")
output_lines <- c()
current_motif <- NULL
reading_matrix <- FALSE
motif_matrix <- c()

for (line in lines) {
if (startsWith(line, "MOTIF")) {
if (!is.null(current_motif) && length(motif_matrix) > 0) {
output_lines <- c(output_lines,
paste0(">", current_motif),
motif_matrix)
}
parts <- strsplit(line, " +")[[1]]
current_motif <- parts[length(parts)]
reading_matrix <- FALSE
motif_matrix <- c()
} else if (grepl("letter-probability matrix", line)) {
reading_matrix <- TRUE
} else if (reading_matrix) {
if (trimws(line) == "" || startsWith(line, "URL") || startsWith(line, "MOTIF")) {
reading_matrix <- FALSE
} else {
motif_matrix <- c(motif_matrix, line)
}
}
}
if (!is.null(current_motif) && length(motif_matrix) > 0) {
output_lines <- c(output_lines,
paste0(">", current_motif),
motif_matrix)
}

#
writeLines(output_lines, con = "wheat.TF.matrix.txt")
#  .cb    motif
awk '
/^>/ {
if (out) close(out)
out = substr($0, 2) ".cb"
print $0 > out
next
```

## Source PDF page 97

```text
}
{
print > out
}
' wheat.TF.matrix.txt

mkdir motif_cb_format
mv *.cb motif_cb_format
#  motif list
grep '^>' wheat.TF.matrix.txt | sed 's/^>//' > wheat.TF.txt

nohup python create_cistarget_motif_databases.py -f wheat.promoter3K.V2.fa -M motif_cb_format/ -m wheat.TF.txt -o wheat.genes_vs_motifs -t 3
&

mkdir cisTarget_databases
mv wheat_tfs.txt wheat_motif2TF.tbl wheat.genes_vs_motifs.ranking.feather cisTarget_databases
```

## Source PDF page 98

```text
cat > auto_filter.R <<EOF
library(Seurat)

# 1.
input_file <- "merged.rds"
if(!file.exists(input_file)) stop(" merged.rds   ")

print(paste0("[", Sys.time(), "]  : ", input_file))
seurat_obj <- readRDS(input_file)

# 2.   TF     (             )
#       TF      cisTarget_databases
tf_file <- "cisTarget_databases/wheat_tfs.txt"
valid_tfs <- c()

if(file.exists(tf_file)){
print("      TF             ...")
#            ID
tf_list <- read.table(tf_file, header=F)\$V1
valid_tfs <- intersect(rownames(seurat_obj), tf_list)
print(paste("        TF   :", length(valid_tfs)))
} else {
print("        wheat_tfs.txt             ")
}

# 3.         (      )
print(paste0("[", Sys.time(), "]             ..."))
#
seurat_obj <- NormalizeData(seurat_obj)
#    20000        (                          )
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 20000)
top_hvgs <- VariableFeatures(seurat_obj)

# 4.         (       +    TF)
keep_genes <- unique(c(top_hvgs, valid_tfs))
print(paste("        :", length(keep_genes)))

seurat_final <- seurat_obj[keep_genes, ]

# 5.
output_file <- "merged_final_20k.rds"
print(paste0("[", Sys.time(), "]           : ", output_file))
```

## Source PDF page 99

```text
saveRDS(seurat_final, output_file)
print("R           ")
EOF



cat > AutoRun.sh <<EOF
#!/bin/bash

echo "============================================="
echo "      SCENIC    (    +    )   "
echo "============================================="

# 1.   R
echo ">>>   1:       R              ..."
Rscript auto_filter.R

#    R
if [ \$? -eq 0 ]; then
echo "  R                RDS         "
echo "---------------------------------------------"
echo ">>>   2:     SCENIC GRN  ..."

# 2.    RunGRN.sh
#                    merged_final_20k.rds
#         12
bash RunGRN.sh -f merged_final_20k.rds -s wheat -t 12 -d cisTarget_databases -o output

else
echo "     : R                        merged.rds       "
exit 1
fi
EOF

nohup bash AutoRun.sh &
```

## Source PDF page 100

```text
install.packages("optparse")
install.packages("Seurat")
devtools::install_github("aertslab/SCopeLoomR")
install.packages("data.table")
install.packages("pbapply")
install.packages("philentropy")
install.packages("dplyr")
library(scPlant)
library(optparse)
library(data.table)
library(pbapply)
library(philentropy)
library(dplyr)
library(pals)
library(ComplexHeatmap)
library(grid)
library(circlize)
library(pheatmap)
library(ComplexHeatmap)
library(CellFunTopic)
library(cowplot)
#          Seurat
rasMat <- readRDS("output/rasMat.rds")
rssMat <- readRDS("output/rssMat.rds")
tf_target <- readRDS("output/tf_target.rds")
merged <- readRDS("merged.rds")
#           debug
#            regulon     TF
ras_exp_hmp_fixed <- function(merged, rasMat, group.by = "seurat_clusters", assay = "SCT",
filename = NULL, width = 10, height = 8){

#
log_msg <- function(msg) message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg))

log_msg("Step 1/8:      ...")
require(Seurat)
require(ComplexHeatmap)
require(circlize)

# 1.    group.by
```

## Source PDF page 101

```text
if (!group.by %in% colnames(merged@meta.data)) {
stop(paste0("           '", group.by, "'"))
}

# 2.                    (           0, 1, 2... 14  )
current_groups <- merged@meta.data[[group.by]]
if (!is.factor(current_groups)) {
log_msg("-->                     ...")
#
num_levels <- sort(as.numeric(unique(as.character(current_groups))))
if (any(is.na(num_levels))) {
#
merged@meta.data[[group.by]] <- factor(current_groups, levels = sort(unique(current_groups)))
} else {
#
merged@meta.data[[group.by]] <- factor(current_groups, levels = num_levels)
}
}
#
target_order <- levels(merged@meta.data[[group.by]])
log_msg(paste0("-->    : ", paste(head(target_order), collapse=","), " ..."))

# 3.       rasMat
if (ncol(rasMat) == ncol(merged) || any(colnames(rasMat) %in% colnames(merged))) {
rasMat <- t(rasMat)
}

# 4.
common_cells <- intersect(rownames(merged@meta.data), rownames(rasMat))
if(length(common_cells) == 0) stop("        ")
merged_sub <- subset(merged, cells = common_cells)
rasMat_sub <- rasMat[common_cells, , drop=FALSE]

# 5.    RAS
log_msg("Step 4/8:  RAS    ...")
groups <- merged_sub@meta.data[[group.by]]
meanRAS <- t(sapply(split(rownames(merged_sub@meta.data), groups), function(cells) {
colMeans(rasMat_sub[cells, , drop=FALSE])
}))
```

## Source PDF page 102

```text
meanRAS_scaled <- scale(meanRAS)

# 6.           (    Seurat V5)
log_msg("Step 5/8:            ...")
merged_sub <- SetIdent(merged_sub, value = groups)
meanExp_obj <- tryCatch({
AverageExpression(merged_sub, assays = assay, slot = "data")[[assay]]
}, error = function(e) stop("AverageExpression       Assay"))

# V5 "g"
log_msg("Step 6/8:       ...")
ras_groups <- rownames(meanRAS_scaled)
exp_groups <- colnames(meanExp_obj)

if (!all(ras_groups %in% exp_groups)) {
cleaned_exp_groups <- gsub("^g", "", exp_groups)
if (all(ras_groups %in% cleaned_exp_groups)) {
colnames(meanExp_obj) <- cleaned_exp_groups
log_msg("-->    Seurat V5  'g'   ")
} else {
#
common_g <- intersect(ras_groups, exp_groups)
if(length(common_g)>0) {
meanRAS_scaled <- meanRAS_scaled[common_g, , drop=FALSE]
meanExp_obj <- meanExp_obj[, common_g, drop=FALSE]
}
}
}

# 7.    TF-Gene
log_msg("Step 7/8:  TF   Gene...")
regulon_names <- colnames(meanRAS_scaled)
clean_tf_names <- gsub("\\(\\+\\)", "", regulon_names)
clean_tf_names <- gsub("g$", "", clean_tf_names)

valid_idx <- which(clean_tf_names %in% rownames(meanExp_obj))

if (length(valid_idx) == 0) {
warning("     TF-Gene      RAS ")
matched_regulons <- regulon_names
meanExp_plot <- matrix(0, nrow=nrow(meanRAS_scaled), ncol=length(matched_regulons))
} else {
matched_regulons <- regulon_names[valid_idx]
```

## Source PDF page 103

```text
matched_genes <- clean_tf_names[valid_idx]

meanExp <- t(meanExp_obj[matched_genes, rownames(meanRAS_scaled)])
meanRAS_plot <- meanRAS_scaled[, matched_regulons]
meanExp_plot <- scale(meanExp)

meanRAS_scaled <- meanRAS_plot
}

# 8.   Target Order       (    0,1...14   )
valid_target_order <- target_order[target_order %in% rownames(meanRAS_scaled)]
meanRAS_scaled <- meanRAS_scaled[valid_target_order, , drop=FALSE]
meanExp_plot <- meanExp_plot[valid_target_order, , drop=FALSE]

#
mat_ras <- t(meanRAS_scaled)
mat_exp <- t(meanExp_plot)

#
ras_max <- max(abs(mat_ras), na.rm=T); if(ras_max==0) ras_max<-1
exp_max <- max(abs(mat_exp), na.rm=T); if(exp_max==0) exp_max<-1

col_fun1 = colorRamp2(c(-ras_max, 0, ras_max), c("navy", "white", "firebrick3"))
col_fun2 = colorRamp2(c(-exp_max, 0, exp_max), c("white", "yellow","purple"))

# ===       ===
log_msg("Step 8/8:         ...")

ht1 <- Heatmap(mat_ras, name = "Regulon Activity", col = col_fun1,
cluster_rows = TRUE,
cluster_columns = FALSE, #
show_row_names = FALSE, # <---
column_title = "Regulon Activity",
border = TRUE)

ht2 <- Heatmap(mat_exp, name = "TF Expression", col = col_fun2,
cluster_rows = TRUE,
cluster_columns = FALSE,
show_row_names = FALSE,
column_title = "TF Expression",
border = TRUE)

ht_list <- ht1 + ht2
```

## Source PDF page 104

```text

if (!is.null(filename)) {
if (grepl("\\.png$", filename, ignore.case = TRUE)) {
png(filename, width = width, height = height, units = "in", res= 300)
} else {
pdf(filename, width = width, height = height)
}
draw(ht_list, ht_gap = unit(1, "cm"))
dev.off()
log_msg(paste0("    : ", filename))
} else {
draw(ht_list, ht_gap = unit(1, "cm"))
log_msg("   ")
}
}

#        PDF
ras_exp_hmp_fixed(merged, rasMat, group.by = "seurat_clusters",
filename = "TF_Heatmap.pdf")

ras_exp_scatter_fixed <- function(merged, rasMat, gene, reduction = 'umap',
assay = "SCT",
filename = NULL, width = 12, height = 6) {
pure_gene_name <- gsub("\\(\\+\\)$", "", gene)
pure_gene_name <- gsub("g$", "", pure_gene_name)

#         Regulon    (       )
target_regulon_name <- paste0(pure_gene_name, "(+)")

message(paste0("    :   [", pure_gene_name, "] | Regulon[", target_regulon_name, "]"))

# === 2.    RAS    ===
#
if (ncol(rasMat) == ncol(merged) || any(colnames(rasMat) %in% colnames(merged))) {
rasMat <- t(rasMat)
}

# === 3.    Regulon      ===
if (!target_regulon_name %in% colnames(rasMat)) {
#
hits <- grep(pure_gene_name, colnames(rasMat), value = TRUE)
if (length(hits) > 0) {
```

## Source PDF page 105

```text
target_regulon_name <- hits[1]
message(paste0("          Regulon: ", target_regulon_name))
} else {
stop(paste0("     rasMat      Regulon '", target_regulon_nam。e, "' "))
}
}

# === 4.            (   Seurat ) ===
#        Assay
if (!pure_gene_name %in% rownames(merged)) {
#      Assay
found <- FALSE
for (a in Assays(merged)) {
if (pure_gene_name %in% rownames(merged[[a]])) {
message(paste0("        Assay '", a, "'           Assay..."))
DefaultAssay(merged) <- a
found <- TRUE
break
}
}
if (!found) stop(paste0("   Seurat        Assay。'", pure_gene_name, "' "))
}

common_cells <- intersect(colnames(merged), rownames(rasMat))
if (length(common_cells) == 0) stop("    ")
sub_obj <- subset(merged, cells = common_cells)
ras_values <- rasMat[common_cells, target_regulon_name]
sub_obj <- AddMetaData(sub_obj, metadata = ras_values, col.name = "RAS_Activity")
p1 <- FeaturePlot(sub_obj, features = pure_gene_name, reduction = reduction, order = TRUE) +
scale_color_gradient(low = "lightgrey", high = "blue") +
ggtitle(paste0("Expression: ", pure_gene_name)) +
theme(plot.title = element_text(hjust = 0.5))

p2 <- FeaturePlot(sub_obj, features = "RAS_Activity", reduction = reduction, order = TRUE) +
scale_color_gradient(low = "lightgrey", high = "red") +
ggtitle(paste0("Activity: ", target_regulon_name)) +
theme(plot.title = element_text(hjust = 0.5))

p_merged <- plot_grid(p1, p2, ncol = 2)
if (!is.null(filename)) {
```

## Source PDF page 106

```text
ggsave(filename = filename, plot = p_merged, width = width, height= height)
message(paste0("    : ", filename))
}
return(p_merged)
}

ras_exp_scatter_fixed(merged, rasMat, gene = 'TraesCS5A02G237900', reduction = 'umap', filename = "TraesCS5A02G237900_scatter.pdf")



######### 8top10
plot_rss_rank <- function(rssMat, cluster_id, topn = 10) {
rss_df <- rssMat %>%
as.data.frame() %>%
rownames_to_column("Regulon") %>%
select(Regulon, RSS = as.character(cluster_id))
rss_df$RSS <- as.numeric(rss_df$RSS)
plot_data <- rss_df %>%
arrange(desc(RSS)) %>%
slice_head(n = topn) %>%
mutate(Regulon = factor(Regulon, levels = rev(Regulon))) #

p <- ggplot(plot_data, aes(x = RSS, y = Regulon)) +
geom_segment(aes(x = 0, xend = RSS, y = Regulon, yend = Regulon), color = "grey") +
geom_point(size = 4, color = "#E41A1C") + #
geom_text(aes(label = round(RSS, 3)), hjust = -0.3, size = 3) + #
labs(title = paste0("Cluster ", cluster_id, " Top ", topn, " Regulons"),
x = "RSS Score", y = "") +
theme_minimal() +
theme(plot.title = element_text(hjust = 0.5, face = "bold"))

return(p)
}
p_rank <- plot_rss_rank(rssMat, cluster_id = "8", topn = 10)
print(p_rank)

ggsave("RSS_Rank_Cluster8.pdf", p_rank, width = 12, height = 8)
```

## Source PDF page 107

```text
ras_exp_hmp_fixed <- function(merged, rasMat, group.by = "celltype", assay = "SCT") {
meanRAS <- sapply(
split(rownames(merged@meta.data), merged@meta.data[[group.by]]),
function(cells) colMeans(rasMat[cells, , drop = FALSE])
)
meanRAS <- t(scale(t(meanRAS), center = TRUE, scale = TRUE))
ph <- pheatmap::pheatmap(meanRAS, cluster_rows = TRUE, cluster_cols= TRUE, silent = TRUE)
rowOrder <- ph$tree_row$labels[ph$tree_row$order]
colOrder <- ph$tree_col$labels[ph$tree_col$order]
Seurat::Idents(merged) <- merged@meta.data[[group.by]]
meanExp <- Seurat::AverageExpression(merged, assays = assay, slot ="data")[[assay]]
colnames(meanExp) <- sub("^[^-]+-", "", colnames(meanExp))
meanExp <- meanExp[colnames(rasMat), , drop = FALSE]
meanExp <- t(scale(t(meanExp), center = TRUE, scale = TRUE))
ht1 <- Heatmap(
meanRAS[rowOrder, colOrder],
name = "Regulon Activity",
col = colorRamp2(c(-2, 0, 2), c("white", "orange", "red")),
cluster_rows  = FALSE,
cluster_columns = FALSE
)
ht2 <- Heatmap(
meanExp[rowOrder, colOrder],
name = "TF Expression",
col = colorRamp2(c(-2, 0, 2), c("white", "blue", "darkblue")),
cluster_rows  = FALSE,
cluster_columns = FALSE
)
ht_list <- ht1 + ht2
grid.newpage()
draw(
ht_list,
newpage = FALSE,
ht_gap  = unit(1, "cm"),
heatmap_legend_side = "left"
)
}
ras_exp_hmp_fixed(merged, rasMat, group.by = "celltype", assay = "SCT")


```

## Source PDF page 108

```text



ras_exp_hmp_fixed <- function (merged, rasMat, group.by = "seurat_clusters", assay = "SCT")
{
meanRAS <- sapply(split(rownames(merged@meta.data), merged@meta.data[[group.by]]),
function(cells) colMeans(rasMat[cells, ]))
meanRAS <- t(scale(t(meanRAS), center = T, scale = T))
col_fun1 = circlize::colorRamp2(seq(min(meanRAS), max(meanRAS),
length.out = 10), c("white", pals::brewer.orrd(9)))
ph <- pheatmap::pheatmap(meanRAS, cluster_rows = T, cluster_cols =T,
silent = T)
rowOrder <- ph$tree_row$labels[ph$tree_row$order]
colOrder <- ph$tree_col$labels[ph$tree_col$order]
Seurat::Idents(merged) <- merged@meta.data[[group.by]]
meanExp <- Seurat::AverageExpression(merged, assays = assay,
slot = "data")[[assay]]
meanExp <- meanExp[colnames(rasMat), ] %>% as.matrix()
colnames(meanExp) <- gsub("-", "_", sub("^g", "", colnames(meanExp)))
meanExp <- t(scale(t(meanExp), center = T, scale = T))
col_fun2 = circlize::colorRamp2(seq(min(meanExp), max(meanExp),
length.out = 10), c("white", pals::brewer.blues(9)))
ht1 <- Heatmap(meanRAS[rowOrder, colOrder], name = "Regulon Activity",
col = col_fun1, column_names_rot = 45, row_names_gp = gpar(fontsize = 1),
cluster_rows = F, cluster_columns = F, column_names_gp = gpar(fontsize = 7),
column_names_centered = F, border = TRUE)
ht2 <- Heatmap(meanExp[rowOrder, colOrder], name = "TF Expression",
col = col_fun2, column_names_rot = 45, row_names_gp = gpar(fontsize = 1),
cluster_rows = F, cluster_columns = F, column_names_gp = gpar(fontsize = 7),
column_names_centered = F, border = TRUE)
ht_list = ht1 + ht2
draw(ht_list, ht_gap = unit(1, "cm"), heatmap_legend_side = "left")
}
#
ras_exp_hmp_fixed(merged, rasMat, group.by = "celltype", assay = "SCT")
```

## Source PDF page 109

```text
library(tidyverse)
library(scales)
library(gtools)

rss_mat <- readRDS("output/rssMat.rds")
targets <- c(
"TraesCS4A02G001300", "TraesCS5D02G491600", "TraesCS2B02G320700",
"TraesCS5B02G490900", "TraesCS1A02G076200", "TraesCS5A02G237900",
"TraesCS5B02G490700", "TraesCS7D02G497400", "TraesCS5D02G411400"
)

matched_names <- grep(paste(targets, collapse = "|"), rownames(rss_mat), value = TRUE)
df <- rss_mat[matched_names, , drop = FALSE] %>%
as.data.frame() %>%
rownames_to_column("Regulon") %>%
pivot_longer(-Regulon, names_to = "Cluster", values_to = "RSS")

df$Cluster <- factor(df$Cluster, levels = gtools::mixedsort(unique(df$Cluster)))
cluster_levels <- levels(df$Cluster)
order_df <- df %>%
group_by(Regulon) %>%
summarise(
max_rss = max(RSS),
max_cluster_idx = which(cluster_levels == Cluster[which.max(RSS)])[1]
) %>%
arrange(max_cluster_idx, desc(max_rss)) # Cluster     Cluster

df$Regulon <- factor(df$Regulon, levels = rev(order_df$Regulon)) # rev

min_val <- min(df$RSS)
max_val <- max(df$RSS)
cols <- c("#313695", "#4575B4", "#FFFFBF", "#D73027", "#A50026")

p <- ggplot(df, aes(x = Cluster, y = Regulon, fill = RSS)) +
geom_tile(color = "white", size = 0.5) +
scale_fill_gradientn(
colors = cols,
#     rescale                    (min  max)
```

## Source PDF page 110

```text
values = scales::rescale(c(min_val, (min_val+max_val)/2, max_val)),
limits = c(min_val, max_val),
name = "RSS Score"
) +
theme_minimal(base_size = 15) +
theme(
axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, color = "black", face = "bold"),
axis.text.y = element_text(color = "black", size = 12),
axis.title = element_blank(),
legend.position = "right",
legend.title = element_text(size = 12),
legend.text = element_text(size = 10),
panel.grid = element_blank()
) +
coord_fixed(ratio = 1)

ggsave("RSS_Heatmap.pdf", p, width = 10, height = 5)
```

## Source PDF page 111

```text
rss_mat <- readRDS("output/rssMat.rds")
message(">>>       Clusters: ")
print(colnames(rss_mat))

target_cluster <- "1"

if (target_cluster %in% colnames(rss_mat)) {
message(paste0("\n>>> [Cluster ", target_cluster, "]    Top 2
Regulons:"))
cluster_rss <- rss_mat[, target_cluster]
top_regulators <- sort(cluster_rss, decreasing = TRUE)[1:20]
print(data.frame(RSS_Score = round(top_regulators, 3)))
} else {
message(paste0("    :        '", target_cluster, "'Clusters    "))
}

target_tf <- "TraesCS6A02G120600"

if (target_tf %in% rownames(rss_mat)) {
message(paste0("\n>>> [", target_tf, "]    RSS        :"))
tf_rss <- rss_mat[target_tf, ]
sorted_tf_rss <- sort(tf_rss, decreasing = TRUE)
print(round(sorted_tf_rss, 3))
best_cluster <- names(sorted_tf_rss)[1]
best_score <- sorted_tf_rss[1]
message(paste0("\n>>> :   TF         Cluster [", best_cluster,"] (RSS=", round(best_score, 3), ")"))
} else {
message(paste0("    :           Regulon: ", target_tf))
message("  :                    TF                  ")
}
```

## Source PDF page 112

```text
library(dplyr)
library(stringr)

# 1.
target_tf <- "TraesCS4A02G001300"
target_gene <- "TraesCS6A02G018200"

print(paste("   :", target_tf, "->", target_gene))

# 2.      (     3    )
# fill=TRUE               , quote=""
df <- read.table("output/reg.tsv", skip = 3, sep = "\t", header = FALSE,
stringsAsFactors = FALSE, quote = "", comment.char ="", fill = TRUE)

# 3.            (             )
#  1   TF   2   MotifID  9       TargetGenes (    )
colnames(df)[1] <- "TF"
colnames(df)[2] <- "MotifID"

# 4.               (            )
tf_subset <- df %>% filter(TF == target_tf)

if(nrow(tf_subset) == 0) {
print("                             (TF)      TF ID")
} else {
print(paste("   TF   Regulon   :", nrow(tf_subset)))

# 5.
#             "[('GeneA', 1.0), ...]"
#                                ID

#
# fixed=TRUE
has_target <- apply(tf_subset, 1, function(x) any(grepl(target_gene, x, fixed = TRUE)))

#
final_result <- tf_subset[has_target, ]

if(nrow(final_result) > 0) {
```

## Source PDF page 113

```text
print("         ")
print("--------------------------------------------------")
print(paste("TF:", target_tf))
print(paste("Target:", target_gene))
print("--------------------------------------------------")
print("     Motif ID  ")
print(unique(final_result$MotifID))

#        Motif URL (           11    )
if(ncol(final_result) >= 11) {
print("Motif      ")
print(unique(final_result[, 11]))
}
} else {
print("         TF                             ")
print("     SCENIC                          ")
}
}
grep -A 20 "MP00290" Tae_TF_binding_motifs.meme
MOTIF TraesCS2A02G413900 MP00290

letter-probability matrix: alength= 4 w= 8 nsites= 158 E= 6.1e-196
0.164557      0.037975      0.126582      0.670886
0.000000      1.000000      0.000000      0.000000
0.000000      1.000000      0.000000      0.000000
0.000000      0.000000      1.000000      0.000000
0.000000      0.000000      0.000000      1.000000
1.000000      0.000000      0.000000      0.000000
0.000000      1.000000      0.000000      0.000000
0.430380      0.183544      0.208861      0.177215
T C C G T A C A
CCGTAC
```

## Source PDF page 114

```text
library(Biostrings)
promoters <- readDNAStringSet("wheat.promoter3K.V2.fa")

# 2.
#    bedtools -name   header      "GeneID::chr:start-end(strand)"
#     grepl
target_id <- "TraesCS4A02G001300"
idx <- grep(target_id, names(promoters))

if(length(idx) == 0) {
stop("    Fasta                    ID          ")
}

seq <- promoters[[idx]]
cat("                    :", length(seq), "bp\n")
# 1.              Motif
# M    A   C
refined_motif <- "MGCCGCC"

# 2.        3000bp
#        fixed = FALSE       IUPAC    M
matches <- matchPattern(refined_motif, seq, fixed = FALSE)

# 3.
cat("\n===        [A/C]GCCGCC      ===\n")
if(length(matches) > 0) {
print(matches)

#                         A    C
res_df <- data.frame(
Start = start(matches),
End = end(matches),
Sequence = as.character(matches)
)
print(res_df)
} else {
cat("              7   Motif \n")
}


dna_site_1 <- subseq(seq, start = 2682, end = 2696)
dna_site_2 <- subseq(seq, start = 2385, end = 2399)

```

## Source PDF page 115

```text
#
cat("\n===  AlphaFold 3         ===\n")
cat("  A (AGCCGCC )    :\n", as.character(dna_site_1), "\n\n")
cat("  B (CGCCGCC )    :\n", as.character(dna_site_2), "\n")
```

## Source PDF page 116

```text
tsv_file <- "output/reg.tsv"
rssMat <- readRDS("output/rssMat.rds")
target_gene <- "TraesCS1B02G276800"

message(">>>           output/reg.tsv ...")
lines <- readLines(tsv_file)
valid_lines <- lines[grep("\\[\\(", lines)]

if (length(valid_lines) == 0) {
stop("         reg.tsv")
}
message(paste0("         ", length(valid_lines), "  Regulon"))
message(paste0(">>>          [ ", target_gene, " ]      ..."))

hits <- grep(target_gene, valid_lines, value = TRUE)
if (length(hits) == 0) {
stop(paste0("   SCENIC             ", target_gene, " \nTF                          "))
}
found_tfs <- sapply(strsplit(hits, "\t"), function(x) x[1])
found_tfs <- unique(found_tfs)
message(paste0("    ", length(found_tfs), "   TF: ", paste(found_tfs, collapse=", ")))
cluster_col <- "1"
if (!cluster_col %in% colnames(rssMat)) cluster_col <- grep("1", colnames(rssMat), value = TRUE)[1]
valid_regulators <- c()
scores <- c()

for (tf in found_tfs) {
pattern <- paste0("^", tf, "(_extended)?\\(\\+\\)")
match_row <- grep(pattern, rownames(rssMat), value = TRUE)
if (length(match_row) > 0) {
#             regulon (   TF     regulon)
best_regulon <- match_row[which.max(rssMat[match_row, cluster_col])]
valid_regulators <- c(valid_regulators, best_regulon)
scores <- c(scores, rssMat[best_regulon, cluster_col])
}
}
if (length(valid_regulators) > 0) {
```

## Source PDF page 117

```text
result_df <- data.frame(
Regulon = valid_regulators,
RSS_Cluster8 = scores
)
result_df <- result_df[order(result_df$RSS_Cluster8, decreasing = TRUE), ]
message("\n ---          8            TF    ---")
print(result_df)
top_tf <- result_df$Regulon[1]
message(paste0("\n            [ ", top_tf, " ]"))
message("    8                                 ")
} else {
message("      TF       RSS                               ")
print(found_tfs)
}
```

## Source PDF page 118

```text
library(ComplexHeatmap)
library(circlize)

# 1.
plot_mat <- rssMat[selected_regulons, , drop = FALSE]

#      Cluster                 1, 10, 2
cluster_names <- colnames(plot_mat)
#
sorted_clusters <- cluster_names[order(as.numeric(gsub("[^0-9]", "", cluster_names)))]
plot_mat <- plot_mat[, sorted_clusters, drop = FALSE]

# 2.
min_val <- min(plot_mat)
max_val <- max(plot_mat)
mid_val <- (min_val + max_val) / 2

col_fun = colorRamp2(
c(min_val, mid_val, max_val),
c("#053061", "#F7F7F7", "#B2182B")
)

# 3.
pdf("RSS_Heatmap_Sorted_Clusters.pdf", width = 12, height = 7)

ht <- Heatmap(plot_mat,
name = "RSS",
col = col_fun,

#
cluster_rows = TRUE,      #      TF             TF
cluster_columns = FALSE,  #                   1, 2,3...

#
rect_gp = gpar(col = "white", lwd = 0.8),
row_names_gp = gpar(fontsize = 9, fontface = "italic"),
column_names_gp = gpar(fontsize = 11, fontface = "bold"),

#
heatmap_legend_param = list(
```

## Source PDF page 119

```text
title = "RSS Score",
at = round(seq(min_val, max_val, length.out = 5), 3),
legend_height = unit(4, "cm")
)
)

draw(ht)
dev.off()
```

## Source PDF page 120

```text
library(dplyr)
library(tibble)
library(data.table)
library(tidyr)
library(stringr)

rss_mat <- readRDS("output/rssMat.rds")
auc_mtx <- read.csv("output/auc_mtx.csv", row.names = 1, check.names =FALSE)

reg_raw <- fread("output/reg.tsv", sep = "\t", header = TRUE)

col_tf <- grep("TF", colnames(reg_raw), ignore.case = TRUE, value = TRUE)
col_target <- grep("Target", colnames(reg_raw), ignore.case = TRUE, value = TRUE)

if (length(col_tf) == 0 || length(col_target) == 0) {
tf_vec <- reg_raw[[1]]
target_col_idx <- which(sapply(reg_raw, function(x) any(grepl("\\('", x[1:min(10, length(x))]))))
if(length(target_col_idx) > 0) {
target_vec <- reg_raw[[target_col_idx[1]]]
} else {
target_vec <- reg_raw[[ncol(reg_raw)]]
}
} else {
tf_vec <- reg_raw[[col_tf[1]]]
target_vec <- reg_raw[[col_target[1]]]
}

reg_df <- data.frame(TF_Raw = tf_vec, Targets_Raw = target_vec)

reg_df$TF_Clean <- str_extract(reg_df$TF_Raw, "TraesCS[0-9A-Z]+G[0-9]+")
reg_df$TF_Clean <- ifelse(is.na(reg_df$TF_Clean), as.character(reg_df$TF_Raw), reg_df$TF_Clean)

all_targets_list <- str_extract_all(reg_df$Targets_Raw, "TraesCS[0-9A-Z]+G[0-9]+")

target_summary <- data.frame(TF = reg_df$TF_Clean) %>%
mutate(Target_List = all_targets_list) %>%
```

## Source PDF page 121

```text
unnest(Target_List) %>%
distinct(TF, Target_List) %>%
group_by(TF) %>%
summarise(
Target_Count = n(),
Top_5_Targets = paste(head(Target_List, 5), collapse = ", "),
All_Targets = paste(Target_List, collapse = ";")
) %>%
rename(Clean_TF = TF)

rss_df <- as.data.frame(rss_mat) %>%
rownames_to_column("Regulon") %>%
pivot_longer(cols = -Regulon, names_to = "Cluster", values_to = "RSS") %>%
group_by(Regulon) %>%
filter(RSS == max(RSS)) %>%
ungroup() %>%
rename(Best_Cluster = Cluster, Max_RSS = RSS) %>%
mutate(Max_RSS = round(Max_RSS, 3)) %>%
mutate(Clean_TF = str_extract(Regulon, "TraesCS[0-9A-Z]+G[0-9]+")) %>%
mutate(Clean_TF = ifelse(is.na(Clean_TF), Regulon, Clean_TF))

final_df <- rss_df %>%
left_join(target_summary, by = "Clean_TF") %>%
left_join(
data.frame(Regulon = names(colMeans(auc_mtx)),
Mean_AUC_Activity = round(colMeans(auc_mtx), 4)),
by = "Regulon"
) %>%
mutate(TF_ID = Clean_TF) %>%
select(
Regulon, TF_ID, Best_Cluster, Max_RSS, Mean_AUC_Activity,
Target_Count, Top_5_Targets, All_Targets
) %>%
arrange(desc(Max_RSS))

write.csv(final_df, "SCENIC_Summary.csv", row.names = FALSE, quote = TRUE)
print(head(final_df[, 1:6]))
```
