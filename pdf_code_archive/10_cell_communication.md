# Plant ligand-receptor and CellChat

PDF transcription only. The PDF may wrap code or omit glyphs;
use the curated scripts for executable analysis.

## Source PDF page 142

```text
cor_thresh = -0.05,
ncores=8
)

#         :
pos_regulon_scores <- GetRegulonScores(seurat_obj, target_type='positive')
neg_regulon_scores <- GetRegulonScores(seurat_obj, target_type='negative')

#         TF
cur_tf <- 'FOXP3'

#  regulon      Seurat metadata
seurat_obj$pos_regulon_score <- pos_regulon_scores[,cur_tf]
seurat_obj$neg_regulon_score <- neg_regulon_scores[,cur_tf]

# plot using FeaturePlot
p1 <- FeaturePlot(seurat_obj, feature=cur_tf) + umap_theme()
p2 <- FeaturePlot(seurat_obj, feature='pos_regulon_score', cols=c('lightgrey', 'red')) + umap_theme()
p3 <- FeaturePlot(seurat_obj, feature='neg_regulon_score', cols=c('lightgrey', 'seagreen')) + umap_theme()

p1 | p2 | p3

# select TF of interest
cur_tf <- 'FOXP3'

# plot with default settings
p <- TFNetworkPlot(seurat_obj, selected_tfs=cur_tf)
p

# plot the FOXP3 network with primary, secondary, and tertiary targets
p1 <- TFNetworkPlot(seurat_obj, selected_tfs=cur_tf, depth=1, no_labels=TRUE)
p2 <- TFNetworkPlot(seurat_obj, selected_tfs=cur_tf, depth=2, no_labels=TRUE)
p3 <- TFNetworkPlot(seurat_obj, selected_tfs=cur_tf, depth=3, no_labels=TRUE)

p1 | p2 | p3
```

## Source PDF page 144

```text
make_lr <- function(ligand, receptor, pathway) {
expand.grid(
ligand = ligand,
receptor = receptor,
pathway_name = pathway,
stringsAsFactors = FALSE
)
}

lr_pairs <- do.call(rbind, list(

make_lr(
"TraesCS1D02G057600",
"TraesCS2B02G251700",
"CIF-SGN"
),

make_lr(
c(
"TraesCS1B02G309600",
"TraesCS2D02G256100",
"TraesCS3A02G346000",
"TraesCS3D02G415500",
"TraesCS4A02G292400",
"TraesCS4B02G277700",
"TraesCS6B02G347500"
),
"TraesCS7D02G166100",
"EPF_EPFL-ER"
),

make_lr(
c(
"TraesCS3A02G176000",
"TraesCS3D02G181700",
"TraesCS5A02G093000",
"TraesCS5B02G099100",
"TraesCS5D02G105300"
),
c(
"TraesCS4A02G441600",
"TraesCS4B02G002300"
),
```

## Source PDF page 145

```text
"RALF-BUPS_ANX"
),

make_lr(
c(
"TraesCS5A02G093000",
"TraesCS5B02G099100",
"TraesCS5D02G105300"
),
c(
"TraesCS1D02G228900",
"TraesCS4A02G133800"
),
"RALF-FER"
),

make_lr(
c(
"TraesCS5A02G093000",
"TraesCS5B02G099100",
"TraesCS5D02G105300"
),
"TraesCS7B02G241900",
"RALF-THE"
),

make_lr(
c(
"TraesCS2A02G276500",
"TraesCS2A02G276600",
"TraesCS2B02G294200",
"TraesCS2B02G294300",
"TraesCS2D02G275500",
"TraesCS2D02G275600"
),
"TraesCS1A02G342600",
"DCC-DCCR"
),

make_lr(
c(
"TraesCS4A02G024100",
"TraesCS4B02G279400"
),
"TraesCS3B02G310000",
"DEP-DEPR"
```

## Source PDF page 146

```text
),

make_lr(
c(
"TraesCS4A02G024100",
"TraesCS4A02G024200",
"TraesCS4A02G024300",
"TraesCS4A02G024400",
"TraesCS4A02G025000",
"TraesCS4B02G279200",
"TraesCS4B02G279300",
"TraesCS4B02G279400",
"TraesCS4D02G277600",
"TraesCS4D02G277700",
"TraesCS4D02G277800",
"TraesCS4D02G277900"
),
"TraesCS3D02G276200",
"TaFIP-TaFIPR"
)
))

#    merged SCT
features <- rownames(merged[["SCT"]])

lr_pairs <- lr_pairs[
lr_pairs$ligand %in% features &
lr_pairs$receptor %in% features,
,
drop = FALSE
]

#
lr_pairs <- unique(lr_pairs)

#   CellChat
lr_pairs$interaction_name <- paste(
lr_pairs$ligand,
lr_pairs$receptor,
sep = "_"
)

lr_pairs$agonist <- ""
lr_pairs$antagonist <- ""
lr_pairs$co_A_receptor <- ""
lr_pairs$co_I_receptor <- ""
```

## Source PDF page 147

```text
lr_pairs$annotation <- "Secreted Signaling"

lr_pairs$interaction_name_2 <- paste(
lr_pairs$ligand,
lr_pairs$receptor,
sep = " - "
)

# CellChat
lr_pairs <- lr_pairs[, c(
"interaction_name",
"pathway_name",
"ligand",
"receptor",
"agonist",
"antagonist",
"co_A_receptor",
"co_I_receptor",
"annotation",
"interaction_name_2"
)]

rownames(lr_pairs) <- lr_pairs$interaction_name

cat("lr_pairs     ", nrow(lr_pairs), " \n")
print(lr_pairs)
```

## Source PDF page 148

```text
############################################################
## Wheat single-cell CellChat analysis
## CK vs FCR
##            seed = 54
##
##
## 1. Seurat  merged
## 2. CellChat           lr_pairs
############################################################

############################################################
## 0.  R
############################################################

suppressPackageStartupMessages({
library(Seurat)
library(CellChat)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(future)
})

options(stringsAsFactors = FALSE)
set.seed(54)

############################################################
## 1.
############################################################

sample_col <- "Samples"
celltype_col <- "celltype"
assay_use  <- "SCT"

condition_order <- c("ck", "fcr")

##
min_cells_per_group <- 10

##
min_expression_pct <- 0.10
```

## Source PDF page 149

```text

## truncated mean
trim_use <- 0.10

## CellChat
nboot_use <- 500

## CellChat
p_cutoff <- 0.05

## BH
fdr_cutoff <- 0.05

##                   CK FCR
balance_cell_numbers <- TRUE

##
workers_use <- 4

##
seed_use <- 54

##
outdir <- "cellchat_wheat_seed54"

dir.create(
outdir,
recursive = TRUE,
showWarnings = FALSE
)

############################################################
## 2.  merged lr_pairs
############################################################

if (!inherits(merged, "Seurat")) {
stop("merged Seurat    ")
}

if (!is.data.frame(lr_pairs)) {
stop("lr_pairs   data.frame ")
}

if (!sample_col %in% colnames(merged@meta.data)) {
stop(
"merged@meta.data         ",
```

## Source PDF page 150

```text
sample_col
)
}

if (!celltype_col %in% colnames(merged@meta.data)) {
stop(
"merged@meta.data            ",
celltype_col
)
}

if (!assay_use %in% Assays(merged)) {
stop(
"    assay ", assay_use,
"\n  assay ",
paste(Assays(merged), collapse = ", ")
)
}

required_lr_columns <- c(
"interaction_name",
"pathway_name",
"ligand",
"receptor"
)

missing_lr_columns <- setdiff(
required_lr_columns,
colnames(lr_pairs)
)

if (length(missing_lr_columns) > 0) {
stop(
"lr_pairs          ",
paste(missing_lr_columns, collapse = ", ")
)
}

############################################################
## 3.    merged
##            Samples celltype
############################################################

merged$.condition <- tolower(
trimws(
as.character(
```

## Source PDF page 151

```text
merged@meta.data[[sample_col]]
)
)
)

merged$.celltype <- trimws(
as.character(
merged@meta.data[[celltype_col]]
)
)

valid_metadata <- (
!is.na(merged$.condition) &
!is.na(merged$.celltype) &
merged$.condition != "" &
merged$.celltype != "" &
merged$.condition %in% condition_order
)

if (!all(valid_metadata)) {
message(
"            ",
sum(!valid_metadata),
"                        "
)
}

analysis_cells <- colnames(merged)[valid_metadata]

merged_analysis <- subset(
merged,
cells = analysis_cells
)

missing_conditions <- setdiff(
condition_order,
unique(merged_analysis$.condition)
)

if (length(missing_conditions) > 0) {
stop(
"             ",
paste(missing_conditions, collapse = ", ")
)
}

```

## Source PDF page 152

```text
############################################################
## 4.  lr_pairs
##                 lr_pairs
############################################################

lr_pairs_analysis <- as.data.frame(
lr_pairs,
stringsAsFactors = FALSE
)

for (column_name in required_lr_columns) {
lr_pairs_analysis[[column_name]] <- trimws(
as.character(
lr_pairs_analysis[[column_name]]
)
)
}

valid_lr_rows <- (
!is.na(lr_pairs_analysis$interaction_name) &
!is.na(lr_pairs_analysis$pathway_name) &
!is.na(lr_pairs_analysis$ligand) &
!is.na(lr_pairs_analysis$receptor) &
lr_pairs_analysis$interaction_name != "" &
lr_pairs_analysis$pathway_name != "" &
lr_pairs_analysis$ligand != "" &
lr_pairs_analysis$receptor != ""
)

if (!all(valid_lr_rows)) {
warning(
"   ",
sum(!valid_lr_rows),
"                      "
)

lr_pairs_analysis <- lr_pairs_analysis[
valid_lr_rows,
,
drop = FALSE
]
}

##
lr_pairs_analysis <- lr_pairs_analysis %>%
distinct(
```

## Source PDF page 153

```text
interaction_name,
ligand,
receptor,
pathway_name,
.keep_all = TRUE
)

## interaction_name
if (anyDuplicated(lr_pairs_analysis$interaction_name)) {
duplicate_names <- unique(
lr_pairs_analysis$interaction_name[
duplicated(
lr_pairs_analysis$interaction_name
)
]
)

stop(
"   interaction_name   ",
paste(duplicate_names, collapse = ", ")
)
}

##   CellChat
optional_columns <- c(
"agonist",
"antagonist",
"co_A_receptor",
"co_I_receptor"
)

for (column_name in optional_columns) {
if (!column_name %in% colnames(lr_pairs_analysis)) {
lr_pairs_analysis[[column_name]] <- ""
}

lr_pairs_analysis[[column_name]][
is.na(lr_pairs_analysis[[column_name]])
] <- ""
}

if (!"annotation" %in% colnames(lr_pairs_analysis)) {
lr_pairs_analysis$annotation <- "Secreted Signaling"
}

lr_pairs_analysis$annotation[
```

## Source PDF page 154

```text
is.na(lr_pairs_analysis$annotation) |
lr_pairs_analysis$annotation == ""
] <- "Secreted Signaling"

if (!"interaction_name_2" %in% colnames(lr_pairs_analysis)) {
lr_pairs_analysis$interaction_name_2 <- paste(
lr_pairs_analysis$ligand,
lr_pairs_analysis$receptor,
sep = " - "
)
}

rownames(lr_pairs_analysis) <-
lr_pairs_analysis$interaction_name

message(
"             ",
nrow(lr_pairs_analysis)
)

message(
"       ",
paste(
unique(lr_pairs_analysis$pathway_name),
collapse = ", "
)
)

############################################################
## 5.
############################################################

assay_features <- rownames(
merged_analysis[[assay_use]]
)

lr_gene_check <- lr_pairs_analysis %>%
transmute(
interaction_name,
pathway_name,
ligand,
receptor,
ligand_found = ligand %in% assay_features,
receptor_found = receptor %in% assay_features,
pair_available = (
ligand %in% assay_features &
```

## Source PDF page 155

```text
receptor %in% assay_features
)
)

write.table(
lr_gene_check,
file = file.path(
outdir,
"LR_gene_availability.tsv"
),
sep = "\t",
quote = FALSE,
row.names = FALSE
)

if (any(!lr_gene_check$pair_available)) {
warning(
sum(!lr_gene_check$pair_available),
"                      ",
assay_use,
" assay           "
)
}

lr_pairs_available <- lr_pairs_analysis[
lr_gene_check$pair_available,
,
drop = FALSE
]

if (nrow(lr_pairs_available) == 0) {
stop("                        ")
}

rownames(lr_pairs_available) <-
lr_pairs_available$interaction_name

all_lr_genes <- unique(
c(
lr_pairs_available$ligand,
lr_pairs_available$receptor
)
)

message(
"      CellChat          ",
```

## Source PDF page 156

```text
nrow(lr_pairs_available)
)

message(
"              ",
length(all_lr_genes)
)

############################################################
## 6.  CK  FCR
############################################################

cell_count_before <- as.data.frame(
table(
celltype = merged_analysis$.celltype,
condition = merged_analysis$.condition
),
stringsAsFactors = FALSE
)

colnames(cell_count_before)[3] <- "n_cells"

cell_count_matrix <- table(
merged_analysis$.celltype,
factor(
merged_analysis$.condition,
levels = condition_order
)
)

common_celltypes <- rownames(cell_count_matrix)[
apply(
cell_count_matrix[
,
condition_order,
drop = FALSE
] >= min_cells_per_group,
1,
all
)
]

if (length(common_celltypes) == 0) {
stop(
"             CK FCR     ",
min_cells_per_group,
```

## Source PDF page 157

```text
"      "
)
}

common_celltypes <- sort(common_celltypes)

message(
"            ",
paste(common_celltypes, collapse = ", ")
)

common_cells <- colnames(merged_analysis)[
merged_analysis$.celltype %in% common_celltypes &
merged_analysis$.condition %in% condition_order
]

merged_common <- subset(
merged_analysis,
cells = common_cells
)

############################################################
## 7.          CK FCR
############################################################

if (balance_cell_numbers) {

selected_balanced_cells <- character(0)

for (
celltype_index in seq_along(common_celltypes)
) {

current_celltype <-
common_celltypes[celltype_index]

ck_cells <- colnames(merged_common)[
merged_common$.celltype == current_celltype &
merged_common$.condition == "ck"
]

fcr_cells <- colnames(merged_common)[
merged_common$.celltype == current_celltype &
merged_common$.condition == "fcr"
]

```

## Source PDF page 158

```text
target_cell_number <- min(
length(ck_cells),
length(fcr_cells)
)

set.seed(
seed_use +
celltype_index * 100 +
1
)

selected_ck_cells <- sample(
ck_cells,
size = target_cell_number,
replace = FALSE
)

set.seed(
seed_use +
celltype_index * 100 +
2
)

selected_fcr_cells <- sample(
fcr_cells,
size = target_cell_number,
replace = FALSE
)

selected_balanced_cells <- c(
selected_balanced_cells,
selected_ck_cells,
selected_fcr_cells
)
}

merged_balanced <- subset(
merged_common,
cells = selected_balanced_cells
)

} else {

merged_balanced <- merged_common
}

```

## Source PDF page 159

```text
merged_balanced$.celltype <- factor(
as.character(merged_balanced$.celltype),
levels = common_celltypes
)

cell_count_after <- as.data.frame(
table(
celltype = merged_balanced$.celltype,
condition = factor(
merged_balanced$.condition,
levels = condition_order
)
),
stringsAsFactors = FALSE
)

colnames(cell_count_after)[3] <- "n_cells"

message("           ")
print(cell_count_after)

############################################################
## 8.
##   Seurat v4 v5
############################################################

get_expression_data <- function(
seu,
assay_name
) {

expression_matrix <- tryCatch(
GetAssayData(
seu,
assay = assay_name,
layer = "data"
),
error = function(e) NULL
)

if (is.null(expression_matrix)) {
expression_matrix <- tryCatch(
GetAssayData(
seu,
assay = assay_name,
slot = "data"
```

## Source PDF page 160

```text
),
error = function(e) NULL
)
}

if (
is.null(expression_matrix) ||
nrow(expression_matrix) == 0 ||
ncol(expression_matrix) == 0
) {
stop(
"    assay ",
assay_name,
"                "
)
}

expression_matrix
}

############################################################
## 9.          CellChat
############################################################

gene_info_wheat <- data.frame(
Symbol = all_lr_genes,
row.names = all_lr_genes,
stringsAsFactors = FALSE
)

CellChatDB_wheat <- updateCellChatDB(
db = lr_pairs_available,
gene_info = gene_info_wheat,
merged = FALSE
)

if (
is.null(CellChatDB_wheat$interaction) ||
nrow(CellChatDB_wheat$interaction) == 0
) {
stop("       CellChat          ")
}

message(
"                      ",
nrow(CellChatDB_wheat$interaction)
```

## Source PDF page 161

```text
)

############################################################
## 10.  CellChat
############################################################

run_cellchat <- function(
seu,
condition_name
) {

message(
"\n========== Processing ",
condition_name,
" =========="
)

set.seed(seed_use)

cells_use <- colnames(seu)[
as.character(seu$.condition) == condition_name
]

if (length(cells_use) == 0) {
stop(
"           ",
condition_name
)
}

seu_sub <- subset(
seu,
cells = cells_use
)

labels <- factor(
as.character(seu_sub$.celltype),
levels = common_celltypes
)

labels <- droplevels(labels)
names(labels) <- colnames(seu_sub)

data_input <- get_expression_data(
seu_sub,
assay_name = assay_use
```

## Source PDF page 162

```text
)

genes_use <- intersect(
all_lr_genes,
rownames(data_input)
)

data_input <- data_input[
genes_use,
colnames(seu_sub),
drop = FALSE
]

if (nrow(data_input) == 0) {
stop(
condition_name,
"                    "
)
}

meta_input <- data.frame(
labels = labels,
samples = rep(
condition_name,
length(labels)
),
row.names = colnames(seu_sub),
stringsAsFactors = FALSE
)

message("Cells: ", ncol(data_input))
message("Genes: ", nrow(data_input))
message(
"Cell types: ",
length(unique(labels))
)

cellchat <- createCellChat(
object = data_input,
meta = meta_input,
group.by = "labels"
)

cellchat@DB <- CellChatDB_wheat

cellchat <- subsetData(cellchat)
```

## Source PDF page 163

```text

if (nrow(cellchat@data.signaling) == 0) {
stop(
condition_name,
" subsetData()          "
)
}

##########################################################
##                    LR
##
##     identifyOverExpressedInteractions()
##                     1 LR
##
##
## 1.
## 2. CellChat P
## 3. BH
##
##########################################################

pairLR_use <- cellchat@DB$interaction

cellchat@LR$LRsig <- pairLR_use

message(
"Candidate LR pairs tested: ",
nrow(pairLR_use)
)

set.seed(seed_use)

cellchat <- computeCommunProb(
object = cellchat,
type = "truncatedMean",
trim = trim_use,
LR.use = pairLR_use,
raw.use = TRUE,
population.size = FALSE,
distance.use = FALSE,
nboot = nboot_use,
seed.use = seed_use
)

cellchat <- filterCommunication(
cellchat,
```

## Source PDF page 164

```text
min.cells = min_cells_per_group
)

if (length(dim(cellchat@net$prob)) != 3) {
stop(
"cellchat@net$prob       ",
"    CellChat            "
)
}

message(
"LR pairs in probability array: ",
dim(cellchat@net$prob)[3]
)

message(
"Non-zero source-target-LR combinations: ",
sum(
cellchat@net$prob > 0,
na.rm = TRUE
)
)

##########################################################
##
##########################################################

group_levels <- levels(labels)

pct_matrix <- do.call(
cbind,
lapply(
group_levels,
function(current_group) {

group_cells <- names(labels)[
labels == current_group
]

Matrix::rowMeans(
data_input[
,
group_cells,
drop = FALSE
] > 0
)
```

## Source PDF page 165

```text
}
)
)

mean_matrix <- do.call(
cbind,
lapply(
group_levels,
function(current_group) {

group_cells <- names(labels)[
labels == current_group
]

Matrix::rowMeans(
data_input[
,
group_cells,
drop = FALSE
]
)
}
)
)

rownames(pct_matrix) <- rownames(data_input)
rownames(mean_matrix) <- rownames(data_input)

colnames(pct_matrix) <- group_levels
colnames(mean_matrix) <- group_levels

list(
cellchat = cellchat,
seurat = seu_sub,
labels = labels,
data_input = data_input,
pct_matrix = pct_matrix,
mean_matrix = mean_matrix,
condition = condition_name
)
}

############################################################
## 11.     CK FCR
############################################################

```

## Source PDF page 166

```text
set.seed(seed_use)

future::plan(
future::multisession,
workers = workers_use
)

result_ck <- run_cellchat(
merged_balanced,
"ck"
)

result_fcr <- run_cellchat(
merged_balanced,
"fcr"
)

future::plan(
future::sequential
)

cellchat_ck <- result_ck$cellchat
cellchat_fcr <- result_fcr$cellchat

############################################################
## 12.
############################################################

matrix_lookup <- function(
input_matrix,
genes,
groups
) {

row_index <- match(
genes,
rownames(input_matrix)
)

column_index <- match(
groups,
colnames(input_matrix)
)

output <- rep(
NA_real_,
```

## Source PDF page 167

```text
length(genes)
)

valid_index <- (
!is.na(row_index) &
!is.na(column_index)
)

output[valid_index] <- input_matrix[
cbind(
row_index[valid_index],
column_index[valid_index]
)
]

output
}

############################################################
## 13.  LR    CellChat
############################################################

extract_communication <- function(
result_object
) {

cellchat <- result_object$cellchat

prob_array <- cellchat@net$prob
pval_array <- cellchat@net$pval

array_names <- dimnames(prob_array)

output <- expand.grid(
source = array_names[[1]],
target = array_names[[2]],
interaction_name = array_names[[3]],
KEEP.OUT.ATTRS = FALSE,
stringsAsFactors = FALSE
)

output$prob <- as.vector(prob_array)
output$pval <- as.vector(pval_array)

database_interactions <-
cellchat@DB$interaction
```

## Source PDF page 168

```text

annotation_index <- match(
output$interaction_name,
database_interactions$interaction_name
)

if (anyNA(annotation_index)) {
stop(
"  interaction_name              "
)
}

output$ligand <-
database_interactions$ligand[
annotation_index
]

output$receptor <-
database_interactions$receptor[
annotation_index
]

output$pathway_name <-
database_interactions$pathway_name[
annotation_index
]

output$ligand_pct <- matrix_lookup(
result_object$pct_matrix,
output$ligand,
output$source
)

output$receptor_pct <- matrix_lookup(
result_object$pct_matrix,
output$receptor,
output$target
)

output$ligand_mean <- matrix_lookup(
result_object$mean_matrix,
output$ligand,
output$source
)

output$receptor_mean <- matrix_lookup(
```

## Source PDF page 169

```text
result_object$mean_matrix,
output$receptor,
output$target
)

output$pval[
is.na(output$pval)
] <- 1

##            source-target-LR    BH
output$permutation_fdr <- p.adjust(
output$pval,
method = "BH"
)

output$pass_expression <- (
!is.na(output$ligand_pct) &
!is.na(output$receptor_pct) &
output$ligand_pct >= min_expression_pct &
output$receptor_pct >= min_expression_pct
)

## retained
##     CK/FCR
output$retained <- (
output$prob > 0 &
output$pass_expression &
output$pval <= p_cutoff &
output$permutation_fdr <= fdr_cutoff
)

output$condition <-
result_object$condition

output %>%
select(
condition,
source,
target,
interaction_name,
pathway_name,
ligand,
receptor,
prob,
pval,
permutation_fdr,
```

## Source PDF page 170

```text
ligand_pct,
receptor_pct,
ligand_mean,
receptor_mean,
pass_expression,
retained
) %>%
arrange(
desc(retained),
permutation_fdr,
desc(prob)
)
}

communication_ck <- extract_communication(
result_ck
)

communication_fcr <- extract_communication(
result_fcr
)

retained_ck <- communication_ck %>%
filter(retained)

retained_fcr <- communication_fcr %>%
filter(retained)

message(
"CK retained source-target-LR routes: ",
nrow(retained_ck)
)

message(
"FCR retained source-target-LR routes: ",
nrow(retained_fcr)
)

############################################################
## 14. CK FCR
############################################################

comparison_ck <- communication_ck %>%
select(
source,
target,
```

## Source PDF page 171

```text
interaction_name,
pathway_name,
ligand,
receptor,
prob_ck = prob,
pval_ck = pval,
fdr_ck = permutation_fdr,
ligand_pct_ck = ligand_pct,
receptor_pct_ck = receptor_pct,
ligand_mean_ck = ligand_mean,
receptor_mean_ck = receptor_mean,
retained_ck = retained
)

comparison_fcr <- communication_fcr %>%
select(
source,
target,
interaction_name,
pathway_name,
ligand,
receptor,
prob_fcr = prob,
pval_fcr = pval,
fdr_fcr = permutation_fdr,
ligand_pct_fcr = ligand_pct,
receptor_pct_fcr = receptor_pct,
ligand_mean_fcr = ligand_mean,
receptor_mean_fcr = receptor_mean,
retained_fcr = retained
)

communication_comparison <- full_join(
comparison_ck,
comparison_fcr,
by = c(
"source",
"target",
"interaction_name",
"pathway_name",
"ligand",
"receptor"
)
) %>%
mutate(
prob_ck = coalesce(prob_ck, 0),
```

## Source PDF page 172

```text
prob_fcr = coalesce(prob_fcr, 0),

pval_ck = coalesce(pval_ck, 1),
pval_fcr = coalesce(pval_fcr, 1),

fdr_ck = coalesce(fdr_ck, 1),
fdr_fcr = coalesce(fdr_fcr, 1),

retained_ck = coalesce(
retained_ck,
FALSE
),

retained_fcr = coalesce(
retained_fcr,
FALSE
),

retained_prob_ck = ifelse(
retained_ck,
prob_ck,
0
),

retained_prob_fcr = ifelse(
retained_fcr,
prob_fcr,
0
),

delta_prob = (
retained_prob_fcr -
retained_prob_ck
),

relative_delta = delta_prob / (
abs(retained_prob_fcr) +
abs(retained_prob_ck) +
1e-12
),

support_status = case_when(
retained_ck & retained_fcr ~
"Retained in both",

retained_ck & !retained_fcr ~
```

## Source PDF page 173

```text
"CK only",

!retained_ck & retained_fcr ~
"FCR only",

TRUE ~
"Not retained"
)
) %>%
arrange(
desc(abs(delta_prob))
)

############################################################
## 15. Sender-receiver
############################################################

sender_receiver_summary <-
communication_comparison %>%
group_by(
source,
target
) %>%
summarise(
strength_ck = sum(
retained_prob_ck,
na.rm = TRUE
),

strength_fcr = sum(
retained_prob_fcr,
na.rm = TRUE
),

delta_strength = (
strength_fcr -
strength_ck
),

n_lr_ck = sum(
retained_ck,
na.rm = TRUE
),

n_lr_fcr = sum(
retained_fcr,
```

## Source PDF page 174

```text
na.rm = TRUE
),

.groups = "drop"
) %>%
arrange(
desc(abs(delta_strength))
)

############################################################
## 16.     pathway
##    computeCommunProbPathway()
############################################################

pathway_summary <-
communication_comparison %>%
group_by(
pathway_name,
source,
target
) %>%
summarise(
strength_ck = sum(
retained_prob_ck,
na.rm = TRUE
),

strength_fcr = sum(
retained_prob_fcr,
na.rm = TRUE
),

delta_strength = (
strength_fcr -
strength_ck
),

n_lr_ck = sum(
retained_ck,
na.rm = TRUE
),

n_lr_fcr = sum(
retained_fcr,
na.rm = TRUE
),
```

## Source PDF page 175

```text

.groups = "drop"
) %>%
arrange(
pathway_name,
desc(abs(delta_strength))
)

pathway_total <- pathway_summary %>%
group_by(pathway_name) %>%
summarise(
total_strength_ck = sum(
strength_ck,
na.rm = TRUE
),

total_strength_fcr = sum(
strength_fcr,
na.rm = TRUE
),

delta_total_strength = (
total_strength_fcr -
total_strength_ck
),

total_lr_routes_ck = sum(
n_lr_ck,
na.rm = TRUE
),

total_lr_routes_fcr = sum(
n_lr_fcr,
na.rm = TRUE
),

.groups = "drop"
) %>%
arrange(
desc(abs(delta_total_strength))
)

############################################################
## 17.        ×
############################################################

```

## Source PDF page 176

```text
build_network_matrix <- function(
communication_table,
value_column,
cell_order
) {

output_matrix <- matrix(
0,
nrow = length(cell_order),
ncol = length(cell_order),
dimnames = list(
cell_order,
cell_order
)
)

if (nrow(communication_table) == 0) {
return(output_matrix)
}

aggregated_table <- communication_table %>%
group_by(
source,
target
) %>%
summarise(
value = sum(
.data[[value_column]],
na.rm = TRUE
),
.groups = "drop"
)

for (
row_number in seq_len(
nrow(aggregated_table)
)
) {

source_index <- match(
aggregated_table$source[row_number],
cell_order
)

target_index <- match(
aggregated_table$target[row_number],
```

## Source PDF page 177

```text
cell_order
)

if (
!is.na(source_index) &&
!is.na(target_index)
) {
output_matrix[
source_index,
target_index
] <- aggregated_table$value[row_number]
}
}

output_matrix
}

strength_matrix_ck <- build_network_matrix(
retained_ck,
"prob",
common_celltypes
)

strength_matrix_fcr <- build_network_matrix(
retained_fcr,
"prob",
common_celltypes
)

count_table_ck <- retained_ck %>%
mutate(lr_count = 1)

count_table_fcr <- retained_fcr %>%
mutate(lr_count = 1)

count_matrix_ck <- build_network_matrix(
count_table_ck,
"lr_count",
common_celltypes
)

count_matrix_fcr <- build_network_matrix(
count_table_fcr,
"lr_count",
common_celltypes
)
```

## Source PDF page 178

```text

delta_strength_matrix <- (
strength_matrix_fcr -
strength_matrix_ck
)

############################################################
## 18.
############################################################

get_group_sizes <- function(
result_object
) {

as.numeric(
table(
factor(
result_object$labels,
levels = common_celltypes
)
)
)
}

group_size_ck <- get_group_sizes(
result_ck
)

group_size_fcr <- get_group_sizes(
result_fcr
)

############################################################
## 19.
############################################################

plot_circle_safe <- function(
network_matrix,
group_sizes,
output_file,
title_text
) {

pdf(
file.path(
outdir,
```

## Source PDF page 179

```text
output_file
),
width = 12,
height = 12,
useDingbats = FALSE
)

if (
length(network_matrix) > 0 &&
any(
network_matrix > 0,
na.rm = TRUE
)
) {

netVisual_circle(
network_matrix,
vertex.weight = pmax(
group_sizes,
1
),
weight.scale = TRUE,
label.edge = FALSE,
edge.weight.max = max(
network_matrix,
na.rm = TRUE
),
title.name = title_text
)

} else {

plot.new()

title(
paste0(
title_text,
"\nNo retained communication"
)
)
}

dev.off()
}

plot_circle_safe(
```

## Source PDF page 180

```text
strength_matrix_ck,
group_size_ck,
"CK_retained_strength_circle.pdf",
"CK retained communication strength"
)

plot_circle_safe(
strength_matrix_fcr,
group_size_fcr,
"FCR_retained_strength_circle.pdf",
"FCR retained communication strength"
)

plot_circle_safe(
count_matrix_ck,
group_size_ck,
"CK_retained_LR_count_circle.pdf",
"CK retained LR-pair count"
)

plot_circle_safe(
count_matrix_fcr,
group_size_fcr,
"FCR_retained_LR_count_circle.pdf",
"FCR retained LR-pair count"
)

############################################################
## 20. CK-FCR
############################################################

heatmap_data <- expand.grid(
source = rownames(
delta_strength_matrix
),
target = colnames(
delta_strength_matrix
),
KEEP.OUT.ATTRS = FALSE,
stringsAsFactors = FALSE
)

heatmap_data$delta_strength <-
as.vector(delta_strength_matrix)

max_abs_delta <- max(
```

## Source PDF page 181

```text
abs(heatmap_data$delta_strength),
na.rm = TRUE
)

if (
!is.finite(max_abs_delta) ||
max_abs_delta == 0
) {
max_abs_delta <- 1
}

p_delta_heatmap <- ggplot(
heatmap_data,
aes(
x = target,
y = source,
fill = delta_strength
)
) +
geom_tile(
color = "grey90",
linewidth = 0.25
) +
scale_fill_gradient2(
low = "#2166AC",
mid = "white",
high = "#B2182B",
midpoint = 0,
limits = c(
-max_abs_delta,
max_abs_delta
),
name = "FCR - CK"
) +
labs(
title = paste0(
"Descriptive change in retained communication\n",
"(no biological-replicate inference)"
),
x = "Receiver",
y = "Sender"
) +
theme_bw(base_size = 11) +
theme(
axis.text.x = element_text(
angle = 45,
```

## Source PDF page 182

```text
hjust = 1,
color = "black"
),
axis.text.y = element_text(
color = "black"
),
panel.grid = element_blank(),
plot.title = element_text(
hjust = 0.5,
face = "bold"
)
)

ggsave(
filename = file.path(
outdir,
"CK_FCR_descriptive_delta_heatmap.pdf"
),
plot = p_delta_heatmap,
width = 12,
height = 10
)

############################################################
## 21. Pathway
############################################################

pathway_plot_data <- pathway_total %>%
select(
pathway_name,
CK = total_strength_ck,
FCR = total_strength_fcr
) %>%
pivot_longer(
cols = c(CK, FCR),
names_to = "condition",
values_to = "strength"
)

p_pathway <- ggplot(
pathway_plot_data,
aes(
x = reorder(
pathway_name,
strength,
FUN = max
```

## Source PDF page 183

```text
),
y = strength,
fill = condition
)
) +
geom_col(
position = position_dodge(
width = 0.8
),
width = 0.7
) +
coord_flip() +
scale_fill_manual(
values = c(
"CK" = "#4477AA",
"FCR" = "#CC6677"
)
) +
labs(
title = "Retained pathway-level communication strength",
subtitle = "Descriptive CK–FCR comparison",
x = "Curated LR pathway",
y = "Summed CellChat probability",
fill = NULL
) +
theme_bw(base_size = 11) +
theme(
panel.grid.major.y = element_blank(),
plot.title = element_text(
face = "bold"
)
)

ggsave(
filename = file.path(
outdir,
"pathway_strength_CK_FCR.pdf"
),
plot = p_pathway,
width = 9,
height = 6
)

############################################################
## 22. LR
############################################################
```

## Source PDF page 184

```text

bubble_data <- communication_comparison %>%
filter(
retained_ck |
retained_fcr
) %>%
mutate(
route = paste(
source,
target,→
sep = "  "
)
) %>%
select(
interaction_name,
pathway_name,
route,
prob_ck,
prob_fcr,
retained_ck,
retained_fcr
) %>%
pivot_longer(
cols = c(
prob_ck,
prob_fcr
),
names_to = "condition",
values_to = "prob"
) %>%
mutate(
condition = recode(
condition,
prob_ck = "CK",
prob_fcr = "FCR"
),

retained = ifelse(
condition == "CK",
retained_ck,
retained_fcr
),

prob = ifelse(
retained,
prob,
```

## Source PDF page 185

```text
0
)
) %>%
filter(prob > 0)

if (nrow(bubble_data) > 0) {

top_routes <- bubble_data %>%
group_by(
interaction_name,
route
) %>%
summarise(
max_prob = max(
prob,
na.rm = TRUE
),
.groups = "drop"
) %>%
arrange(desc(max_prob)) %>%
slice_head(n = 40)

bubble_plot_data <- bubble_data %>%
semi_join(
top_routes,
by = c(
"interaction_name",
"route"
)
) %>%
mutate(
display_label = paste(
interaction_name,
route,
sep = " | "
)
)

p_bubble <- ggplot(
bubble_plot_data,
aes(
x = condition,
y = reorder(
display_label,
prob
),
```

## Source PDF page 186

```text
size = prob,
color = pathway_name
)
) +
geom_point(alpha = 0.85) +
scale_size_continuous(
range = c(2, 9)
) +
labs(
title = "Top retained ligand–receptor routes",
subtitle = paste0(
"Internal CellChat support; ",
"not replicate-level treatment inference"
),
x = NULL,→
y = "LR pair | sender receiver",
size = "Probability",
color = "Pathway"
) +
theme_bw(base_size = 10) +
theme(
panel.grid.major.y = element_line(
color = "grey90"
),
legend.position = "right",
plot.title = element_text(
face = "bold"
)
)

ggsave(
filename = file.path(
outdir,
"top_retained_LR_routes.pdf"
),
plot = p_bubble,
width = 13,
height = 12
)
}

############################################################
## 23.
############################################################

overview_table <- data.frame(
```

## Source PDF page 187

```text
metric = c(
"Cells",
"Cell types",
"Candidate LR pairs tested",
"Retained source-target-LR routes",
"Total retained communication strength"
),

CK = c(
ncol(result_ck$seurat),
length(unique(result_ck$labels)),
dim(cellchat_ck@net$prob)[3],
nrow(retained_ck),
sum(
retained_ck$prob,
na.rm = TRUE
)
),

FCR = c(
ncol(result_fcr$seurat),
length(unique(result_fcr$labels)),
dim(cellchat_fcr@net$prob)[3],
nrow(retained_fcr),
sum(
retained_fcr$prob,
na.rm = TRUE
)
),

stringsAsFactors = FALSE
)

print(overview_table)

############################################################
## 24.  Excel
############################################################

output_workbook <- list(
analysis_overview = overview_table,

cell_counts_before =
cell_count_before,

cell_counts_after =
```

## Source PDF page 188

```text
cell_count_after,

input_LR_pairs =
lr_pairs_available,

LR_gene_check =
lr_gene_check,

CK_all_tests =
communication_ck,

FCR_all_tests =
communication_fcr,

CK_retained =
retained_ck,

FCR_retained =
retained_fcr,

CK_FCR_comparison =
communication_comparison,

sender_receiver =
sender_receiver_summary,

pathway_routes =
pathway_summary,

pathway_total =
pathway_total,

CK_strength_matrix =
as.data.frame(
strength_matrix_ck
),

FCR_strength_matrix =
as.data.frame(
strength_matrix_fcr
),

delta_matrix =
as.data.frame(
delta_strength_matrix
)
```

## Source PDF page 189

```text
)

openxlsx::write.xlsx(
output_workbook,
file = file.path(
outdir,
"wheat_CellChat_complete_results.xlsx"
),
overwrite = TRUE,
rowNames = TRUE
)

############################################################
## 25.  CellChat
############################################################

saveRDS(
cellchat_ck,
file = file.path(
outdir,
"cellchat_CK_seed54.rds"
)
)

saveRDS(
cellchat_fcr,
file = file.path(
outdir,
"cellchat_FCR_seed54.rds"
)
)

complete_cellchat_results <- list(
result_ck = result_ck,
result_fcr = result_fcr,

communication_ck =
communication_ck,

communication_fcr =
communication_fcr,

retained_ck =
retained_ck,

retained_fcr =
```

## Source PDF page 190

```text
retained_fcr,

comparison =
communication_comparison,

sender_receiver_summary =
sender_receiver_summary,

pathway_summary =
pathway_summary,

pathway_total =
pathway_total,

strength_matrix_ck =
strength_matrix_ck,

strength_matrix_fcr =
strength_matrix_fcr,

delta_strength_matrix =
delta_strength_matrix,

parameters = list(
seed = seed_use,
assay = assay_use,

condition_order =
condition_order,

min_cells_per_group =
min_cells_per_group,

min_expression_pct =
min_expression_pct,

trim =
trim_use,

nboot =
nboot_use,

p_cutoff =
p_cutoff,

fdr_cutoff =
```

## Source PDF page 191

```text
fdr_cutoff,

balance_cell_numbers =
balance_cell_numbers,

biological_replicates =
FALSE
)
)

saveRDS(
complete_cellchat_results,
file = file.path(
outdir,
"complete_CellChat_analysis_seed54.rds"
)
)

############################################################
## 26.
############################################################

capture.output(
sessionInfo(),
file = file.path(
outdir,
"sessionInfo.txt"
)
)

############################################################
## 27.
############################################################

writeLines(
c(
"Statistical interpretation",
"",
"1. CK and FCR each contain only one biological sample.",
"2. Individual cells are subsamples, not independent biological replicates.",
"3. CellChat permutation P values quantify internal cell-label-level support.",
"4. BH-adjusted values do not test replicate-level CK versus FCR effects.",
"5. CK-FCR differences must be interpreted as descriptive pattern
```

## Source PDF page 192

```text
"6. Do not describe CK-FCR differences as statistically significant",
"  treatment effects.",
"7. Independent biological replicates are required for population-level",
"  differential communication inference."
),
con = file.path(
outdir,
"STATISTICAL_INTERPRETATION.txt"
)
)

############################################################
## 28.
############################################################

message("\n==============================================")
message("CellChat analysis completed.")
message(
"Output directory: ",
normalizePath(outdir)
)
message("")
message("Main output files:")
message(" wheat_CellChat_complete_results.xlsx")
message(" cellchat_CK_seed54.rds")
message(" cellchat_FCR_seed54.rds")
message(" complete_CellChat_analysis_seed54.rds")
message(" CK_retained_strength_circle.pdf")
message(" FCR_retained_strength_circle.pdf")
message(" CK_retained_LR_count_circle.pdf")
message(" FCR_retained_LR_count_circle.pdf")
message(" CK_FCR_descriptive_delta_heatmap.pdf")
message(" pathway_strength_CK_FCR.pdf")

if (nrow(bubble_data) > 0) {
message(" top_retained_LR_routes.pdf")
}

message(" LR_gene_availability.tsv")
message(" STATISTICAL_INTERPRETATION.txt")
message(" sessionInfo.txt")
message("==============================================")
```
