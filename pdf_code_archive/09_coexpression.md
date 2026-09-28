# hdWGCNA coexpression

PDF transcription only. The PDF may wrap code or omit glyphs;
use the curated scripts for executable analysis.

## Source PDF page 126

```text
df_BP_dot$GeneRatio_Num <- sapply(strsplit(df_BP_dot$GeneRatio,"/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
df_BP_dot$Description <- sapply(df_BP_dot$Description, function(x) {
x_clean <- gsub(",", "", x)
paste(head(strsplit(x_clean, " ")[[1]], 5), collapse = " ")
})

p_go_dot <- ggplot(df_BP_dot, aes(x = GeneRatio_Num, y = reorder(Description, GeneRatio_Num), size = Count, color = p.adjust)) +
geom_point() +
scale_color_gradientn(colors = dot_cols, trans = "reverse", name ="adj.P") +
scale_size(range = c(3, 8)) +
theme_bw() +
theme(panel.grid = element_blank(),
plot.title = element_text(hjust = 0.5),
axis.text.y = element_text(size = 10)) +
labs(x = "Gene Ratio", y = NULL, title = "GO Enrichment (Biological Process)")

ggsave(filename = file.path(out_dir, "Knockout_GO_BP_dotplot.pdf"), plot = p_go_dot, width = 10, height = 8)
print(">>>     Dotplot")
}
```

## Source PDF page 127

```text
# ============================================================
# hdWGCNA
# ============================================================

# ---- 1.         ----
library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(Seurat)

# ---- 2.            /           ----
target_group <- "Epidermis cells I"
target_gene <- "TraesCS3D02G094200"
n_threads  <- 8

# ---- 3.    ----
set.seed(54)
enableWGCNAThreads(nThreads = n_threads)
Idents(merged) <- "celltype"

# ---- 4.      ----
merged <- SeuratObject::UpdateSeuratObject(merged)
seurat_obj <- SetupForWGCNA(
merged,
gene_select = "fraction",
fraction = 0.05,
wgcna_name = "fcr"
)

# ---- 5.       ----
seurat_obj <- MetacellsByGroups(
seurat_obj = seurat_obj,
group.by = c("celltype", "orig.ident"),
reduction = 'pca',
k = 25,
max_shared = 10,
ident.group = 'celltype',
min_cells = 100
)
seurat_obj <- NormalizeMetacells(seurat_obj)

```

## Source PDF page 128

```text
# ---- 6.                    ----
seurat_obj <- SetDatExpr(
seurat_obj,
group_name = target_group,
group.by = 'celltype',
assay = 'RNA',
slot = 'data'
)

# ---- 7.       ----
seurat_obj <- TestSoftPowers(
seurat_obj,
powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
networkType = 'unsigned'
)

#
plot_list <- PlotSoftPowers(seurat_obj)
pdf("soft_power_analysis.pdf", width = 10, height = 8)
print(wrap_plots(plot_list, ncol = 2))
dev.off()

#
power_table <- GetPowerTable(seurat_obj)
print(head(power_table))

# ---- 8.                       ----
seurat_obj <- ConstructNetwork(
seurat_obj,
setDatExpr = FALSE,
tom_outdir = "TOM",
tom_name = target_group,
overwrite_tom = TRUE
)

# ---- 9.            ----
pdf(paste0(target_group, "_hdWGCNA_dendrogram.pdf"), width = 10, height = 6)
PlotDendrogram(seurat_obj, main = paste(target_group, 'hdWGCNA Dendrogram'))
dev.off()

# ---- 10.             kME ----
seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)

```

## Source PDF page 129

```text
modules <- GetModules(seurat_obj)

# ---- 11.                 ----
target_info <- modules[modules$gene_name == target_gene, ]
print(target_info)

target_module <- target_info$module
cat("\n    ", target_gene, " :", target_module, " \n")

if (target_module != "grey") {
#             kME
kme_col <- grep(paste0("kME.*", target_module), colnames(modules),
value = TRUE, ignore.case = TRUE)

#               kME
co_expressed <- modules %>%
filter(module == target_module) %>%
select(gene_name, module, all_of(kme_col)) %>%
arrange(desc(.data[[kme_col]]))

cat("\n---  20          ---\n")
print(head(co_expressed, 20))

write.csv(co_expressed,
paste0("co_expressed_", target_module, "_module.csv"),
row.names = FALSE)
} else {
cat("\n             grey                    \n")
}

# ---- 12. Hub     ----
hub_genes <- GetHubGenes(seurat_obj, n_hubs = 20)
print(head(hub_genes))
write.csv(hub_genes, "hub_genes_per_module.csv", row.names = FALSE)
```

## Source PDF page 130

```text
install.packages(c("WGCNA", "devtools"))
devtools::install_github('smorabit/hdWGCNA', ref='dev')
library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(Seurat)

set.seed(54)
enableWGCNAThreads(nThreads = 8)
Idents(merged) <- "celltype"

merged <- SeuratObject::UpdateSeuratObject(merged)
seurat_obj <- SetupForWGCNA(
merged,
gene_select = "fraction",
fraction = 0.05,
wgcna_name = "fcr"
)
seurat_obj <- MetacellsByGroups(
seurat_obj = seurat_obj,
group.by = c("celltype", "orig.ident"),
reduction = 'pca',
k = 25,
max_shared = 10,
ident.group = 'celltype',
min_cells = 100
)

seurat_obj <- NormalizeMetacells(seurat_obj)

#####4.
seurat_obj <- SetDatExpr(
seurat_obj,
group_name = "Epidermis cells I", # the name of the group of interest in the group.by column
group.by='celltype', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
assay = 'RNA', # using RNA assay
slot = 'data' # using normalized data
)

```

## Source PDF page 131

```text
# #      group_names
# seurat_obj <- SetDatExpr(
#  seurat_obj,
#  group_name = c("Treg", "Tm"),
#  group.by='celltype'
# )

#####4.1
# Test different soft powers:
seurat_obj <- TestSoftPowers(
seurat_obj,
powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
networkType = 'unsigned' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)
# 1   1 0.4193109 -2.829398   0.9708088 551.0532441 532.1816874 1079.59624
# 2   2 0.7722457 -2.846229   0.9935455 90.8116702 80.0299812 331.23789
# 3   3 0.8469701 -2.813601   0.9812133 19.9509707 15.0395890 127.98068
# 4   4 0.8914758 -2.532259   0.9884110 5.4419966  3.277504357.54649
# 5   5 0.9024030 -2.199419   0.9834531 1.7820619  0.792562628.86319
# 6   6 0.9511793 -1.813044   0.9753076 0.6868952  0.208334115.72666

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

# check
power_table <- GetPowerTable(seurat_obj)
head(power_table)
#  Power  SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
# 1   1 0.02536182 3.273051    0.9541434 4370.9149 4379.0629 4736.

# 2   2 0.11091306 -3.571441   0.8008960 2322.5480 2286.2454 2871.

# 3   3 0.50454728 -4.960822   0.8035027 1286.6453 1241.8414 1898.

# 4   4 0.79569568 -4.812735   0.9183803 740.0525 697.1193 1338.
```

## Source PDF page 132

```text

# 5   5 0.86641323 -4.110731   0.9517671 440.6141 402.5530 985.

# 6   6 0.88593187 -3.582879   0.9624951 270.9020 237.8831 750.


#####4.2
#               construcNetwork
# construct co-expression network:
seurat_obj <- ConstructNetwork(
seurat_obj,
soft_power = 4, #     4                    3
tom_outdir = "TOM",
tom_name = 'Epidermis cells I', # name of the topoligical overlap matrix written to disk
overwrite_tom = TRUE #
)

#    WGCNA
# “  ”
pdf("Epidermis_cells_hdWGCNA_dendrogram.pdf", width = 10, height = 6)
PlotDendrogram(seurat_obj, main = 'Epidermis cells hdWGCNA Dendrogram')
dev.off()

#                (topoligcal overlap matrix TOM)
# TOM <- GetTOM(seurat_obj)
# TOM

#####4.3
#       ScaleData    harmony   :
# seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

#                    MEs
seurat_obj <- ModuleEigengenes(
seurat_obj,
group.by.vars="orig.ident"
)

#          :
#       MEs   Harmony                     (hMEs)
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
#MEs <- GetMEs(seurat_obj, harmonized=FALSE)
```

## Source PDF page 133

```text

#####4.4
#                 (kME)
#
seurat_obj <- ModuleConnectivity(
seurat_obj,
group.by = 'celltype',
group_name = 'Epidermis cells I'
)

#
seurat_obj <- ResetModuleNames(
seurat_obj,
new_name = "Treg_NEW"
)

#          kME
p <- PlotKMEs(seurat_obj, ncol=4)
p

#####4.5
#            :
#
#
modules <- GetModules(seurat_obj) %>%
subset(module != 'grey')

#       :
head(modules[,1:6])
#        gene_name  module    color  kME_grey kME_Treg_NEW1 kME_Treg_NEW2
# ISG15     ISG15 Treg_NEW1 turquoise 0.09485063 0.316180060.2177907
# TNFRSF18 TNFRSF18 Treg_NEW1 turquoise 0.12119087 0.398862460.4605542
# TNFRSF4  TNFRSF4 Treg_NEW1 turquoise 0.08844463 0.359223370.3728684
# SDF4       SDF4 Treg_NEW2   black 0.11518097 0.112121550.1883993
# B3GALT6  B3GALT6 Treg_NEW3 purple 0.03314139 0.086108110.1067775
# AURKAIP1 AURKAIP1 Treg_NEW1 turquoise 0.09062613 0.262448270.1252306

#
#       kME      N           ,       10
```

## Source PDF page 134

```text
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
head(hub_df)
#  gene_name  module     kME
# 1   GAPDH Treg_NEW1 0.6160237
# 2  S100A4 Treg_NEW1 0.5886924
# 3    MYL6 Treg_NEW1 0.5558792
# 4  TMSB10 Treg_NEW1 0.5371290
# 5    IL32 Treg_NEW1 0.5161320
# 6  ARPC1B Treg_NEW1 0.5138853

#
qsave(seurat_obj, 'hdWGCNA_object.qs')

#####4.6 hub   siganture
#          25         kME
#   UCell
library(UCell)
seurat_obj <- ModuleExprScore(
seurat_obj,
n_genes = 25,
method='UCell' # Seurat (AddModuleScore)
)

#####5.
#####
#         hMEs
plot_list <- ModuleFeaturePlot(
seurat_obj,
features='hMEs', # plot the hMEs
order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=4)

#####       hub
#         hub scores
plot_list <- ModuleFeaturePlot(
seurat_obj,
features='scores', # plot the hub gene scores
order='shuffle', # order so cells are shuffled
ucell = TRUE # depending on Seurat vs UCell for gene scoring
)
# stitch together with patchwork
wrap_plots(plot_list, ncol=4)

```

## Source PDF page 135

```text
#####
#
seurat_obj$cluster <- do.call(rbind, strsplit(as.character(seurat_obj$orig.ident), ' '))[,1]

ModuleRadarPlot(
seurat_obj,
group.by = 'cluster',
barcodes = seurat_obj@meta.data %>%
subset(celltype == 'Treg') %>%
rownames(),
axis.label.size=4,
grid.label.size=4
)

#
ModuleCorrelogram(seurat_obj)

#####
# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'celltype')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
RotatedAxis() +
scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p

#   ModuleNetworkPlot         50(       )  hub gene
ModuleNetworkPlot(
seurat_obj,
outdir='ModuleNetworks', # new folder name
n_inner = 20, # number of genes in inner ring
n_outer = 30, # number of genes in outer ring
n_conns = Inf, # show all of the connections
plot_size=c(10,10), # larger plotting area
```

## Source PDF page 136

```text
vertex.label.cex=1 # font size
)

options(future.globals.maxSize = 5 * 1024^3) # 5GB

# hubgene network(      )
HubGeneNetworkPlot(
seurat_obj,
n_hubs = 2,
n_other=2,
edge_prop = 0.75,
mods = 'all'
)

#
g <- HubGeneNetworkPlot(seurat_obj, return_graph=TRUE)
# get the list of modules:
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
# hubgene network
HubGeneNetworkPlot(
seurat_obj,
n_hubs = 2,
n_other= 2,
edge_prop = 0.75,
mods = mods[1:5] # only select 5 modules
)

seurat_obj <- RunModuleUMAP(
seurat_obj,
n_hubs = 10, # number of hub genes to include for the UMAP embedding
n_neighbors=15, # neighbors parameter for UMAP
min_dist=0.1 # min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
geom_point(
color=umap_df$color, # color each point by WGCNA module
size=umap_df$kME*2 # size of each point based on intramodular connectivity
) +
umap_theme()
```

## Source PDF page 137

```text


ModuleUMAPPlot(
seurat_obj,
edge.alpha=0.25,
sample_edges=TRUE,
edge_prop=0.1, # proportion of edges to sample (20% here)
label_hubs=2 ,# how many hub genes to plot per module?
keep_grey_edges=FALSE
)

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(GeneOverlap)
library(qs)

#dir.create("14-hdWGCNA")
#setwd("14-hdWGCNA")

seurat_obj <- qread("hdWGCNA_object.qs")

#   enrichr databases
dbs <- c('GO_Biological_Process_2023',
'GO_Cellular_Component_2023',
'GO_Molecular_Function_2023')

#
seurat_obj <- RunEnrichr(
seurat_obj,
dbs=dbs,
max_genes = 100 # use max_genes = Inf to choose all genes
)

#
enrich_df <- GetEnrichrTable(seurat_obj)

#
head(enrich_df)

# make GO term plots:
EnrichrBarPlot(
```

## Source PDF page 138

```text
seurat_obj,
outdir = "enrichr_plots", # name of output directory
n_terms = 10, # number of enriched terms to show (sometimes more areshown if there are ties)
plot_size = c(5,7), # width, height of the output .pdfs
logscale=TRUE # do you want to show the enrichment as a log scale?
)

# enrichr dotplot
EnrichrDotPlot(
seurat_obj,
mods = "all", # use all modules (default)
database = "GO_Biological_Process_2023", # this must match one of the dbs used previously
n_terms=2, # number of terms per module
term_size=8, # font size for the terms
p_adj = FALSE # show the p-val or adjusted p-val?
) + scale_color_stepsn(colors=rev(viridis::magma(256)))

library(fgsea)

# load the GO Biological Pathways file (downloaded from EnrichR website)
pathways <- fgsea::gmtPathways('GO_Biological_Process_2023.txt')

# optionally, remove the GO term ID from the pathway names to make thedownstream plots look cleaner
names(pathways) <- stringr::str_replace(names(pathways), " \\s*\\([^\\)]+\\)", "")


# get the modules table and remove grey genes
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# rank by Treg_NEW1 genes only by kME
cur_mod <- 'Treg_NEW1'
modules <- GetModules(seurat_obj) %>% subset(module == cur_mod)
cur_genes <- modules[,(c('gene_name', 'module', paste0('kME_', cur_mod)))]
ranks <- cur_genes$kME; names(ranks) <- cur_genes$gene_name
ranks <- ranks[order(ranks)]

# run fgsea to compute enrichments
gsea_df2 <- fgsea::fgsea(
pathways = pathways,
stats = ranks,
```

## Source PDF page 139

```text
minSize = 3,
maxSize = 500
)

#
top_pathways <- gsea_df2 %>%
subset(pval < 0.05) %>%
slice_max(order_by=NES, n=25) %>%
.$pathway

plotGseaTable(
pathways[top_pathways],
ranks,
gsea_df2,
gseaParam=0.5,
colwidths = c(10, 4, 1, 1, 1)
)

# name of the pathway to plot
selected_pathway <- 'Cellular Respiration'
plotEnrichment(
pathways[[selected_pathway]],
ranks
) + labs(title=selected_pathway)

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(magrittr)
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(JASPAR2020)
library(JASPAR2024)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(xgboost)
library(JASPAR2024)
library(RSQLite)
library(EnsDb.Hsapiens.v86)
library(qs)

```

## Source PDF page 140

```text
#dir.create("14-hdWGCNA")
#setwd("14-hdWGCNA")

seurat_obj <- qread("hdWGCNA_object.qs")

# JASPAR 2020
pfm_core <- TFBSTools::getMatrixSet(
x = JASPAR2020,
opts = list(collection = "CORE",
tax_group = 'vertebrates',
all_versions = FALSE)
)

# JASPAR 2024 (not used for this tutorial)
# JASPAR2024 <- JASPAR2024()
# sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))
# pfm_core <- TFBSTools::getMatrixSet(
#  x = sq24,
#  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
# )

#   motif
seurat_obj <- MotifScan(
seurat_obj,
species_genome = 'hg38',
pfm = pfm_core,
EnsDb = EnsDb.Hsapiens.v86
)

#   motif df:
motif_df <- GetMotifs(seurat_obj)

#      TFs,
tf_genes <- unique(motif_df$gene_name)
modules <- GetModules(seurat_obj)
nongrey_genes <- subset(modules, module != 'grey') %>% .$gene_name
genes_use <- c(tf_genes, nongrey_genes)

# update the gene list and re-run SetDatExpr
seurat_obj <- SetWGCNAGenes(seurat_obj, genes_use)
seurat_obj <- SetDatExpr(seurat_obj, group.by = 'celltype', group_name='Treg')

# define model params:
model_params <- list(
```

## Source PDF page 141

```text
objective = 'reg:squarederror',
max_depth = 1,
eta = 0.1,
nthread=16,
alpha=0.5
)

#
seurat_obj <- ConstructTFNetwork(seurat_obj, model_params)
results <- GetTFNetwork(seurat_obj)
head(results)
#      tf  gene     Gain     Cover Frequency       Cor
# 1 ZKSCAN1 FOXD1 0.11708119 0.04695019 0.04695019 -0.19391353
# 2  NFIL3 FOXD1 0.10931756 0.05379589 0.05379589 0.19738680
# 3 ZNF652 FOXD1 0.09789632 0.07897635 0.07897635 -0.17072678
# 4  NR4A1 FOXD1 0.09624640 0.04498028 0.04498028 0.18091337
# 5  ZNF24 FOXD1 0.05250378 0.02255133 0.02255133 -0.09566174
# 6  NFKB2 FOXD1 0.05148296 0.06481095 0.06481095 0.16717336

#   “A”           10 TF
#   “B”
#   “C”                  TF-
seurat_obj <- AssignTFRegulons(
seurat_obj,
strategy = "A", # B  C
reg_thresh = 0.01,
n_tfs = 10
)

#
#             TF    (   )      (  )
p1 <- RegulonBarPlot(seurat_obj, selected_tf='ZNF652')
p2 <- RegulonBarPlot(seurat_obj, selected_tf='NFKB2', cutoff=0.15)
p1 | p2

#   regulons
seurat_obj <- RegulonScores(
seurat_obj,
target_type = 'positive',
ncores=8
)

#   regulons
seurat_obj <- RegulonScores(
seurat_obj,
target_type = 'negative',
```
