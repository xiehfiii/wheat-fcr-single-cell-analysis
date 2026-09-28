# RNA velocity

PDF transcription only. The PDF may wrap code or omit glyphs;
use the curated scripts for executable analysis.

## Source PDF page 80

```text
mamba install -c conda-forge -c bioconda loompy velocyto.py -y --channel-priority flexible

umap_coords <- data.frame(
cell_id = colnames(merged),
umap_1 = merged@reductions$umap@cell.embeddings[,1],
umap_2 = merged@reductions$umap@cell.embeddings[,2]
)

colnames(umap_coords)
write.csv(umap_coords, "umap_coords_merged.csv", row.names = FALSE)
cluster_labels <- data.frame(
cell_id = colnames(merged),
cluster = merged$celltype
)

table(cluster_labels$cluster)
write.csv(cluster_labels, "seurat_clusters.csv", row.names = FALSE)
head(umap_coords$cell_id)
velocyto run10x /home/xiehangfei/scRNAseq/merged_15000cells/ ./Triticum_aestivum.gtf
```

## Source PDF page 81

```text
#  loompy (python)
python - <<'PY'
import loompy
ds = loompy.connect("merged_15000cells.loom", mode="r")
print("Layers:", ds.layers.keys()) # spliced/unspliced/ambiguous
print("Shape:", ds.shape)
ds.close()
PY



vim scvelo_pipeline.py

#!/usr/bin/env python3
"""
Robust scVelo pipeline for velocyto .loom -> scVelo analysis.
Supports importing Seurat UMAP/clustering for accurate cell mapping.

Usage:
# Auto-compute UMAP/louvain:
python3 scvelo_pipeline.py --loom data.loom --outdir results --plot

# Use Seurat UMAP/clusters:
python3 scvelo_pipeline.py --loom data.loom --outdir results --plot \
--umap seurat_umap_coords.csv --clusters seurat_clusters.csv
"""

import os, sys, argparse
import numpy as np
import warnings
warnings.filterwarnings('ignore')

#
try:
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import anndata
import pandas as pd
except ImportError as e:
print(f"Error: Missing required package: {e}")
print("Please install: conda install -c conda-forge scanpy scvelo matplotlib anndata pandas")
```

## Source PDF page 82

```text
sys.exit(1)

try:
import igraph
except ImportError:
print("Warning: python-igraph not found. Please run: conda install-c conda-forge python-igraph")
sys.exit(1)

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo')

parser = argparse.ArgumentParser()
parser.add_argument("--loom", required=True, help="input loom file from velocyto")
parser.add_argument("--outdir", default="scvelo_result", help="output directory")
parser.add_argument("--min_shared_counts", type=int, default=20, help="filter_and_normalize min_shared_counts")
parser.add_argument("--n_top_genes", type=int, default=2000, help="number of HVG for scVelo steps")
parser.add_argument("--min_unspliced_per_cell", type=int, default=100,help="min unspliced counts per cell to keep")
parser.add_argument("--ambiguous_gene_threshold", type=float, default=0.3, help="max fraction of ambiguous reads per gene to keep (<=)")
parser.add_argument("--plot", action="store_true", help="generate final velocity stream plot")
parser.add_argument("--umap", default=None, help="CSV Seurat UMAP    cell_id, umap_1, umap_2 ")
parser.add_argument("--clusters", default=None, help="CSV Seuratcell_id, cluster ")
args = parser.parse_args()

loomF = args.loom
outdir = args.outdir
os.makedirs(outdir, exist_ok=True)

print("="*60)
print(f"Reading loom: {loomF}")
print("="*60)

#   loom
adata = sc.read_loom(loomF)

#
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
```

## Source PDF page 83

```text
for layer in ['spliced','unspliced','ambiguous']:
if layer not in adata.layers:
print(f"Warning: layer '{layer}' not found in loom. Exiting.")
sys.exit(1)

#
for layer in ['spliced','unspliced','ambiguous']:
if not hasattr(adata.layers[layer], 'tocsr'):
adata.layers[layer] = scv.utils.as_sparse(adata.layers[layer])

# QC
spliced_sum = np.sum(adata.layers['spliced'], axis=1).A1
unspliced_sum = np.sum(adata.layers['unspliced'], axis=1).A1
amb_sum = np.sum(adata.layers['ambiguous'], axis=1).A1

adata.obs['spliced_counts'] = spliced_sum
adata.obs['unspliced_counts'] = unspliced_sum
adata.obs['ambiguous_counts'] = amb_sum
adata.obs['unspliced_ratio'] = adata.obs['unspliced_counts'] / (adata.obs['spliced_counts'] + 1e-9)

#   QC
print(f"Median spliced per cell: {np.median(spliced_sum):.1f}")
print(f"Median unspliced per cell: {np.median(unspliced_sum):.1f}")
print(f"Median unspliced/spliced ratio: {np.median(adata.obs['unspliced_ratio']):.4f}")

#    unspliced
keep_cells = adata.obs['unspliced_counts'] >= args.min_unspliced_per_cell
print(f"Keeping cells with >={args.min_unspliced_per_cell} unspliced:{adata.n_obs} -> {np.sum(keep_cells)}")
adata = adata[keep_cells].copy()

#    ambiguous
gene_spliced = np.array(adata.layers['spliced'].sum(axis=0)).ravel()
gene_unspliced = np.array(adata.layers['unspliced'].sum(axis=0)).ravel()
gene_amb = np.array(adata.layers['ambiguous'].sum(axis=0)).ravel()
gene_total = gene_spliced + gene_unspliced + gene_amb
gene_amb_frac = np.divide(gene_amb, gene_total, out=np.zeros_like(gene_amb, dtype=float), where=gene_total>0)

keep_genes = gene_amb_frac <= args.ambiguous_gene_threshold
print(f"Genes kept after ambiguous filter (<={args.ambiguous_gene_threshold}): {np.sum(keep_genes)} / {len(keep_genes)}")
```

## Source PDF page 84

```text
adata = adata[:, keep_genes].copy()

# scVelo
print("\n" + "="*60)
print("scVelo  ...")
print("="*60)
scv.pp.filter_and_normalize(adata, min_shared_counts=args.min_shared_counts, n_top_genes=args.n_top_genes)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

#   Seurat
if args.umap and args.clusters:
print("\n" + "="*60)
print("  Seurat UMAP      ...")
print("="*60)

#        UMAP
umap_df = pd.read_csv(args.umap, index_col=0)
umap_df = umap_df.loc[adata.obs_names]
adata.obsm['X_umap'] = umap_df[['umap_1', 'umap_2']].values
print(f"    UMAP  : {adata.obsm['X_umap'].shape}")

#
cluster_df = pd.read_csv(args.clusters, index_col=0)
cluster_df = cluster_df.loc[adata.obs_names]
adata.obs['seurat_clusters'] = cluster_df['cluster'].astype('category')
print(f"         : {adata.obs['seurat_clusters'].nunique()}")

#
color_col = 'seurat_clusters'
cluster_title = "Seurat Clusters"
print("   scvelo          Seurat  ")
else:
#      UMAP Louvain
print("\n" + "="*60)
print("  UMAP Louvain  ...")
print("="*60)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep='X_pca')
sc.tl.umap(adata, min_dist=0.3)
sc.tl.louvain(adata, resolution=0.6)

color_col = 'louvain'
cluster_title = "Louvain Clusters"
print(f" UMAP     : {adata.obsm['X_umap'].shape}")
```

## Source PDF page 85

```text
print(f" Louvain   : {adata.obs['louvain'].nunique()}")

# Stochastic velocity
print("\n" + "="*60)
print(" Stochastic     ...")
print("="*60)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

#   QC           scvelo save
qc_path = os.path.join(outdir, "velocity_stochastic_qc.png")
print(f"    QC ...")
plt.figure(figsize=(10, 8))
scv.pl.velocity_embedding_stream(adata, basis='umap', color=color_col,show=False)
plt.savefig(qc_path, dpi=300, bbox_inches="tight")
plt.close()
print(f" Stochastic QC   : {qc_path}")

# Dynamical model
print("\n" + "="*60)
print(" Dynamical           ...")
print("="*60)
try:
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.tl.latent_time(adata)

dyn_path = os.path.join(outdir, "velocity_dynamical_latent_time.png")
print(f"     Dynamical ...")
plt.figure(figsize=(10, 8))
scv.pl.velocity_embedding_stream(adata, basis='umap', color='latent_time', show=False)
plt.savefig(dyn_path, dpi=300, bbox_inches="tight")
plt.close()
print(f" Dynamical       : {dyn_path}")
except Exception as e:
print(f" Dynamical     : {e}")
print("      Stochastic  ")

# Save results
output_h5ad = os.path.join(outdir, "scvelo_analysis.h5ad")
adata.write(output_h5ad)
print(f"\n          : {output_h5ad}")
```

## Source PDF page 86

```text

# Export metadata
meta_path = os.path.join(outdir, "cell_metadata_scvelo.csv")
adata.obs.to_csv(meta_path)
print(f"            : {meta_path}")

#
if args.plot:
print("\n" + "="*60)
print("           ...")
print("="*60)

final_path = os.path.join(outdir, "velocity_stream_by_cluster.png")
print(f"         ...")
plt.figure(figsize=(14, 10))
scv.pl.velocity_embedding_stream(
adata,
basis="umap",
color=color_col,
palette="tab20",
legend_loc="right margin",
title=f"RNA Velocity Stream Plot (Colored by {cluster_title})",
linewidth=1.5,
arrow_size=2,
density=2,
show=False,
size=50
)
plt.savefig(final_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"         : {final_path}")

print("\n" + "="*60)
print("                 :", outdir)
print("="*60)
print("    :")
print(f" ├──     : {output_h5ad}")
print(f" ├──    :  {meta_path}")
print(f" ├── QC :   {qc_path}")
if os.path.exists(os.path.join(outdir, "velocity_dynamical_latent_time.png")):
print(f" ├── Dynamical : {outdir}/velocity_dynamical_latent_time.png")
if args.plot:
print(f" └──    :   {final_path}")

```
