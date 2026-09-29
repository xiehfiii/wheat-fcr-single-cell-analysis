# Wheat leaf-sheath FCR single-cell analysis code

Code associated with the CK and Fusarium crown rot (FCR) wheat leaf-sheath
single-cell analysis. This repository contains code only: no sequencing data,
Seurat objects, source tables, rendered figures, or credentials.

## Organization

| Directory | Analysis block |
| --- | --- |
| `analysis/00_preprocessing_QC` | QC and DoubletFinder preparation |
| `analysis/01_atlas_and_annotation` | Atlas source export, cell-type labels, and annotation stability |
| `analysis/02_regulatory_network` | SCENIC-source export, Figure 2 assembly, and SCENIC/virtual-knockout robustness |
| `analysis/03_pseudotime_BEAM` | Epidermal extraction, Monocle2 sensitivity analyses, trajectory, BEAM, and enrichment |
| `analysis/04_cell_communication` | Wheat ligand-receptor/CellChat analysis and communication figures |
| `analysis/05_protoplasting` | Protoplasting-impact analysis and plotting |

## Environments and inputs

- The observed `scrna` software and package versions are listed in
  [`environment/README.md`](environment/README.md).
- `scrna`: Seurat, annotation, regulatory-network inputs, and CellChat.
- `monocle2`: trajectory and BEAM.
- Stochastic finalized workflows use seed 54 where applicable.
- Inputs such as `merged.rds`, count matrices, SCENIC AUCell matrices, and
  ligand-receptor evidence tables are deliberately not included here.
- Scripts inherited from the original server may contain paths under
  `/home/xiehangfei/scRNAseq`; adjust those paths before rerunning. Some
  figure-assembly scripts accept a project directory as an argument.

The scripts were reorganized and checked for non-English text without
changing analytic thresholds, gene IDs, or statistical definitions.
CK and FCR each represent one library in the source
analysis; cell-level resampling is a sensitivity check, not a biological
replicate-based treatment test.

## Example commands

```bash
Rscript analysis/01_atlas_and_annotation/prepare_submission_source_data.R \
  /path/to/project /path/to/output
python analysis/02_regulatory_network/extract_scenic_submission.py \
  --input /path/to/auc_mtx.csv --output /path/to/auc.tsv.gz
Rscript analysis/01_atlas_and_annotation/build_suppfigure2_annotation_stability.R \
  /path/to/project /path/to/output
```
