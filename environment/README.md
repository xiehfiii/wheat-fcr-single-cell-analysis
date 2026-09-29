# Software and package inventory

Inventory date: 2026-09-29. Sources: the `scrna` environment on the analysis
server and package calls in `/home/xiehangfei/scRNAseq`, excluding the bundled
Cell Ranger source code. These are versions **observed at inventory time**,
not necessarily the historical versions used to generate every result.

## Platform and command-line software

| Software | Observed version | Usage evidence |
| --- | --- | --- |
| Ubuntu | 20.04.6 LTS, x86_64 | Server operating system |
| R | 4.4.3 | `scrna/bin/Rscript`; R analysis and figure scripts |
| Python | 3.10.19 | `scrna/bin/python`; Python preparation and RBH scripts |
| BLAST+ (`blastp`, `makeblastdb`) | 2.17.0+ | Executables in `scrna`; explicitly called by protein reciprocal-best-hit scripts |
| Cell Ranger | 9.0.1 | Bundled under `scRNAseq/cellranger-9.0.1`; a drought-sample pipeline invocation was found, **not** proof that this version generated the CK/FCR matrices |
| STAR | 2.7.11b | Installed in `scrna`; CK/FCR-specific invocation not located in the scanned scripts |
| SAMtools | 1.22.1 | Installed in `scrna`; CK/FCR-specific invocation not located |
| FastQC | 0.12.1 | Installed in `scrna`; CK/FCR-specific invocation not located |
| BEDTools | 2.31.1 | Installed in `scrna`; CK/FCR-specific invocation not located |

Presence of a binary is not evidence that it was used for the manuscript's
CK/FCR libraries. In particular, this inventory did not establish a
CK/FCR-specific `cellranger count` invocation.

## R packages referenced in project scripts and installed in `scrna`

The scan includes `library()`, `require()`, `requireNamespace()`, and `pkg::`
references across the project, including exploratory and supplementary scripts.
It is broader than the dependencies of any single final figure.

| Package | Version | Package | Version |
| --- | --- | --- | --- |
| AnnotationDbi | 1.68.0 | AUCell | 1.28.0 |
| batchelor | 1.22.0 | CellChat | 2.2.0.9001 |
| circlize | 0.4.16 | cluster | 2.1.8.1 |
| clusterProfiler | 4.14.0 | clustree | 0.5.1 |
| ComplexHeatmap | 2.22.0 | data.table | 1.18.0 |
| DoubletFinder | 2.0.6 | dplyr | 1.1.4 |
| edgeR | 4.4.2 | future | 1.67.0 |
| ggalluvial | 0.12.6 | ggplot2 | 4.0.0 |
| ggraph | 2.2.2 | ggrepel | 0.9.6 |
| glmGamPoi | 1.18.0 | grid | 4.4.3 |
| GSEABase | 1.68.0 | gtools | 3.9.5 |
| httr | 1.4.7 | igraph | 2.2.1 |
| jsonlite | 2.0.0 | magick | 2.9.0 |
| magrittr | 2.0.4 | Matrix | 1.7-4 |
| monocle | 2.34.0 | openxlsx | 4.2.8.1 |
| org.Ta.eg.db | 0.1 | patchwork | 1.3.2 |
| pdftools | 3.9.0 | pheatmap | 1.0.13 |
| plotly | 4.11.0 | png | 0.1-8 |
| ragg | 1.5.0 | RColorBrewer | 1.1-3 |
| readr | 2.1.5 | readxl | 1.4.5 |
| remotes | 2.5.0 | rlang | 1.1.6 |
| rtracklayer | 1.66.0 | scales | 1.4.0 |
| Seurat | 5.3.1 | SeuratObject | 5.2.0 |
| SoupX | 1.6.2 | stringr | 1.6.0 |
| svglite | 2.2.2 | tibble | 3.3.0 |
| tidydr | 0.0.6 | tidyHeatmap | 1.13.1 |
| tidyr | 1.3.1 | tidyverse | 2.0.0 |
| UCell | 2.10.1 |  |  |

## Python packages referenced in project scripts and installed in `scrna`

| Package | Version |
| --- | --- |
| pandas | 2.3.3 |
| requests | 2.32.5 |

Python standard-library imports are omitted. Trajectory/BEAM analyses may
also have run in the separate `monocle2` environment, and SCENIC in the
separate `scenic` environment. This page documents **`scrna` only** and must
not be cited as an audit of those other environments.
