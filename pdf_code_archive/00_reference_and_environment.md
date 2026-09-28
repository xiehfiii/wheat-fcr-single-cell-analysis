# Reference preparation and environments

PDF transcription only. The PDF may wrap code or omit glyphs;
use the curated scripts for executable analysis.

## Source PDF page 1

```text
nohup cellranger mkref --memgb=200 --genome=wheat --fasta=./Triticum_aestivum.fasta --genes=./Triticum_aestivum.gtf &

nohup cellranger count --id=ck_fcr_10000cells --transcriptome=/home/xiehangfei/scRNAseq/ChineseSpring_mkref_ALL --fastqs=/home/xiehangfei/scRNAseq/CK-FCR --sample=CK-FCR --expect-cells=10000 --localcores=20 --localmem=200 --include-introns=true --create-bam=true &

nohup cellranger count --id=fcr_10000cells --transcriptome=/home/xiehangfei/scRNAseq/ChineseSpring_mkref_ALL --fastqs=/home/xiehangfei/scRNAseq/FCR --sample=FCR --expect-cells=10000 --localcores=20 --localmem=200--include-introns=true --create-bam=true &

#fastqs:              ak58_S1_L001_R1_001.fastq.gz;ak58_S1_L001_R2_001.fastq.gz
mamba create -n scrna r-base=4.4.3 r-seurat r-tidyverse -c conda-forge-c bioconda -c defaults r-rcpp r-essentials git wget unzip -y
```

## Source PDF page 3

```text
mamba install bioconductor-rtracklayer -c conda-forge -c bioconda -y
```

## Source PDF page 4

```text
BiocManager::install("clustree")
BiocManager::install("tidyHeatmap")
install.packages("remotes")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
BiocManager::install('glmGamPoi')
install.packages('SoupX')
install.packages("Matrix")
BiocManager::install("batchelor")
remotes::install_github("cole-trapnell-lab/monocle3")
BiocManager::install("AUCell")

library(remotes)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(batchelor)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(clustree)
library(plotly)
library(rtracklayer)
library(tidyHeatmap)
library(future)
library(DoubletFinder)
library(SoupX)
library(glmGamPoi)
library(ggrepel)
library(tidydr)
library(RColorBrewer)
library(AUCell)
library(GSEABase)
library(pheatmap)
library(tibble)
library(clusterProfiler)
library(edgeR)
library(org.Ta.eg.db)
library(igraph)
library(ComplexHeatmap)
library(circlize)
```
