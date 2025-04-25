# Cellular and Molecular Insights into Asymptomatic and Typical Alzheimer's Disease Using Single-Cell Transcriptomics

## Project Overview
This project investigates molecular and cellular mechanisms underlying cognitive resilience and vulnerability during Alzheimer's disease progression. Single-nucleus RNA-seq data from human prefrontal cortex samples are analyzed to compare Asymptomatic Alzheimer's Disease (AsymAD) and Typical Alzheimer's Disease (TAD) groups. Analyses include differential gene expression (DEG), pathway enrichment, cell-type and subtype composition analysis, pseudotime trajectory modeling, and cognitive association studies.

## Repository Structure
- `/` — R and Python scripts for preprocessing, DEG analysis, enrichment, pseudotime inference, and figure generation.
- `/data/` — Processed metadata, pseudotime scores, and summarized DEG results will be uploaded here
- `/figures/` — Final figures including UMAPs, heatmaps, Circos plots, and pseudotime visualizations will be uploaded here
- `/notebooks/` — Jupyter notebooks for reproducible trajectory analysis will be uploaded here

## Analysis Workflow
- Differential expression analysis (AsymAD vs TAD) across major brain cell types.
- Gene Ontology (GO), SynGO, and CORUM protein complex enrichment analyses.
- Cross-comparison of DEG signatures across pathology and cognition measures.
- Validation of DEGs through external datasets (De Jager DLPFC, SEA-AD MTG) and permutation testing.
- Cell-type and fine-grained neuronal subtype composition analysis.
- Trajectory inference using Partition-based Graph Abstraction (PAGA) and pseudotime calculation.
- Correlation of gene expression with cognitive performance metrics.

## Installation
Required software:
- **R** (≥ 4.2) — packages: `Seurat`, `edgeR`, `limma`, `muscat`, `ggplot2`

You can install necessary R packages with:
```
install.packages(c("Seurat", "edgeR", "limma", "clusterProfiler", "ggplot2"))
```
