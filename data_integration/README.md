# Data Integration

## Overview

Integrates multimodal single-cell data (RNA, ATAC, mitochondrial variants) across samples and merges with gRNA assignments to generate a unified Seurat object.

**Input:**
- Cell Ranger multi output (RNA + ATAC)
- mgatk mitochondrial variant calls
- gRNA assignment metadata (from gRNA_assignment workflow)

**Output:**
- Integrated Seurat object with RNA, ATAC, and mito assays
- Combined metadata with gRNA/gene annotations
- Peak sets and dimensionality reductions

---

## Workflow

| Script | Description | Parameters |
|:-------|:------------|:-----------|
| [`data_integration.R`](./data_integration.R) | Integrates RNA and ATAC data across samples, generates combined peak sets, performs dimensionality reduction | Algorithm: SLM (3), resolution: 2, RNA dims: 1:20, ATAC dims: 2:40 |
| [`data_integration_add_guide_info.R`](./data_integration_add_guide_info.R) | Adds gRNA and target gene annotations to integrated object | Requires gRNA assignment completion (Steps 1-3) |

---

## Analysis Details

**Integration Steps:**

1. **Peak Calling:**
   - Generate combined peak sets using MACS2
   - Create common reference for ATAC integration

2. **RNA Integration:**
   - Process multiple samples
   - Normalize and scale data
   - Dimensionality reduction (PCA, UMAP)
   - Clustering: SLM algorithm, resolution 2.0

3. **ATAC Integration:**
   - Create fragment objects
   - TF-IDF normalization
   - Latent semantic indexing (LSI)
   - Dimensionality reduction (dims 2-40)

4. **Mitochondrial Variants:**
   - Load mgatk output
   - Calculate heteroplasmy levels
   - Add to Seurat metadata

5. **gRNA Annotation:**
   - Match cell barcodes to gRNA assignments
   - Add guide and target gene columns
   - Export combined metadata

**Key Parameters:**
- Clustering algorithm: 3 (SLM)
- Clustering resolution: 2.0
- RNA PCA dimensions: 1:20
- RNA UMAP dimensions: 1:50
- ATAC LSI dimensions: 2:40
- Peak feature cutoff: ≥10

---

## Dependencies

| Tool | Version | Purpose |
|:-----|:--------|:--------|
| [R](https://www.r-project.org/) | ≥4.0 | Statistical computing |
| [Seurat](https://satijalab.org/seurat/) | ≥4.0 | Single-cell integration |
| [Signac](https://stuartlab.org/signac/) | ≥1.6 | ATAC-seq analysis |
| [MACS2](https://github.com/macs3-project/MACS) | ≥2.2 | Peak calling |
| [tidyverse](https://www.tidyverse.org/) | ≥1.3 | Data manipulation |

---

## Prerequisites

**Before running data_integration.R:**
- Cell Ranger multi output directories for all samples
- mgatk analysis completed

**Before running data_integration_add_guide_info.R:**
- Completed gRNA assignment workflow (Steps 1-3)
- Integrated Seurat object from `data_integration.R`

---

## Expected Output

**Integrated Seurat Objects:**
- `project_RNA_ATAC_MGATK_combined.rds` - Multimodal object (RNA + ATAC + mito)
- `project_RNA_ATAC_MGATK_GUIDES_combined.rds` - Multimodal object with gRNA annotations

**Metadata:**
- `project_mtDNA_heteroplasmy_with_guides_results.csv` - Combined metadata with heteroplasmy and gRNA assignments

**Analysis Files:**
- Combined peak sets (BED format)
- Dimensionality reductions (PCA, LSI, UMAP)
- Clustering assignments
