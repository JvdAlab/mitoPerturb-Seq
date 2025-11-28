# Tricycle Cell Cycle Analysis

## Overview

Infers continuous cell cycle positions for single cells using the tricycle algorithm and validates with canonical cell cycle markers.

**Input:**
- Integrated Seurat object with RNA assay and gRNA annotations

**Output:**
- Cell cycle positions (θ, in radians: 0-2π)
- Phase assignments (G1, G1/S, S, G2, G2/M, M)
- Expression profiles of cell cycle markers along θ

---

## Scripts

| Script | Description | Parameters |
|:-------|:------------|:-----------|
| [`tricycle_analysis.R`](./tricycle_analysis.R) | Estimates cell cycle position using tricycle and validates with 36 canonical markers | Species: mouse, gene name type: SYMBOL |

---

## Analysis Details

**Cell Cycle Position Estimation:**
- Algorithm: tricycle (transfer learning + periodic LOESS)
- Output: θ (theta) values in radians (0-2π)
- Reference: Pre-trained model for mouse scRNA-seq

**Cell Cycle Markers (36 genes):**
- **G1:** Ccnd1, Cdk6, Cdkn1a
- **G1/S:** Cdc6, Cdt1, E2f1, Gins2, Mcm6, Orc1
- **S:** Ccne1, Ccne2, Dhfr, Pcna, Rfc4, Rpa1, Rrm2, Tk1
- **G2:** Anln, Ccnf, Cdk1, Kpna2, Melk, Pbk, Smc2, Top2a
- **G2/M:** Aurka, Bub1b, Ccna2, Ccnb1, Cenpe, Nusap1, Plk1, Tpx2
- **M:** Ccnb2, Cdc20, Cdkn3

**Validation:**
- Periodic LOESS regression of marker expression along θ
- Phase boundaries:
  - 0.5π: Start of S phase
  - 1.0π: Start of G2/M
  - 1.5π: Middle of M phase
  - 1.75π: Start of G1

**Key Calculations:**
- `tricyclePosition`: Cell cycle position (0-2π radians)
- `pi`: Normalized position (0-2, where 2 = 2π)

---

## Dependencies

| Tool | Version | Purpose |
|:-----|:--------|:--------|
| [R](https://www.r-project.org/) | ≥4.0 | Statistical computing |
| [tricycle](https://bioconductor.org/packages/tricycle/) | ≥1.0 | Cell cycle inference |
| [Seurat](https://satijalab.org/seurat/) | ≥4.0 | Single-cell analysis |
| [SeuratWrappers](https://github.com/satijalab/seurat-wrappers) | ≥0.3 | tricycle integration |
| [ggplot2](https://ggplot2.tidyverse.org/) | ≥3.4 | Visualization |

---

## Expected Output

**Metadata Additions:**
- `tricyclePosition` - Cell cycle position (0-2π radians)
- `pi` - Normalized position (0-2)
- `tricycleEmbedding` - Dimensionality reduction coordinates

**Analysis Files:**
- `{project}_tricycle_genes_to_plot_dataframe.csv` - Expression matrix with θ and marker genes

**Visualizations:**
- Periodic LOESS fit plots per marker gene
- Faceted expression profiles across cell cycle (36 markers)
- Phase boundaries annotated with vertical lines

---

## Quick Start

```R
# Load required packages
library(Seurat)
library(SeuratWrappers)
library(tricycle)

# Load integrated Seurat object
seurat_obj <- readRDS("path/to/integrated_object.rds")

# Run tricycle
seurat_obj <- Runtricycle(
  object = seurat_obj,
  slot = "data",
  gname.type = "SYMBOL",
  species = "mouse"
)

# Add normalized position
seurat_obj$pi <- seurat_obj$tricyclePosition / 3.14

# Extract marker expression
markers <- FetchData(
  seurat_obj,
  vars = c("tricyclePosition", "pi", "Top2a", "Ccna2", "Mcm6")
)

# Visualize periodic expression
fit_top2a <- fit_periodic_loess(
  seurat_obj$tricyclePosition,
  seurat_obj[["RNA"]]$data["Top2a", ],
  plot = TRUE
)
```

---

## Key Features

- **Continuous cell cycle states:** Unlike discrete phase assignment, tricycle provides continuous θ values
- **Transfer learning:** Leverages pre-trained reference for accurate inference
- **Species-specific:** Trained on mouse reference (customizable for human)
- **Validation framework:** Built-in periodic regression for marker validation
