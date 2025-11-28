# pySCENIC Gene Regulatory Network Analysis

## Overview

Infers gene regulatory networks (GRNs) and transcription factor (TF) activities from single-cell RNA-seq data using the pySCENIC workflow.

**Input:**
- Integrated Seurat object with RNA assay
- Human TF list and motif databases (hg38)

**Output:**
- Gene-TF adjacency matrices
- Regulon activity scores (AUCell) per cell
- Loom file with regulon annotations

---

## Workflow

| Script | Description | Parameters |
|:-------|:------------|:-----------|
| [`01_prep_pyscenic.R`](./01_prep_pyscenic.R) | Converts Seurat object to loom format, filters genes present in pySCENIC databases | RNA assay, gene filtering |
| [`02_pyscenic_arboreto_ctx_aucell.sh`](./02_pyscenic_arboreto_ctx_aucell.sh) | Runs pySCENIC pipeline: GRNBoost2, cisTarget motif enrichment, AUCell scoring | Method: grnboost2, seed: 777, 32 workers |
| [`03_post_pyscenic_add_metadata.py`](./03_post_pyscenic_add_metadata.py) | Adds Seurat metadata to pySCENIC output loom file | Python h5py |
| [`04_pyscenic_analysis.R`](./04_pyscenic_analysis.R) | Loads loom results, extracts regulons and AUCell scores for downstream analysis | AUCell thresholding |

---

## Analysis Details

**Step 1: GRN Inference (GRNBoost2)**
- Algorithm: Gradient boosting regression trees
- Input: Expression matrix + TF list (hg38)
- Output: Gene-TF adjacency matrix with importance scores

**Step 2: Regulon Pruning (cisTarget)**
- Databases:
  - `hg38_10kbp_up_10kbp_down_full_tx_v10_clust`
  - `hg38_500bp_up_100bp_down_full_tx_v10_clust`
- Motif annotations: `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`
- Filters: TF-target links with motif support in promoter regions

**Step 3: Activity Scoring (AUCell)**
- Calculates regulon enrichment scores per cell
- Output: AUCell matrix (regulons × cells)

**Step 4: Integration**
- Adds Seurat metadata (clusters, gRNA, heteroplasmy) to loom file
- Enables regulon analysis across perturbations

---

## Dependencies

| Tool | Version | Purpose |
|:-----|:--------|:--------|
| [Python](https://www.python.org/) | ≥3.7 | pySCENIC pipeline |
| [pySCENIC](https://pyscenic.readthedocs.io/) | ≥0.12 | GRN inference |
| [arboreto](https://arboreto.readthedocs.io/) | ≥0.1.6 | GRNBoost2 algorithm |
| [loompy](http://loompy.org/) | ≥3.0 | Loom file handling |
| [R](https://www.r-project.org/) | ≥4.0 | Pre/post-processing |
| [Seurat](https://satijalab.org/seurat/) | ≥4.0 | Data conversion |
| [loomR](https://github.com/mojaveazure/loomR) | ≥0.2 | Loom import to R |
| [SCopeLoomR](https://github.com/aertslab/SCopeLoomR) | ≥0.13 | SCENIC loom utilities |

**Reference Databases:**
- Human TF list: `allTFs_hg38.txt`
- Motif rankings: `hg38_*_clust.genes_vs_motifs.rankings.feather`
- Motif scores: `hg38_*_clust.genes_vs_motifs.scores.feather`

---

## Expected Output

**pySCENIC Results:**
- `{project}_rna_assay_filtered.loom` - Expression matrix in loom format
- `{project}_rna_assay_filtered_adjacencies.csv` - Gene-TF co-expression network
- `{project}_rna_assay_filtered_ctxreg.csv` - Pruned regulons with motif support
- `{project}_rna_assay_pyscenic_output_filtered.loom` - Final loom with AUCell scores and metadata

**Regulon Data:**
- Regulon-gene incidence matrix
- AUCell enrichment scores per cell
- TF activity profiles across perturbations

---

## Quick Start

```bash
# Step 1: Prepare loom file from Seurat object
Rscript 01_prep_pyscenic.R

# Step 2: Run pySCENIC pipeline (HPC/SLURM)
sbatch 02_pyscenic_arboreto_ctx_aucell.sh

# Step 3: Add metadata to output loom
python 03_post_pyscenic_add_metadata.py

# Step 4: Analyze regulons in R
Rscript 04_pyscenic_analysis.R
```

---

## Notes

- Bash script assumes SLURM workload manager
- Runtime: ~24-96 hours depending on dataset size and available cores
- Requires 64+ GB RAM for large datasets
- Uses hg38 databases (human); modify for mouse (mm10) if needed
