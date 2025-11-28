# gRNA Assignment

## Overview

Assigns CRISPR guide RNAs to single cells from 10X Chromium FASTQ files and classifies perturbation status using Seurat's mixscape framework.

**Input:**
- Raw FASTQ files (R1: cell barcode/UMI, R2: guide sequences)
- Seurat object with scRNA-seq data
- Guide library FASTA reference

**Output:**
- Seurat object with guide assignments and mixscape classifications
- Guide count matrices and QC plots

---

## Workflow

| Script | Description | Parameters |
|:-------|:------------|:-----------|
| **01: Barcode Extraction** |
| [`01_umitools_extract_cellbarcode.sh`](./01_umitools_extract_cellbarcode.sh) | Extracts 16-bp cell barcodes from R1 and appends them to R2 read headers | [UMI-tools](https://github.com/CGATOxford/UMI-tools) ≥1.0.0, `--bc-pattern=NNNNNNNNNNNNNNNN` |
| **02: Guide Alignment** |
| [`02_bwa_align_guides_cellbarcode.sh`](./02_bwa_align_guides_cellbarcode.sh) | Aligns R2 reads containing guide sequences to guide library reference | [BWA](https://github.com/lh3/bwa) ≥0.7.17, [SAMtools](https://github.com/samtools/samtools) ≥1.10 |
| **03: Guide Assignment** |
| [`03_guide_detection.R`](./03_guide_detection.R) | Counts guide UMIs per cell and assigns primary guide based on threshold | Min UMI threshold: 2 |
| [`03_guide_detection_multiple_guides.R`](./03_guide_detection_multiple_guides.R) | Handles cells with multiple guide assignments (doublets or multi-guide designs) | Optional: For multi-guide libraries |
| **04: Perturbation Classification** |
| [`04_mixscape_analysis.R`](./04_mixscape_analysis.R) | Classifies cells as KO/NP/NTC using mixscape based on perturbation signatures | Min cells per guide: 50, requires NT controls |

---

## Analysis Details

**Barcode Extraction (UMI-tools):**
- Pattern: 16-bp cell barcode from R1 (10X Chromium v2/v3)
- Appends barcode to R2 read headers for downstream tracking

**Guide Alignment (BWA):**
- Reference: Guide library FASTA (requires BWA index)
- Filters: Aligned reads only (`samtools view -F 4`)
- Output: BAM file with guide alignments per cell barcode

**Guide Assignment:**
- Counts UMIs per guide per cell
- Assignment threshold: Typically ≥2 UMIs
- Handles single or multiple guide assignments

**Mixscape Classification:**
- **KO (Knock-out):** Cells showing perturbation signature
- **NP (Non-perturbed):** Cells with guide but no transcriptional response
- **NTC (Non-targeting control):** Cells with control guides
- Requires: ≥50 cells per guide, non-targeting control guides in library

---

## Dependencies

| Tool | Version | Purpose |
|:-----|:--------|:--------|
| [UMI-tools](https://github.com/CGATOxford/UMI-tools) | ≥1.0.0 | Barcode extraction |
| [BWA](https://github.com/lh3/bwa) | ≥0.7.17 | Guide alignment |
| [SAMtools](https://github.com/samtools/samtools) | ≥1.10 | BAM processing |
| [R](https://www.r-project.org/) | ≥4.0 | Guide detection & classification |
| [Seurat](https://satijalab.org/seurat/) | ≥4.0 | Single-cell analysis & mixscape |
| [tidyverse](https://www.tidyverse.org/) | ≥1.3 | Data manipulation |

Bash scripts assume execution on an HPC cluster with SLURM workload manager.

---

## Expected Output

**FASTQ Files:**
- `Sample_R2_UMI_extract_*.fastq.gz` - R2 reads with cell barcodes in headers

**Alignment Files:**
- `Sample_bwa_R2_aligned_guides.bam` - Aligned guide reads
- `Sample_bwa_R2_aligned_guides_F4.sam` - Filtered alignments (mapped only)

**Seurat Objects:**
- `Sample_with_guides.rds` - Seurat object with guide assignments
- `Sample_with_mixscape.rds` - Seurat object with KO/NP/NTC classifications

**Guide Metrics:**
- `guide_counts_per_cell.csv` - UMI count matrix (guides × cells)
- QC plots: UMI distribution, guides per cell, assignment statistics
- UMAP plots colored by guide, gene, and perturbation status

---

## Quick Start

```bash
# Step 1: Extract cell barcodes
sbatch 01_umitools_extract_cellbarcode.sh

# Step 2: Align guides (requires BWA index: bwa index guide_library.fasta)
sbatch 02_bwa_align_guides_cellbarcode.sh

# Step 3: Assign guides to cells
Rscript 03_guide_detection.R

# Step 4 (optional): Handle multiple guides per cell
Rscript 03_guide_detection_multiple_guides.R

# Step 5: Classify perturbations with mixscape
Rscript 04_mixscape_analysis.R
```
