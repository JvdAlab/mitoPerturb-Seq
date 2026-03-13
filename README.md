# mitoPerturb-Seq

[![DOI](https://zenodo.org/badge/1088389701.svg)](https://doi.org/10.5281/zenodo.19008664)

Analysis code for:

> **"MitoPerturb-Seq identifies gene-specific single-cell responses to mitochondrial DNA depletion and heteroplasmy**
> [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2025.07.08.663208v1.full)

---

## Overview

Computational pipelines for the MitoPerturb-Seq study: single-cell multiome (scRNA-seq + scATAC-seq) with CRISPR perturbation screening across 13 mitochondrial target genes, mitochondrial heteroplasmy calling and variance modelling, and ATF4 chromatin binding analysis using DamID-seq. ~6,500 mouse embryonic fibroblasts (MEFs) from the m.5024C>T heteroplasmic mouse line.

---

## Repository Structure

```
mitoPerturb-Seq/
├── gRNA_assignment/                  # UMI-tools barcode extraction → BWA alignment → mixscape classification
├── data_integration/                 # RNA + ATAC + mgatk multimodal integration (Seurat/Signac)
├── pyscenic/                         # GRN inference (GRNBoost2 → cisTarget → AUCell)
├── tricycle/                         # Continuous cell cycle annotation
├── differential_gene_expression/     # FindMarkers per perturbation + GO enrichment
├── heteroplasmy_variance_analysis/   # Depth-dependent variance modelling (Quarto report + data)
├── heteroplasmy_calling_simulations/ # In silico depth simulations and SNV correlation analysis
├── DamID-seq-peak-calling/           # ATF4 binding: Bowtie2 → MACS3 → IDR → HOMER
└── bulk_rnaseq/                      # Bulk RNA-seq validation (Trimmomatic → RUM → HTSeq → edgeR)
```

See each module's `README.md` for full parameters, dependencies, and usage.

---

## Requirements

R ≥ 4.3.3, Python 3.12.2, and an HPC environment with SLURM. See each module's `README.md` for the full list of packages and command-line tools.

---

## Data Availability

Raw sequencing data are deposited in GEO:

| Accession | Contents |
|:----------|:---------|
| [GSE297416](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297416) | 10X Multiome ATAC + Gene Expression; gRNA enrichment (TAP-seq); mtDNA enrichment (xGEN Hybrid Capture) |
| [GSE297418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297418) | DamID-seq (Dam-only and Dam-ATF4) |
| [GSE297491](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297491) | Bulk RNA-seq validation |

---

## Citation

If you use this code, please cite:

```
Burr et. al. "MitoPerturb-Seq identifies gene-specific single-cell responses to
mitochondrial DNA depletion and heteroplasmy."

Preprint: bioRxiv (2025). https://doi.org/10.1101/2025.07.08.663208
```

---

## License

This code is released under the [MIT License](LICENSE).

---

## Contact

- Corresponding Author: [Jelle van den Ameele](mailto:jv361@cam.ac.uk)
- Maintainer: [Abhilesh Dhawanjewar](mailto:ad2347@cam.ac.uk)
