# mitoPerturb-Seq: Multi-Modal Analysis of Mitochondrial Perturbations

This repository contains the complete analysis code for the mitoPerturb-Seq manuscript, a comprehensive study combining CRISPR perturbation screening, single-cell RNA sequencing, and chromatin profiling to investigate mitochondrial function and stress response.

## Overview

mitoPerturb-Seq integrates multiple experimental and computational approaches:

- **CRISPR guide detection and classification** - Identifying perturbed cells using mixscape
- **Single-cell RNA-seq integration** - Multi-sample integration with perturbation annotations
- **Gene regulatory network inference** - pySCENIC analysis of transcription factor activity
- **Cell cycle profiling** - Tricycle-based continuous cell cycle scoring
- **Differential gene expression** - Perturbation-specific transcriptional responses
- **Mitochondrial heteroplasmy analysis** - Measuring and modeling mtDNA variant frequencies
- **Chromatin profiling** - DamID-seq analysis of ATF4 binding sites

## Repository Structure

```
mitoPerturb-Seq/
├── guide_detection/              # CRISPR guide detection and mixscape (Steps 1-4)
├── data_integration/             # Single-cell dataset integration
├── pyscenic/                     # Gene regulatory network inference
├── tricycle/                     # Cell cycle phase annotation
├── differential_gene_expression/ # DGE and GO enrichment analysis
├── heteroplasmy_variance_analysis/ # Mitochondrial heteroplasmy modeling
├── DamID-seq-peak-calling/       # ATF4 chromatin binding analysis
├── Heteroplasmic_SNV_Correlations.R  # SNV correlation analysis
└── Heteroplasmy_Simulations.R    # Sequencing depth simulations
```

## Quick Start

### Prerequisites

- **R** (≥4.0) with Seurat, Signac, and Bioconductor packages
- **Python** (≥3.8) for pySCENIC
- **Command-line tools**: BWA, UMI-tools, SAMtools, Bowtie2, MACS3, HOMER
- **HPC environment**: Most analyses require SLURM job submission

### Installation

1. Clone this repository:
```bash
git clone https://github.com/[username]/mitoPerturb-Seq.git
cd mitoPerturb-Seq
```

2. Install dependencies (see individual module READMEs for specific requirements):
```bash
# For R packages, install from individual R scripts
# For DamID-seq pipeline
cd DamID-seq-peak-calling
micromamba env create -f environment.yml
micromamba activate damid_pipeline
```

3. Configure paths and parameters as needed in each module

## Analysis Workflow

### Core scRNA-seq Analysis

```
Raw FASTQ files
    ↓
1. Guide Detection (guide_detection/)
    - Extract cell barcodes with UMI-tools
    - Align guide sequences with BWA
    - Detect guides and classify with mixscape
    ↓
2. Data Integration (data_integration/)
    - Integrate Seurat objects across samples
    - Add guide detection and classification results
    ↓
3. Downstream Analyses (parallel)
    ├─→ pySCENIC (pyscenic/)
    │   - GRN inference and TF activity scoring
    ├─→ Tricycle (tricycle/)
    │   - Continuous cell cycle annotation
    └─→ Differential Expression (differential_gene_expression/)
        - FindMarkers and GO enrichment
```

### Parallel Analysis Tracks

**Mitochondrial Heteroplasmy** (heteroplasmy_variance_analysis/)
- Independent analysis of mtDNA variant frequencies
- Variance modeling and quality control
- Uses integrated dataset metadata

**ATF4 Chromatin Binding** (DamID-seq-peak-calling/)
- Independent DamID-seq pipeline
- Comprehensive peak calling with reproducibility filtering
- Motif enrichment and gene annotation

## Module Descriptions

### guide_detection/
Detects CRISPR guide RNAs in single cells and classifies perturbations using Seurat's mixscape framework.

**Key outputs:** Guide assignments per cell, perturbation classifications

**Runtime:** ~2-4 hours per sample

**README:** [guide_detection/README.md](guide_detection/README.md)

### data_integration/
Integrates multiple single-cell RNA-seq samples with guide detection results to create a unified dataset.

**Key outputs:** Integrated Seurat object with perturbation annotations

**Dependencies:** Requires guide_detection (steps 1-3)

**README:** [data_integration/README.md](data_integration/README.md)

### pyscenic/
Performs gene regulatory network (GRN) inference using pySCENIC to identify transcription factor regulons and their activity.

**Key outputs:** TF regulon activities per cell, regulon gene sets

**Runtime:** ~6-12 hours (GRNBoost2 is rate-limiting)

**README:** [pyscenic/README.md](pyscenic/README.md)

### tricycle/
Applies tricycle to project cells onto a continuous cell cycle space without discrete phase assignments.

**Key outputs:** Cell cycle coordinates (theta), pseudotime

**README:** [tricycle/README.md](tricycle/README.md)

### differential_gene_expression/
Performs differential gene expression analysis between perturbation groups with Gene Ontology enrichment.

**Key outputs:** Marker genes per perturbation, GO enrichment results

**README:** [differential_gene_expression/README.md](differential_gene_expression/README.md)

### heteroplasmy_variance_analysis/
Quantifies mitochondrial DNA heteroplasmy levels and models measurement variance as a function of sequencing depth.

**Key outputs:** Heteroplasmy estimates, variance models, QC metrics

**Key analysis:** Quarto document in src/

**Data:** Includes simulation results and metadata

### DamID-seq-peak-calling/
Complete pipeline for processing DamID-seq data to identify ATF4 chromatin binding sites with reproducibility filtering.

**Key outputs:** 5,789 reproducible ATF4 peaks, gene annotations, enriched motifs

**Status:** Most mature module with comprehensive documentation (already well-organized internally with numbered subdirectories)

**README:** [DamID-seq-peak-calling/README.md](DamID-seq-peak-calling/README.md)

## Key Dependencies

### R Packages
- **Single-cell:** Seurat, Signac, SingleCellExperiment, SeuratWrappers
- **Genomics:** GenomicRanges, BSgenome.Mmusculus.UCSC.mm10, EnsDb.Mmusculus.v79
- **Analysis:** ComplexHeatmap, clusterProfiler, fgsea, DESeq2
- **Specialized:** tricycle, loomR, presto, gprofiler2

### Python Tools
- **pySCENIC:** arboreto, pyscenic (for GRNBoost2, CTX pruning, AUCell)
- **General:** pandas, numpy, scanpy, loompy

### Command-Line Tools
- **Alignment:** BWA, Bowtie2, SAMtools
- **Peak calling:** MACS3, IDR, bedtools
- **Annotation:** HOMER
- **Preprocessing:** UMI-tools, Trimmomatic
- **Visualization:** deepTools, pyGenomeTracks

See individual module READMEs for specific version requirements.

## Data Availability

**Raw sequencing data:** [GEO accession to be added]

**Processed data:** [Zenodo/Figshare link to be added]

**Reference genomes:**
- Mouse genome: mm10 (GRCm38)
- Mitochondrial genome: NC_005089.1

## Computational Requirements

- **HPC cluster** with SLURM job scheduler recommended
- **Memory:** 32-128 GB RAM depending on module
- **Storage:** ~500 GB for raw data, intermediates, and outputs
- **Runtime:** Complete workflow ~24-48 hours walltime

See individual module READMEs for specific requirements.

## Citation

If you use this code, please cite:

```
[Author list]. (Year). [Manuscript title].
[Journal]. doi: [DOI]
```

For individual tools used in this analysis, please cite:
- **Seurat:** Hao et al., Cell 2021
- **Signac:** Stuart et al., Nat Methods 2021
- **pySCENIC:** Van de Sande et al., Nat Protoc 2020
- **mixscape:** Papalexi et al., Nat Genet 2021
- **tricycle:** Zheng et al., Genome Biol 2022
- **MACS3:** Zhang et al., Genome Biol 2008
- **HOMER:** Heinz et al., Mol Cell 2010

## License

[License type to be added - e.g., MIT, GPL-3]

## Contact

For questions about the code or analysis:
- [Primary contact name and email]
- Open an issue on GitHub

For questions about the manuscript:
- [Corresponding author name and email]

## Repository Organization

This repository is organized for maximum reproducibility:
- Each module contains a detailed README with usage instructions
- Scripts are numbered to indicate execution order
- Configuration files specify all parameters
- Environment files document exact software versions (where applicable)

For a complete description of the repository organization and development plan, see [REPOSITORY_ORGANIZATION_PLAN.md](REPOSITORY_ORGANIZATION_PLAN.md).

## Troubleshooting

### Common Issues

**Issue:** R package installation fails
- **Solution:** Check R version (≥4.0 required). Install Bioconductor packages first: `BiocManager::install(c("GenomicRanges", "BSgenome.Mmusculus.UCSC.mm10"))`

**Issue:** pySCENIC fails with memory error
- **Solution:** Increase memory allocation in SLURM script. GRNBoost2 may need 64-128 GB for large datasets.

**Issue:** Guide detection finds no guides
- **Solution:** Verify guide reference FASTA matches your library. Check that FASTQ files contain guide sequences (usually in R2).

**Issue:** DamID-seq pipeline fails at alignment
- **Solution:** Ensure Bowtie2 index was built successfully. Check that GATC fragment file was generated for mm10.

For module-specific troubleshooting, see individual READMEs.

## Acknowledgments

This work was supported by [funding information].

We thank [acknowledgments].

## Version History

- **v1.0.0** (2025-XX-XX): Initial release with manuscript publication
- Development version: Code refactoring and documentation

---

**Last updated:** 2025-11-26
**Maintainer:** [Name]
**Status:** Active development
