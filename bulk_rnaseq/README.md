# Bulk RNA-seq Analysis

## Overview

Bulk RNA-seq validation of single-gRNA CRISPR perturbations. RFP-positive transduced cells were sorted at day 6 post-transduction, RNA extracted, and libraries sequenced on the NovaSeq X (50 bp paired-end).

**Input:**
- Paired-end FASTQ files (50 bp, NovaSeq X)

**Output:**
- Gene-level read counts per sample
- Differential expression results (edgeR)

> **Note:** No scripts are provided for this module. The commands below document the exact tools and parameters used.

---

## Workflow

### 1. Adapter Trimming — Trimmomatic v0.39

Paired-end trimming with TruSeq3 adapters:

```bash
trimmomatic PE \
  -threads 24 \
  Sample_Read_1.fastq Sample_Read_2.fastq \
  Sample_Read_1_trimmed_paired.fastq \
  Sample_Read_1_trimmed_unpaired.fastq \
  Sample_Read_2_trimmed_paired.fastq \
  Sample_Read_2_trimmed_unpaired.fastq \
  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True \
  LEADING:3 TRAILING:3 MINLEN:36
```

Adapter sequences: [TruSeq3-PE-2.fa](https://github.com/usadellab/Trimmomatic/tree/main/adapters)

### 2. Alignment — RUM v2.0.4

Align trimmed reads to mm10:

```bash
rum_runner align \
  --index-dir ../rum_indexes/mm10 \
  --output Sample_trimmed_aligned \
  --name Sample_trimmed \
  --chunks 32 \
  --variable-length-reads \
  Sample_Read_1_trimmed_paired.fastq Sample_Read_2_trimmed_paired.fastq
```

### 3. Read Counting — HTSeq-count v2.0.3

Count reads in union mode:

```bash
htseq-count \
  -f bam \
  -r name \
  -m union \
  --additional-attr gene_name \
  -o Sample_counts.sam \
  -s no \
  Sample_trimmed_aligned/RUM.bam \
  Reference_genes.gtf > Sample_counts.txt
```

### 4. Differential Expression — edgeR v4.4.2

Performed in R following the [edgeR user guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf).

---

## Dependencies

| Tool | Version | Purpose |
|:-----|:--------|:--------|
| [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) | v0.39 | Adapter trimming |
| [RUM](https://github.com/itmat/rum) | v2.0.4 | Alignment (mm10) |
| [HTSeq-count](https://htseq.readthedocs.io/) | v2.0.3 | Read counting |
| [edgeR](https://bioconductor.org/packages/edgeR/) | v4.4.2 | Differential expression |

---

## Expected Output

- `Sample_counts.txt` — Gene-level read counts per sample
- edgeR differential expression tables per perturbation vs. NT control