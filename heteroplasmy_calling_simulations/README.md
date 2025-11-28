# Heteroplasmy Calling Simulations

## Overview

Simulates heteroplasmy sampling and quantifies SNV correlation accuracy at varying sequencing depths using the MitoPerturb-Seq dataset.

**Input:**
- Seurat object with mgatk mitochondrial variant calls
- Simulated heteroplasmy distributions (normal: mean=0.58, sd=0.11)

**Output:**
- Heteroplasmy sampling distributions at specified depths
- SNV correlation matrices across depth thresholds

---

## Scripts

| Script | Description | Parameters |
|:-------|:------------|:-----------|
| [`Heteroplasmy_Simulations.R`](./Heteroplasmy_Simulations.R) | Simulates heteroplasmy sampling from mtDNA populations at varying sequencing depths | Copy number: 1750×7, depth cutoff: 20× |
| [`Heteroplasmic_SNV_Correlations.R`](./Heteroplasmic_SNV_Correlations.R) | Calculates pairwise correlations between m.5024C>T and other heteroplasmic SNVs across depth thresholds (0-150×) | Depth increments: 5×, reference SNV: m.5024C>T |

---

## Analysis Details

**Heteroplasmy Simulation:**
- mtDNA copy number: 1750 genomes × 7 SNV sites = 12,250 sites
- Population mean heteroplasmy: 0.58
- Population standard deviation: 0.11
- Minimum depth threshold: 20×

**SNV Correlation Analysis:**
- **Reference SNV:** m.5024C>T
- **Compared SNVs:** m.13715C>T, m.13614C>T, m.1781C>T, m.1866A>G, m.3009G>T, m.3823C>T
- **Depth range:** 0× to 150× (5× increments)
- **Outputs:** Pearson correlation coefficients and % cells retained at each threshold

---

## Dependencies

| Tool | Version | Purpose |
|:-----|:--------|:--------|
| [R](https://www.r-project.org/) | ≥4.0 | Statistical computing |
| [Seurat](https://satijalab.org/seurat/) | ≥4.0 | Single-cell analysis |
| [tidyverse](https://www.tidyverse.org/) | ≥1.3 | Data manipulation |

---

## Expected Output

**Simulation Results:**
- Heteroplasmy call distributions at specified sequencing depths
- Comparison of true vs. sampled heteroplasmy values

**SNV Correlation Results:**
- Correlation matrix: m.5024C>T vs. 6 other SNVs across depth thresholds
- Cell retention percentage at each depth cutoff
