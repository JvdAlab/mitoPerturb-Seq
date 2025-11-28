# Heteroplasmy Variance Analysis

## Overview

Comprehensive statistical analysis of heteroplasmy variance across genetic perturbations, accounting for technical measurement variance and biological variability.

**Input:**
- Integrated metadata CSV with heteroplasmy calls and mtDNA depth
- Mixscape classifications (KO/NP/NTC)
- Measurement variance function (from depth simulations)

**Output:**
- Statistical models quantifying heteroplasmy variance per perturbation
- Weighted linear regression accounting for technical variance
- Interactive HTML report with figures

---

## Scripts

| Script | Description | Parameters |
|:-------|:------------|:-----------|
| [`MitoPerturbSeq_heteroplasmy_analysis.qmd`](./src/MitoPerturbSeq_heteroplasmy_analysis.qmd) | Quarto document performing heteroplasmy variance analysis with weighted linear models | Depth threshold: 20×, measurement variance weighting |

---

## Analysis Details

**Data Requirements:**
- `integrated_dataset_metadata.csv` - Cell metadata with heteroplasmy and depth
- `integrated_mixscape_metadata.csv` - Perturbation classifications
- `het_depth_sim.rds` - Heteroplasmy-depth simulation results
- `measurement_variance_function.rds` - Technical variance predictor

**Statistical Approach:**
- Weighted linear regression using empirical measurement variance
- Variance partitioning: biological vs. technical
- Quality threshold: mtDNA depth ≥ 20×

**Perturbations Analyzed:**
- Mitochondrial biogenesis: Akap1, Nnt, Polg, Tfam
- Mitochondrial dynamics: Dnm1l, Mtfp1, Prkn, Mfn1, Mfn2, Opa1
- Mitophagy: Pink1, Bnip3l
- Non-targeting controls (NT)

---

## Dependencies

| Tool | Version | Purpose |
|:-----|:--------|:--------|
| [R](https://www.r-project.org/) | ≥4.0 | Statistical computing |
| [Quarto](https://quarto.org/) | ≥1.3 | Reproducible reporting |
| [tidyverse](https://www.tidyverse.org/) | ≥1.3 | Data manipulation |
| [ggplot2](https://ggplot2.tidyverse.org/) | ≥3.4 | Visualization |
| [broom](https://broom.tidymodels.org/) | ≥1.0 | Model summaries |
| [car](https://cran.r-project.org/package=car) | ≥3.0 | Statistical tests |

---

## Expected Output

**HTML Report:**
- `MitoPerturbSeq_heteroplasmy_analysis.html` - Interactive analysis with embedded figures

**Key Results:**
- Heteroplasmy variance estimates per perturbation (weighted)
- Statistical significance tests (ANOVA, post-hoc comparisons)
- Cell count summaries passing depth threshold
- Diagnostic plots for model assumptions

---

## Quick Start

```bash
# Render Quarto document
quarto render src/MitoPerturbSeq_heteroplasmy_analysis.qmd

# Output: src/MitoPerturbSeq_heteroplasmy_analysis.html
```
