# Differential Gene Expression

## Overview

Identifies differentially expressed genes between perturbations and non-targeting controls, followed by Gene Ontology enrichment analysis.

**Input:**
- Integrated Seurat object with gRNA/gene annotations (from data_integration)
- Mixscape classifications (KO/NP/NTC)

**Output:**
- Differential expression tables per perturbation
- GO enrichment results with treeplots

---

## Scripts

| Script | Description | Parameters |
|:-------|:------------|:-----------|
| [`find_markers_gene_ontology.R`](./find_markers_gene_ontology.R) | Performs FindMarkers for each perturbation vs. NT controls and enrichGO analysis | p_adj < 0.05, \|log2FC\| > 0.25, pct.1 > 0.05, pct.2 > 0.01 |

---

## Analysis Details

**Differential Expression (FindMarkers):**
- Assay: SCT (sctransform normalized)
- Comparison: Each perturbation vs. non-targeting (NT) controls
- Filtering thresholds:
  - Adjusted p-value < 0.05
  - \|Average log2 fold change\| > 0.25
  - Expression in ident.1 (pct.1) > 5%
  - Expression in ident.2 (pct.2) > 1%

**Gene Ontology Enrichment (clusterProfiler):**
- Database: org.Mm.eg.db (mouse)
- Ontology: Biological Process (BP)
- Adjustment method: Benjamini-Hochberg (BH)
- p-value cutoff: 0.05
- Similarity metric: Wang semantic similarity
- Visualization: Treeplots with consistent scaling across comparisons

**Databases Queried (enrichR):**
- GO_Biological_Process_2023
- GO_Cellular_Component_2018
- GO_Molecular_Function_2018
- KEGG_2019_Mouse
- Mouse_Gene_Atlas
- WikiPathways_2019_Mouse

---

## Dependencies

| Tool | Version | Purpose |
|:-----|:--------|:--------|
| [R](https://www.r-project.org/) | ≥4.0 | Statistical computing |
| [Seurat](https://satijalab.org/seurat/) | ≥4.0 | Single-cell DE analysis |
| [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/) | ≥4.0 | GO enrichment |
| [enrichplot](https://bioconductor.org/packages/enrichplot/) | ≥1.12 | Enrichment visualization |
| [org.Mm.eg.db](https://bioconductor.org/packages/org.Mm.eg.db/) | ≥3.13 | Mouse gene annotations |
| [GOSemSim](https://bioconductor.org/packages/GOSemSim/) | ≥2.18 | GO semantic similarity |

---

## Expected Output

**Differential Expression Results:**
- `{project}_{GENE}_NT_findmarkers_no_filtering.csv` - Unfiltered DE results per gene vs. NT

**GO Enrichment Results:**
- `{project}_enrichGO_{Gene}_v_NT.rds` - Enriched GO terms per perturbation
- Treeplots showing hierarchical GO term relationships with standardized scaling

**Key Columns in FindMarkers Output:**
- `p_val` - Unadjusted p-value
- `avg_log2FC` - Average log2 fold change
- `pct.1` - % cells expressing gene in perturbation
- `pct.2` - % cells expressing gene in controls
- `p_val_adj` - Adjusted p-value (Bonferroni)

---

## Quick Start

```R
# Load required packages
library(Seurat)
library(clusterProfiler)

# Load integrated object with mixscape classifications
seurat_obj <- readRDS("path/to/mixscape_object.rds")

# Run FindMarkers (example: Opa1 vs. NT)
DefaultAssay(seurat_obj) <- "SCT"
seurat_obj <- PrepSCTFindMarkers(seurat_obj)
Idents(seurat_obj) <- seurat_obj$gene

opa1_markers <- FindMarkers(seurat_obj, ident.1 = "Opa1", ident.2 = "NT")

# Run GO enrichment
opa1_enrichGO <- enrichGO(
  gene = opa1_markers$gene,
  OrgDb = "org.Mm.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  keyType = "SYMBOL"
)

# Visualize
treeplot(pairwise_termsim(opa1_enrichGO))
```
