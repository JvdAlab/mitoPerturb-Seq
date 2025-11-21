##### Find markers
##### Gene ontology

# load packages
library(clusterProfiler)
library(ComplexHeatmap)
library(cowplot)
library(DESeq2)
library(dplyr)
library(enrichplot)
library(enrichR)
library(future)
library(genomation)
library(GenomicRanges)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(glmGamPoi)
library(GOSemSim)
library(gprofiler2)
library(irlba)
library(matrixStats)
library(mixtools) 
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(readr)
library(reshape2)
library(R.utils)
library(scales)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(Signac)
library(SingleCellExperiment)
library(SingleR)
library(stringr)
library(tidyverse)
library(viridis)

library(AnnotationDbi)
library(AnnotationHub)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(msigdbr)
library(org.Mm.eg.db)
library(biovizBase)

#### Set Project Variables ####

# enter project name
project <- "project"

# set base directory
baseDir <- "/../../"
setwd(baseDir)


##### STEP ONE - Load integrated data #####

# set integrated rds results folder
integrate <- "/../../mixscape_dir/"

# load the Seurat object following integration and addition of guide annotation
# loaded the mixscape object here as this contains additional columns such as guide/gene and mixscape classifications
MultiOme_obj_no_unkmul <- readRDS(paste0(integrate,"project_no_unkmul_polg6_RNA_ATAC_MGATK_GUIDES_mixscape.rds"))
MultiOme_obj_no_unkmul

# set results directory
markers_folder <- "/../../markers/"

# list databases
enrichR_DB <- as.data.frame(listEnrichrDbs())

# select databases of interest - these will be used to query data
db <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018", 
        "KEGG_2019_Mouse", "Mouse_Gene_Atlas", "WikiPathways_2019_Mouse")   

# check if the website is live
websiteLive <- TRUE

d <- godata("org.Mm.eg.db", ont = "BP")


##### STEP TWO - FindMarkers #####

# NB: ran for all genes, one shown as an example

# set default assay
DefaultAssay(MultiOme_obj_no_unkmul) <- "SCT"

# run prep sct find markers 
MultiOme_obj_no_unkmul <- PrepSCTFindMarkers(MultiOme_obj_no_unkmul, assay = "SCT", verbose = TRUE)

# set idents - gene column
Idents(MultiOme_obj_no_unkmul) <- MultiOme_obj_no_unkmul@meta.data$gene

# genes v NT
opa1 <- FindMarkers(MultiOme_obj_no_unkmul, ident.1 = "Opa1", ident.2 = "NT")
write.csv(opa1, paste(markers_folder, project, "_OPA1_NT_findmarkers_no_filtering.csv", sep = ""))


##### STEP THREE - clusterprofiler ##### 

# NB: ran for all genes, one shown as an example

# load data
opa1_data <- read_csv("/../../markers/project_OPA1_NT_findmarkers_no_filtering.csv")

# look at data
head(opa1_data)
# create gene column
opa1_data$gene <- opa1_data$...1
# convert first column to rownames
opa1_data <- opa1_data %>% column_to_rownames(var = "...1")

# extract only those that pass p_val_adj, avg_log2FC and pct thresholds
opa1_data_plog <- opa1_data[opa1_data$p_val_adj < 0.05 & abs(opa1_data$avg_log2FC) > 0.25 &
                              opa1_data$pct.1 > 0.05 & opa1_data$pct.2 > 0.01,]

### clusterprofiler enrichGO

# create gene list
opa1_CP_genes <- opa1_data_plog$gene
# run CP enrichGO
opa1_CP_enrichGO <- enrichGO(gene = opa1_CP_genes, 
                             OrgDb = "org.Mm.eg.db", 
                             ont = "BP", 
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             keyType = "SYMBOL", 
                             readable=TRUE)

opa1_enrichGO <- data.frame(opa1_CP_enrichGO)
saveRDS(opa1_enrichGO, file = "/../../markers/project_enrichGO_Opa1_v_NT.rds")

opa1_CP_enrichGO2 <- pairwise_termsim(opa1_CP_enrichGO, method="Wang", semData = d)

treeplot(opa1_CP_enrichGO2)


##### STEP FOUR - Combining results #####

# generate multiple tree plots with consistent scaling
# example with three genes
ego_list <- list(opa1_CP_enrichGO, polg_CP_enrichGO, tfam_CP_enrichGO)

# generate pairwise similarity for each enrichGO object
pairwise_sim_list <- lapply(ego_list, pairwise_termsim, method = "Wang", semData = d)

# extract gene count range
gene_count_range <- range(sapply(pairwise_sim_list, function(x) x@result$Count), na.rm = TRUE)

# create treeplots with consistent scaling
plots_coloured <- lapply(pairwise_sim_list, function(sim_result) {
  treeplot(sim_result) +
    scale_size_continuous(limits = gene_count_range, name = "Gene Count") +  # Standardize size
    theme_minimal()  + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank()) # remove background lines and axis ticks
})

# plots side-by-side
plots_coloured[[1]] + plots_coloured[[2]] + plots_coloured[[3]]

