##### pySCENIC analysis 
##### Load loom file from step 3

# load packages
library(arrow)
library(AUCell)
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(igraph)
library(loomR)
library(Matrix)
library(RColorBrewer)
library(readr)
library(reshape2)
library(R.utils)
library(scales)
library(SCENIC)
library(SCopeLoomR)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratWrappers)
library(Signac)
library(SingleCellExperiment)
library(tidyverse)


#### Project Variables ####

# enter project name
project <- "project"

# set directory 
baseDir <- "/../../"
setwd(baseDir)

# set seurat directory path
seurat_dir <- paste0(baseDir, "seurat/")

# set pyscenic results directory path 
pyscenic_dir2 <- paste0(seurat_dir, "pyscenic_analyses/")
pyscenic_dir <- paste0(seurat_dir, "pyscenic_analyses/figures/")


##### STEP ONE - Read in metadata #####

# read in metadata, generated during step 1
metadata <- read.csv(paste0(pyscenic_dir2, project, "_metadata.csv")) 
head(metadata)
rownames(metadata) <- metadata$X
metadata$X <- NULL
head(metadata)


##### STEP TWO - Load Loom file #####

# file path to loom file
scenic_loom_path <- file.path("/../../project_rna_assay_pyscenic_output_filtered.loom")

# open loom file 
loom <- open_loom(scenic_loom_path)

### read and extract information from loom file:
# expression matrix
exprMat <- get_dgem(loom)
# normalise expression matrix
exprMat_log <- log2(exprMat+1) 

# regulons
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
# convert the regulons (stored as an incidence matrix) to a list of genes - this function is from SCENIC
regulons <- regulonsToGeneLists(regulons_incidMat)

# AUCell matrix 
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")

# removes regulons that are named as regulon_extended - this function is from SCENIC
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

# AUC thresholds
regulonAucThresholds <- get_regulon_thresholds(loom)

# embeddings (t-SNE, UMAP)
embeddings <- get_embeddings(loom)

# close loom file connection
close_loom(loom)


##### STEP THREE - Average Regulon Activity by Gene #####

# set column metadata column to be used for plotting
GROUPS <- "gene"

# split the cells by gene and calculate average expression
regulonActivity_byGENE <- sapply(split(rownames(metadata), metadata[[GROUPS]]), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

# scale expression
regulonActivity_byGENE_Scaled <- t(scale(t(regulonActivity_byGENE), center = TRUE, scale=TRUE))

# remove NA's
regulonActivity_byGENE_Scaled_remove_NA <- na.omit(regulonActivity_byGENE_Scaled)

# save data 
write.csv(regulonActivity_byGENE_Scaled_remove_NA, paste(pyscenic_dir, project, "_regulonActivity_byGENE_Scaled_remove_NA.csv", sep = ""), row.names = TRUE)

# complex heatmap
f2 = circlize::colorRamp2( c(-2, 0, 2.5), c("blue", "white","red"), space = "RGB")
heatmap_gene <- ComplexHeatmap::Heatmap(regulonActivity_byGENE_Scaled_remove_NA, name="Regulon activity",
                                        col = f2, 
                                        clustering_distance_rows = "pearson", clustering_distance_columns = "pearson", 
                                        cluster_columns = TRUE, show_row_names = TRUE,
                                        row_names_gp=grid::gpar(fontsize=4)) # row font size

pdf(paste(pyscenic_dir, project, "_pyscenic_gene_regulons_heatmap.pdf", sep=""), width=15, height=30)
par(bg=NA)
heatmap_gene
dev.off()

