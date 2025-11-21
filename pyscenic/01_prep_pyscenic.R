##### Prep pySCENIC object 
##### Load seurat object 
##### Extract RNA assay 
##### Filter for only those genes in pyscenic databases 
##### Convert to loom object and save 

# load packages
library(arrow)
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

# set rds directory path
guide_dir <- paste0(seurat_dir, "guide_detection/")

# set pseudo results directory path
pseudo_dir <- paste0(seurat_dir, "pseuodtime_analyses/")

# set pyscenic results directory path - create
pyscenic_dir <- paste0(seurat_dir, "pyscenic_analyses/")


##### STEP ONE - Read in Seurat object #####

print("Reading in integrated Seurat object")

# read in multiome object
MultiOme <- readRDS(paste0(guide_dir, project, "_RNA_ATAC_MGATK_GUIDES_combined.rds"))
MultiOme

head(MultiOme[[]])
colnames(MultiOme[[]])

# set default assay
DefaultAssay(MultiOme) <- "RNA"

# check layers and join if necessary
Layers(MultiOme[["RNA"]])

# join layers
MultiOme[["RNA"]] <- JoinLayers(MultiOme[["RNA"]])

# check layers 
Layers(MultiOme[["RNA"]])


##### STEP TWO - Extract data #####

# extract metadata
meta <- MultiOme@meta.data
write.csv(meta, paste(pyscenic_dir, project, "_metadata.csv", sep = ""), row.names = TRUE)

# set default assay
DefaultAssay(MultiOme) <- "RNA"

# extract expression matrix (cells × genes)
expr_matrix <- as.matrix(GetAssayData(MultiOme, assay = "RNA", slot = "data"))
dim(expr_matrix)

# read in the four feather files
feather_dir <- "/../../pyscenic_databases/"
feather_files <- c(
  "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
  "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather",
  "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
  "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"
)

all_genes <- unique(unlist(lapply(file.path(feather_dir, feather_files), function(f) {
  message("Reading genes from: ", f)
  cols <- colnames(read_feather(f))
  cols[cols != "motif"]  # remove the motif column
})))

length(all_genes)  

genes_to_keep <- intersect(rownames(expr_matrix), all_genes)
expr_matrix_filtered <- expr_matrix[genes_to_keep, , drop = FALSE]
dim(expr_matrix_filtered)

# extract metadata
cell_metadata <- MultiOme@meta.data
cell_ids <- rownames(cell_metadata)
metadata_to_keep <- c("orig.ident", "guide", "gene") 

# convert metadata columns to character and replace NAs
metadata_list <- lapply(metadata_to_keep, function(col) {
  x <- cell_metadata[[col]]
  x <- as.character(x)                  
  x[is.na(x)] <- ""                     
  return(x)
})

names(metadata_list) <- metadata_to_keep

# combine with cell_ids
col_attrs <- c(list(CellID = cell_ids), metadata_list)

# row metadata
gene_names <- rownames(expr_matrix_filtered)
row_attrs <- list(Gene = gene_names)

# loom output path
loom_file <- paste0(pyscenic_dir, project, "_rna_assay_filtered.loom")

# loomR::create will never overwrite an existing file
# remove existing file
if (file.exists(loom_file)) {
  file.remove(loom_file)
}

# create loom file (genes × cells)
lp <- loomR::create(
  filename = loom_file,
  data = expr_matrix_filtered,  
  row.attrs = row_attrs,
  col.attrs = col_attrs
)

lp$close_all()

