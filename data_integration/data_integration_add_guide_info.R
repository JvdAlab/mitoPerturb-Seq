##### Addition of guide/gene information

# load packages
library(clusterProfiler)
library(clustree)
library(ComplexHeatmap)
library(cowplot)
library(DESeq2)
library(dplyr)
library(fgsea)
library(future)
library(genomation)
library(GenomicRanges)
library(ggalluvial)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(glmGamPoi)
library(gprofiler2)
library(irlba)
library(matrixStats)
library(mixtools)
library(patchwork)
library(pheatmap)
library(presto)
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
library(viridis)

library(AnnotationDbi)
library(AnnotationHub)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(msigdbr)
library(biovizBase)


#### Set Project Variables ####

# enter project name
project <- "project"

# set base directory
baseDir <- "/../../"
setwd(baseDir)


##### STEP ONE - Load Sample1 data #####

# set guide results folder
guide_folder_sampleone <- "/../../Sample1_guide_detection/"

# load the Seurat object 
MultiOme_obj_sampleone <- readRDS(paste0(guide_folder_sampleone,"Sample1_GuideEnrich_Detection.rds"))

# extract cell metadata 
cells_sampleone <- MultiOme_obj_sampleone@meta.data

# create new column called cell_barcode 
cells_sampleone$cell <- rownames(cells_sampleone)

# count number of rows
nrow(cells_sampleone) 

# add project ID to cell names 
cells_sampleone$cell <-paste0("Sample1_", cells_sampleone$cell)

# remove object
rm(MultiOme_obj_sampleone)


##### STEP TWO - Load Sample2 data #####

# set guide results folder
guide_folder_sampletwo <- "/../../Sample2_guide_detection/"

# load the Seurat object 
MultiOme_obj_sampletwo <- readRDS(paste0(guide_folder_sampletwo,"Sample2_GuideEnrich_Detection.rds"))

# extract cell metadata 
cells_sampletwo <- MultiOme_obj_sampletwo@meta.data

# create new column called cell_barcode 
cells_sampletwo$cell <- rownames(cells_sampletwo)

# count number of rows
nrow(cells_sampletwo)  

# add project ID to cell names - this should then match metadata from integrated analysis
cells_sampletwo$cell <-paste0("Sample2_", cells_sampletwo$cell)

# remove object
rm(MultiOme_obj_sampletwo)


##### STEP THREE - Join guide information together #####

# join together sampleone and sampletwo metadata 
p1_p2_guide_info <- rbind(cells_sampleone, cells_sampletwo)
head(p1_p2_guide_info)
tail(p1_p2_guide_info)
nrow(p1_p2_guide_info)


##### STEP FOUR - Load integrated data #####

# set integrated rds results folder
integrate <- "/../../rds_files/"

# load the Seurat object following integration
MultiOme_obj_int <- readRDS(paste0(integrate,"project_RNA_ATAC_MGATK_combined.rds"))
MultiOme_obj_int


##### STEP FIVE - Add guide information to integrated data #####

# add guide and gene columns to integrated metadata
MultiOme_obj_int@meta.data$guide <- p1_p2_guide_info[match(rownames(MultiOme_obj_int@meta.data) , p1_p2_guide_info$cell),]$top_guide_final_enrich
MultiOme_obj_int@meta.data$gene <- p1_p2_guide_info[match(rownames(MultiOme_obj_int@meta.data) , p1_p2_guide_info$cell),]$top_guide_final_enrich_combined

# examine data 
table(MultiOme_obj_int@meta.data$guide) 

# label those with no guide information as unknown 
MultiOme_obj_int@meta.data$guide[is.na(MultiOme_obj_int@meta.data$guide)] <- 'unknown'

# examine data 
table(MultiOme_obj_int@meta.data$gene) 

# label those with no gene information as unknown 
MultiOme_obj_int@meta.data$gene[is.na(MultiOme_obj_int@meta.data$gene)] <- 'unknown'


##### STEP SIX - Save data #####

# set results folder 
guide_folder <- "/../../guide_detection/"

# save data
write.csv(MultiOme_obj_int@meta.data[,c("nCount_mito", "Heteroplasmy_all_5024", "Depth_all_5024", "guide", "gene")],  
          paste(guide_folder, project, "_mtDNA_heteroplasmy_with_guides_results.csv", sep = ""))

# save file as RDS
saveRDS(MultiOme_obj_int, file="/../../guide_detection/project_RNA_ATAC_MGATK_GUIDES_combined.rds")

