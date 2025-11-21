##### Mixscape analysis

# load packages
library(ComplexHeatmap)
library(cowplot)
library(dplyr)
library(future)
library(genomation)
library(GenomicRanges)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(glmGamPoi)
library(gprofiler2)
library(irlba)
library(matrixStats)
library(mixtools) 
library(patchwork)
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


##### STEP ONE - Load integrated data #####

# set integrated rds results folder
integrate <- "/../../guide_detection/"

# load the Seurat object following integration and addition of guide annotation
MultiOme_obj_no_unkmul <- readRDS(paste0(integrate,"project_RNA_ATAC_MGATK_GUIDES_combined.rds"))
MultiOme_obj_no_unkmul

# remove unknown guides, multiple guides and polg-6
# this needs to be done in three steps
MultiOme_obj_no_unkmul <- subset(MultiOme_obj_no_unkmul, top_guide_final_enrich != "Polg-6")
MultiOme_obj_no_unkmul <- subset(MultiOme_obj_no_unkmul, top_guide_final_enrich_combined != "unknown")
MultiOme_obj_no_unkmul <- subset(MultiOme_obj_no_unkmul, top_guide_final_enrich_combined != "multiple_guides")

# create guide column
MultiOme_obj_no_unkmul@meta.data$guide <- MultiOme_obj_no_unkmul@meta.data$top_guide_final_enrich
MultiOme_obj_no_unkmul@meta.data$guide <- ifelse(MultiOme_obj_no_unkmul@meta.data$top_guide_final_enrich %in% c("Eomes-1", "Eomes-2", "Eomes-3", "Neurod1-2", "Neurod1-3", "Neurod1-4", 
                                                                                                                "Nontargeting-1", "Nontargeting-2", "Nontargeting-3", "Olig1-2", "Olig1-3", 
                                                                                                                "Olig1-4", "Opn4-2", "Opn4-3", "Opn4-6", "Rgr-1", "Rgr-2", "Rgr-3", "Rrh-1", 
                                                                                                                "Rrh-3", "Rrh-6"), "NT", MultiOme_obj_no_unkmul@meta.data$top_guide_final_enrich)

# create gene column
MultiOme_obj_no_unkmul@meta.data$gene <- gsub( "-.*", "", MultiOme_obj_no_unkmul@meta.data$guide)

# create crispr column
MultiOme_obj_no_unkmul@meta.data$crispr <- ifelse(MultiOme_obj_no_unkmul@meta.data$top_guide_final_enrich %in% c("Eomes-1", "Eomes-2", "Eomes-3", "Neurod1-2", "Neurod1-3", "Neurod1-4", 
                                                                                                                 "Nontargeting-1", "Nontargeting-2", "Nontargeting-3", "Olig1-2", "Olig1-3", 
                                                                                                                 "Olig1-4", "Opn4-2", "Opn4-3", "Opn4-6", "Rgr-1", "Rgr-2", "Rgr-3", "Rrh-1", 
                                                                                                                 "Rrh-3", "Rrh-6"), "NT", "Perturbed")  

# set results directory 
mixscape_dir <- "/../../mixscape/"


##### STEP TWO - RNA clustering #####

# set default assay
DefaultAssay(MultiOme_obj_no_unkmul) <- 'RNA'

# need to join layers in order to use RNA assay 
MultiOme_obj_no_unkmul[["RNA"]] <- JoinLayers(MultiOme_obj_no_unkmul[["RNA"]])

# prepare RNA assay for dimensionality reduction: Normalize data, find variable features and scale data
MultiOme_obj_no_unkmul <- NormalizeData(MultiOme_obj_no_unkmul) %>% FindVariableFeatures() %>% ScaleData()

# run PCA 
MultiOme_obj_no_unkmul <- RunPCA(MultiOme_obj_no_unkmul, reduction.name = "mixscape_pca")

# run UMAP 
MultiOme_obj_no_unkmul <- RunUMAP(MultiOme_obj_no_unkmul, dims = 1:40, reduction = "mixscape_pca", reduction.name = "umap_mixscape", reduction.key = "UMAPmixscape_")


##### STEP THREE - Calculate local perturbation signatures #####

# Calculate perturbation signature (PRTB)
MultiOme_obj_no_unkmul <- CalcPerturbSig(
  object = MultiOme_obj_no_unkmul, # An object of class Seurat
  assay = "RNA", # Name of Assay PRTB signature is being calculated on
  slot = "data", # Data slot to use for PRTB signature calculation
  gd.class ="gene", # Metadata column containing target gene classification
  nt.cell.class = "NT", # Non-targeting gRNA cell classification identity
  reduction = "mixscape_pca", # Reduction method used to calculate nearest neighbors
  ndims = 40, # Number of dimensions to use from dimensionality reduction method
  num.neighbors = 20, # Number of nearest neighbors to consider
  split.by = "orig.ident", # Provide metadata column if multiple biological replicates exist to calculate PRTB signature for every replicate separately
  new.assay.name = "PRTB") # Name for the new assay

# set default assay
DefaultAssay(MultiOme_obj_no_unkmul) <- "PRTB"

# use variable features from RNA assay
VariableFeatures(MultiOme_obj_no_unkmul) <- VariableFeatures(MultiOme_obj_no_unkmul[["RNA"]])
MultiOme_obj_no_unkmul <- ScaleData(MultiOme_obj_no_unkmul, do.scale = F, do.center = T)

# run PCA 
MultiOme_obj_no_unkmul <- RunPCA(MultiOme_obj_no_unkmul, reduction.name = "prtb_pca", reduction.key = "prtbpca_")

# Run UMAP 
MultiOme_obj_no_unkmul <- RunUMAP(object = MultiOme_obj_no_unkmul, dims = 1:40, reduction = "prtb_pca", reduction.name = "umap_prtb", reduction.key = "UMAPprtb_")


##### STEP FOUR - Run mixscape #####

# set default assay
DefaultAssay(MultiOme_obj_no_unkmul) <- "PRTB"

# run mixscape
MultiOme_obj_no_unkmul <- RunMixscape(
  object = MultiOme_obj_no_unkmul, # An object of class Seurat
  assay = "PRTB", # Assay to use for mixscape classification
  slot = "scale.data", # Assay data slot to use
  labels = "gene", # Metadata column with target gene labels
  nt.class.name = "NT", # Classification name of non-targeting gRNA cells
  min.de.genes = 5, # Required number of genes that are differentially expressed for method to separate perturbed and non-perturbed cells
  iter.num = 10, # Number of normalmixEM iterations to run if convergence does not occur
  de.assay = "RNA", # Assay to use when performing differential expression analysis. Usually RNA
  split.by = "orig.ident", # Metadata column with experimental condition/cell type classification information
  verbose = T, # Display messages
  prtb.type = "KO") # specify type of CRISPR perturbation expected for labeling mixscape classifications. Default is KO


##### STEP FIVE - Plot results #####

# Explore the perturbation scores of cells - only those with KO can be plotted 
for(i in c("Atg5", "Opa1", "Polg", "Tfam")){ 
  print(i)
  pdf(paste0(mixscape_dir, project, "_PlotPerturbScore_" , i,  ".pdf"), onefile=FALSE, width=8, height=6) 
  par(bg=NA)
  print({PlotPerturbScore(object = MultiOme_obj_no_unkmul,  
                          target.gene.ident = i, 
                          mixscape.class = "mixscape_class",  
                          col = "coral2") + labs(fill = "mixscape class")})
  dev.off()
}

# perform DE
Idents(MultiOme_obj_no_unkmul) <- "gene"

# balanced = T 
MixscapeHeatmap(object = MultiOme_obj_no_unkmul, # An object of class Seurat
                ident.1 = "NT", # Identity class to define markers for
                ident.2 = "Atg5", # A second identity class for comparison
                balanced = T, # Plot an equal number of genes with both groups of cells
                assay = "RNA", # Assay to use in differential expression testing
                max.genes = 30, angle = 0, # Total number of DE genes to plot
                group.by = "mixscape_class", # Metadata column with mixscape classifications, originally: group.by (split densities based on mixscape classification)
                max.cells.group = 300, # Number of cells per identity to plot
                size=3) + NoLegend() + theme(axis.text.y = element_text(size = 10))

# extract metadata and save to csv 
cells_meta_all <- MultiOme_obj_no_unkmul@meta.data
write.csv(cells_meta_all, paste(mixscape_dir, project, "_mixscape_metadata.csv", sep = ""), row.names = TRUE)


##### STEP SIX - Visualise perturbation responses with LDA #####

# set idents
Idents(MultiOme_obj_no_unkmul) <- "mixscape_class.global" # KO, NP, NT

# subset data to only those with KO or NT classification
KO_NT_sub <- subset(MultiOme_obj_no_unkmul, idents = c("KO", "NT")) 
KO_NT_sub

# Run LDA 
KO_NT_sub <- MixscapeLDA(
  object = KO_NT_sub, # An object of class Seurat
  assay = "RNA", # Assay to use for performing Linear Discriminant Analysis (LDA)
  pc.assay = "PRTB", # Assay to use for running Principle components analysis
  labels = "gene", # Meta data column with target gene class labels
  nt.label = "NT", # Name of non-targeting cell class
  npcs = 10, # Number of principle components to use
  logfc.threshold = 0.25, # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells (default 0.25)
  verbose = T) # Print progress bar

# run UMAP 
KO_NT_sub <- RunUMAP(object = KO_NT_sub, dims = 1:4, reduction = "lda", reduction.name = "umap_lda", reduction.key = "UMAPlda_")

# extract metadata and save to csv so can be loaded in and compared to p1/p2 mixscape results
cells_meta_sub <- KO_NT_sub@meta.data
write.csv(cells_meta_sub, paste(mixscape_dir, project, "_mixscape_metadata_KO_NT_sub.csv", sep = ""), row.names = TRUE)


##### STEP SEVEN - Save data as RDS #####

# save file as RDS
saveRDS(MultiOme_obj_no_unkmul, file="/../../mixscape_dir/project_no_unkmul_polg6_RNA_ATAC_MGATK_GUIDES_mixscape.rds")

