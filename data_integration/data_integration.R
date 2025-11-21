##### Data integration 

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

# enter sample names 
sample_list <- c("Sample1", "Sample2")

# set algorithm to use for FindClusters 
# 1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm
algorithm <- 3

# set resolution to use for FindClusters
resolution <- 2

# set no dimensions to use when calculating RNA reductions
RNA_dims <- 1:20

# set no dimensions to use when calculating RNA umaps
RNA_umap_dims <- 1:50

# set no dimensions to use when calculating ATAC reductions
ATAC_dims <- 2:40

# set no dimensions to use when integrating ATAC data
ATAC_integrate_dims <- 1:40

# set minimum cutoff for FindTopFeatures
min_cutoff <- 10

# set path to MACS2
macs2_path <- "/../../miniconda/envs/PeakCalling_analysis/bin/macs2"

# set directory containing cellranger and mgatk output directories
baseDir <- "/../../"
setwd(baseDir)

# create seurat directory
dir.create(paste0(baseDir, "seurat/"))

# set seurat directory path
seurat_dir <- paste0(baseDir, "seurat/")

# create rds_files directory
dir.create(paste0(seurat_dir, "rds_files/"))

# set rds_files directory path
rds_dir <- paste0(seurat_dir, "rds_files/")

# create peak_calling directory
dir.create(paste0(seurat_dir, "peak_calling/"))

# set peak_calling directory path
peaks_dir <- paste0(seurat_dir, "peak_calling/")

# create integrated directory
dir.create(paste0(seurat_dir, "integrated/"))

# set integrated directory path
integrated_dir <- paste0(seurat_dir, "integrated/")

# create mgatk directory
dir.create(paste0(integrated_dir, "mgatk/"))

# set mgatk directory path
mgatk_dir <- paste0(integrated_dir, "mgatk/")


#### Generate combined_peaks ####

# Generate combined_peaks bed file from cell ranger output to use as common reference for creating ATAC assays, helps with merging/integrating (https://stuartlab.org/signac/articles/merging.html)
gr_list_cr <- c()

for (sample in sample_list) {

  peaks <- read.table(file = paste0(baseDir, "cell_ranger/", sample, "/outs/atac_peaks.bed"),
                      col.names = c("chr", "start", "end"))
  peaks <- makeGRangesFromDataFrame(peaks)
  assign(paste0("gr_", sample), peaks)
  gr_list_cr <- append(gr_list_cr, get(paste0("gr_", sample)))

}

combined_peaks_cr <- reduce(x = gr_list_cr)

peakwidths_cr <- width(combined_peaks_cr)
combined_peaks_cr <- combined_peaks_cr[peakwidths_cr < 10000 & peakwidths_cr >20]

export(combined_peaks_cr, paste0(peaks_dir, "combined_cr_peaks.bed"))

# Generate combined_peaks bed file using MACS2 to use as common reference for creating ATAC assays, helps with merging/integrating (https://stuartlab.org/signac/articles/merging.html)
gr_list_macs <- c()

for (sample in sample_list) {

  fragpath <- paste0(baseDir, "cell_ranger/", sample, "/outs/atac_fragments.tsv.gz")
  frags <- CreateFragmentObject(fragpath)
  assign(paste0("frags_", sample), frags)
  peaks <- CallPeaks(fragpath, macs2.path = macs2_path)
  assign(paste0("gr_", sample), peaks)
  gr_list_macs <- append(gr_list_macs, get(paste0("gr_", sample)))

}

combined_peaks_macs <- reduce(x = gr_list_macs)

peakwidths_macs <- width(combined_peaks_macs)
combined_peaks_macs <- combined_peaks_macs[peakwidths_macs < 10000 & peakwidths_macs >20]

export(combined_peaks_macs, paste0(peaks_dir, "combined_macs_peaks.bed"))

for (sample in sample_list) {
  # set sample id
  sample_id <- sample
  print(sample_id)

  ##### STEP ONE - Create ATAC Multiome Objects ####

  # read in 10x data
  MultiOme_data <- Read10X_h5(paste0(baseDir,"cell_ranger/", sample_id, "/outs/filtered_feature_bc_matrix.h5"))

  # get ensembl gene annotations for mouse
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotation) <- "UCSC"

  # read in ATAC fragments file
  fragpath <- paste0(baseDir,"cell_ranger/", sample_id, "/outs/atac_fragments.tsv.gz")

  # load in cellranger metadata (only use for combined_peaks method)
  md <- read.table(paste0(baseDir,"cell_ranger/", sample_id, "/outs/per_barcode_metrics.csv"), stringsAsFactors = FALSE, sep = ",", header = TRUE)

  # filter on true cells only (only use for combined_peaks method)
  md <- md[md$is_cell == 1, ]

  combined_macs_peaks <- readGeneric(paste0(peaks_dir, "combined_macs_peaks.bed"), chr = 1, start = 2, end = 3)

  # extract fragments for cells and calculate counts (only use for combined_peaks method)
  frags <- CreateFragmentObject(fragpath)

  counts <- FeatureMatrix(fragments = frags, features = combined_macs_peaks, cells = md$barcode)

  # generate Seurat object (only use for combined_peaks method)
  chromatin_assay <- CreateChromatinAssay(counts, fragments = fragpath, annotation = annotation)
  MultiOme_ATAC <- CreateSeuratObject(chromatin_assay, assay = "ATAC", project = sample_id)

  atac_cells <- Cells(MultiOme_ATAC)
  MultiOme_ATAC[["cell_id"]] <- atac_cells

  ##### STEP TWO - Create RNA Multiome object #####

  # create Seurat object for RNA
  MultiOme_RNA <- CreateSeuratObject(counts = MultiOme_data$`Gene Expression`, assay = "RNA", project = sample_id)

  rna_cells <- Cells(MultiOme_RNA)
  MultiOme_RNA[["cell_id"]] <- rna_cells

  ##### STEP THREE - Quality control #####

  # create QC results folder
  dir.create(paste0(seurat_dir, sample_id, "/seurat_qc/"), recursive = TRUE)

  # set QC results folder
  qc_folder <- paste0(seurat_dir, sample_id, "/seurat_qc/")

  ### RNA Quality control ###

  # add mitochondrial DNA column to meta data
  grep("^mt-",rownames(MultiOme_RNA),value = TRUE)
  MultiOme_RNA[["percent.mt"]] <- PercentageFeatureSet(MultiOme_RNA, pattern = "^mt-")
  head(MultiOme_RNA[["percent.mt"]])

  # add ribosomal column to meta data
  grep("^Rp[ls]",rownames(MultiOme_RNA),value = TRUE)
  MultiOme_RNA[["percent.Ribosomal"]] <- PercentageFeatureSet(MultiOme_RNA, pattern = "^Rp[ls]")
  head(MultiOme_RNA[["percent.Ribosomal"]])

  # Normalise and scale data for cell-cycle scoring
  MultiOme_RNA <- NormalizeData(MultiOme_RNA)

  MultiOme_RNA <- FindVariableFeatures(MultiOme_RNA, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(MultiOme_RNA), 10)

  all.genes <- rownames(MultiOme_RNA)
  MultiOme_RNA <- ScaleData(MultiOme_RNA, features = all.genes)

  # run PCA
  # creates "pca" dimensional reduction
  MultiOme_RNA <- RunPCA(MultiOme_RNA, features = VariableFeatures(object = MultiOme_RNA))

  # jackstraw plot
  MultiOme_RNA <- JackStraw(MultiOme_RNA, num.replicate = 100)
  MultiOme_RNA <- ScoreJackStraw(MultiOme_RNA, dims = 1:20)

  # find nearest neighbours - change dims based on elbow plot
  MultiOme_RNA <- FindNeighbors(MultiOme_RNA, dims = 1:20)

  # find clusters
  MultiOme_RNA <- FindClusters(MultiOme_RNA, resolution = c(0.5,1,1.5))

  # run UMAP - change dims based on elbow plot
  # creates "umap" dimensional reduction
  MultiOme_RNA <- RunUMAP(MultiOme_RNA, dims = 1:20)

  mmus_s.genes = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens",
                       target_organism = "mmusculus")$ortholog_name

  mmus_g2m.genes = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens",
                         target_organism = "mmusculus")$ortholog_name

  # set default assay as RNA 
  DefaultAssay(MultiOme_RNA) <- "RNA"

  # run cell cycle scoring
  MultiOme_RNA <- CellCycleScoring(MultiOme_RNA, s.features = mmus_s.genes, g2m.features = mmus_g2m.genes, set.ident = TRUE)

  # run a pca on cell cycle genes to reveal if cells separate by phase
  # creates "pca" dimensional reduction
  MultiOme_RNA <- RunPCA(MultiOme_RNA, approx = FALSE, features=c(mmus_s.genes, mmus_g2m.genes))

  ### ATAC quality control for cell ranger peaks ###

  # set default assay as ATAC
  DefaultAssay(MultiOme_ATAC) <- "ATAC"

  # quality control steps - nucleosome signal and transcription start site enrichment score
  MultiOme_ATAC <- NucleosomeSignal(MultiOme_ATAC)
  MultiOme_ATAC <- TSSEnrichment(MultiOme_ATAC)

  ### RNA and ATAC Quality control plots ###

  # plot quality control metrics
  pdf(paste0(qc_folder, sample_id, "_RNA_QC_metrics_violin.pdf"), onefile=FALSE, width=10, height=10)
  par(bg=NA)
  x <- VlnPlot(object = MultiOme_RNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.Ribosomal"), ncol = 4, pt.size = 0.0001)
  print(x)
  dev.off()

  pdf(paste0(qc_folder, sample_id, "_ATAC_QC_metrics_violin.pdf"), onefile=FALSE, width=10, height=10)
  par(bg=NA)
  x <- VlnPlot(object = MultiOme_ATAC, features = c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "TSS.enrichment"), ncol = 4, pt.size = 0.0001)
  print(x)
  dev.off()

  ##### STEP FOUR - Filter out low quality cells #####

  # Filter low quality cells - count RNA, mtDNA
  MultiOme_RNA <- subset(x = MultiOme_RNA,
                         subset = nCount_RNA > 1000 & nCount_RNA < 60000 &
                           nFeature_RNA > 1000 & nFeature_RNA < 10000 &
                           percent.mt < 10 & percent.Ribosomal < 30
  )

  # Filter low quality cells for cr peaks - count ATAC, nucleosome signal and TSS enrichment
  MultiOme_ATAC <- subset(x = MultiOme_ATAC,
                          subset = nCount_ATAC > 1000 & nCount_ATAC < 150000 &
                          nFeature_ATAC > 500 & nFeature_ATAC < 55000 &
                          nucleosome_signal < 2 & TSS.enrichment > 1
  )

  # Extract cell_id post-filtering
  rna_cells_filtered <- Cells(MultiOme_RNA)
  atac_cells_filtered <- Cells(MultiOme_ATAC)

  # Filter RNA based on ATAC cells
  MultiOme_RNA <- subset(x = MultiOme_RNA,
                         subset = cell_id %in% atac_cells_filtered)

  # Filter ATAC based on RNA cells
  MultiOme_ATAC <- subset(x = MultiOme_ATAC,
                          subset = cell_id %in% rna_cells_filtered)

  #Extract RNA & ATAC filtered cell_id
  rna_atac_cells_filtered <- Cells(MultiOme_RNA)

  # Process ATAC_cr object
  MultiOme_ATAC <- FindTopFeatures(MultiOme_ATAC, min.cutoff = min_cutoff)
  MultiOme_ATAC <- RunTFIDF(MultiOme_ATAC)
  MultiOme_ATAC <- RunSVD(MultiOme_ATAC)

  ##### STEP FIVE - Create mgatk Seurat Object #####

  # create mgatk results folder
  dir.create(paste0(seurat_dir, sample_id, "/seurat_mgatk/"))

  # set mgatk results folder
  mgatk_folder <- paste0(seurat_dir, sample_id, "/seurat_mgatk/")

  # load mgatk output
  mito.data <- ReadMGATK(paste0(baseDir, "mgatk/", sample_id, "/final"))

  # create Seurat object for mgatk (mtdna variant) data
  MultiOme_MGATK <- CreateSeuratObject(counts = mito.data$counts, assay = "mito", project = sample_id)

  # Subset to only those cells present in the scATAC-seq assay
  MultiOme_MGATK <- subset(MultiOme_MGATK, cells = rna_atac_cells_filtered)

  # add metadata to the seurat object
  MultiOme_MGATK <- AddMetaData(MultiOme_MGATK, metadata = mito.data$depth, col.name = "mtDNA_depth")

  # save objects as RDS
  saveRDS(MultiOme_RNA, file=paste0(rds_dir, sample_id, "_RNA_processed.rds"))
  saveRDS(MultiOme_ATAC, file=paste0(rds_dir, sample_id, "_ATAC_processed.rds"))
  saveRDS(MultiOme_MGATK, file=paste0(rds_dir, sample_id, "_MGATK_filtered.rds"))

  rm(MultiOme_RNA, MultiOme_ATAC, MultiOme_MGATK, mito.data, MultiOme_data, frags, counts, chromatin_assay)
  gc()
}


##### Merge & Integrate Data   #####

#### Process RNA files ####

# Based on https://satijalab.org/seurat/articles/integration_introduction.html

# load in filtered RNA rds files and create list of RNA objects
RNA_rds_list <- c()

for (sample in sample_list) {
  assign(paste0("MultiOme_RNA_", sample), readRDS(paste0(rds_dir, sample, "_RNA_processed.rds")))
  RNA_rds_list <- append(RNA_rds_list, get(paste0("MultiOme_RNA_", sample)))
}

# remove first object from list, as this is the object that the others will be merged to
y <- RNA_rds_list[-1]

# merge into a single Seurat object 
MultiOme_RNA_merged <- merge(RNA_rds_list[[1]],
                             y = y,
                             add.cell.ids = sample_list,
                             project = paste0(project, "_RNA"))

# Normalize & cluster with no integration
MultiOme_RNA_merged <- SCTransform(MultiOme_RNA_merged, vars.to.regress = c("S.Score", "G2M.Score"), return.only.var.genes = FALSE)
MultiOme_RNA_merged <- RunPCA(MultiOme_RNA_merged, reduction.name = "merged_pca")
MultiOme_RNA_merged <- FindNeighbors(MultiOme_RNA_merged, dims = RNA_dims, reduction = "merged_pca")
MultiOme_RNA_merged <- FindClusters(MultiOme_RNA_merged, algorithm = algorithm, resolution = resolution, cluster.name = "rna_unintegrated_clusters")

MultiOme_RNA_merged <- RunUMAP(MultiOme_RNA_merged, dims = RNA_umap_dims, reduction = "merged_pca", reduction.name = "rna_unintegrated", reduction.key = "UMAPrnaunint_")
u1 <- DimPlot(MultiOme_RNA_merged, reduction = "rna_unintegrated", group.by = c("orig.ident", "rna_unintegrated_clusters"))
ggsave(paste0(integrated_dir, "all_RNA_unintegrated_UMAP.pdf"), plot = u1, width = 16, height = 8)

saveRDS(MultiOme_RNA_merged, file=paste0(rds_dir, project, "_RNA_unintegrated.rds"))

# Integrate layers
MultiOme_RNA_integrated <- IntegrateLayers(object = MultiOme_RNA_merged, method = CCAIntegration, normalization.method = "SCT", orig.reduction = "merged_pca",
                                       new.reduction = "integrated.cca", scale.layer = "scale.intdata")

# Cluster after integration
MultiOme_RNA_integrated <- FindNeighbors(MultiOme_RNA_integrated, dims = RNA_dims, reduction = "integrated.cca")
MultiOme_RNA_integrated <- FindClusters(MultiOme_RNA_integrated, algorithm = algorithm, resolution = resolution, cluster.name = "rna_integrated_clusters")
MultiOme_RNA_integrated <- RunUMAP(MultiOme_RNA_integrated, dims = RNA_umap_dims, reduction = 'integrated.cca', reduction.name = "rna_integrated", reduction.key = "UMAPrnaint_")
u2 <- DimPlot(MultiOme_RNA_integrated, reduction = "rna_integrated", group.by = c("orig.ident", "rna_integrated_clusters"))
ggsave(paste0(integrated_dir, "all_RNA_integrated_UMAP.pdf"), plot = u2, width = 16, height = 8)

# save integrated Seurat object
saveRDS(MultiOme_RNA_integrated, file=paste0(rds_dir, project, "_RNA_integrated.rds"))

#### Process ATAC files ####

# Based on https://stuartlab.org/signac/articles/integrate_atac

# load in processed ATAC rds files and create list of ATAC objects
ATAC_rds_list <- c()

for (sample in sample_list) {
  assign(paste0("MultiOme_ATAC_", sample), readRDS(paste0(rds_dir, sample, "_ATAC_processed.rds")))
  ATAC_rds_list <- append(ATAC_rds_list, get(paste0("MultiOme_ATAC_", sample)))
}

# remove first object from list, as this is the object that the others will be merged to
y <- ATAC_rds_list[-1]

# merge into a single Seurat object
MultiOme_ATAC_merged <- merge(ATAC_rds_list[[1]],
                             y = y,
                             add.cell.ids = sample_list,
                             project = paste0(project, "_ATAC"),
                             merge.data = TRUE)

# set default assay as ATAC
DefaultAssay(MultiOme_ATAC_merged) <- "ATAC"

# process merged Seurat object without integration
MultiOme_ATAC_merged <- FindTopFeatures(MultiOme_ATAC_merged, min.cutoff = min_cutoff)
MultiOme_ATAC_merged <- RunTFIDF(MultiOme_ATAC_merged)
MultiOme_ATAC_merged <- RunSVD(MultiOme_ATAC_merged)

# drop first dimension, as this is usually correlated with sequencing depth
MultiOme_ATAC_merged <- FindNeighbors(MultiOme_ATAC_merged, reduction = 'lsi', dims = ATAC_dims)
MultiOme_ATAC_merged <- FindClusters(MultiOme_ATAC_merged, algorithm = algorithm, resolution = resolution, cluster.name = "atac_unintegrated_clusters")
MultiOme_ATAC_merged <- RunUMAP(MultiOme_ATAC_merged, dims = ATAC_dims, reduction = "lsi", reduction.name = "atac_unintegrated", reduction.key = "UMAPatacunint_")
u3 <- DimPlot(MultiOme_ATAC_merged, reduction = "atac_unintegrated", group.by = c("orig.ident", "atac_unintegrated_clusters"))
ggsave(paste0(integrated_dir, "all_ATAC_unintegrated_UMAP.pdf"), plot = u3, width = 16, height = 8)

saveRDS(MultiOme_ATAC_merged, file=paste0(rds_dir, project, "_ATAC_unintegrated.rds"))

# create empty list for ATAC objects
atac_object_list <- list()

# add cell ids to original objects to prevent cell names clashing & add objects to object_list
for (sample in sample_list) {
  assign(paste0("MultiOme_ATAC_", sample), RenameCells(get(paste0("MultiOme_ATAC_", sample)), add.cell.id = sample))
  atac_object_list <- append(atac_object_list, get(paste0("MultiOme_ATAC_", sample)))
}

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = atac_object_list,
  reduction = "rlsi",
  dims = ATAC_dims
)

MultiOme_ATAC_integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = MultiOme_ATAC_merged[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = ATAC_integrate_dims
)

# set default assay as ATAC
DefaultAssay(MultiOme_ATAC_integrated) <- "ATAC"

# process merged Seurat object post integration
MultiOme_ATAC_integrated <- FindTopFeatures(MultiOme_ATAC_integrated, min.cutoff = min_cutoff)
MultiOme_ATAC_integrated <- RunTFIDF(MultiOme_ATAC_integrated)
MultiOme_ATAC_integrated <- RunSVD(MultiOme_ATAC_integrated)

# drop first dimension, as this is usually correlated with sequencing depth
MultiOme_ATAC_integrated <- FindNeighbors(MultiOme_ATAC_integrated, reduction = 'integrated_lsi', dims = ATAC_dims)
MultiOme_ATAC_integrated <- FindClusters(MultiOme_ATAC_integrated, algorithm = algorithm, resolution = resolution, cluster.name = "atac_integrated_clusters")
MultiOme_ATAC_integrated <- RunUMAP(MultiOme_ATAC_integrated, dims = ATAC_dims, reduction = "integrated_lsi", reduction.name = "atac_integrated", reduction.key = "UMAPatacint_")
u4 <- DimPlot(MultiOme_ATAC_integrated, reduction = "atac_integrated", group.by = c("orig.ident", "atac_integrated_clusters"))
ggsave(paste0(integrated_dir, "all_ATAC_integrated_UMAP.pdf"), plot = u4, width = 16, height = 8)

saveRDS(MultiOme_ATAC_integrated, file=paste0(rds_dir, project, "_ATAC_integrated.rds"))

#### Combine RNA & ATAC ####

MultiOme_RNA_integrated <- readRDS(paste0(rds_dir, project, "_RNA_integrated.rds"))
MultiOme_ATAC_integrated <- readRDS(paste0(rds_dir, project, "_ATAC_integrated.rds"))

MultiOme_RNA_integrated[["ATAC"]] <- MultiOme_ATAC_integrated[["ATAC"]]
MultiOme_RNA_integrated[["integrated_lsi"]] <- MultiOme_ATAC_integrated[["integrated_lsi"]]
MultiOme_RNA_integrated[["atac_integrated"]] <- MultiOme_ATAC_integrated[["atac_integrated"]] 

# extract and import ATAC metadata
ATAC_meta <- MultiOme_ATAC_integrated@meta.data
MultiOme_RNA_integrated <- AddMetaData(MultiOme_RNA_integrated, metadata = ATAC_meta)

# build a joint neighbor graph using both assays
MultiOme_RNA_integrated <- FindMultiModalNeighbors(MultiOme_RNA_integrated,
  reduction.list = list("integrated.cca", "integrated_lsi"),
  dims.list = list(RNA_umap_dims, ATAC_dims),
  modality.weight.name = c("RNA.weight", "ATAC.weight"),
  verbose = TRUE
)

# build a joint UMAP visualisation
# creates "umap_wnn" dimensional reduction
MultiOme_RNA_integrated <- RunUMAP(
  object = MultiOme_RNA_integrated,
  nn.name = "weighted.nn",
  reduction.name = "umap_wnn",
  reduction.key = "UMAPwnn_",
  assay = "RNA",
  verbose = TRUE
)

# find clusters 
MultiOme_RNA_integrated <- FindClusters(MultiOme_RNA_integrated,
  graph.name = "wsnn",
  cluster.name = "wsnn_clusters",
  algorithm = algorithm,
  resolution = resolution,
  verbose = TRUE
)

# save combined/integrated Seurat object
saveRDS(MultiOme_RNA_integrated, file=paste0(rds_dir, project, "_combined_integrated.rds"))

u2 <- DimPlot(MultiOme_RNA_integrated, reduction = "rna_integrated", group.by = c("orig.ident", "rna_integrated_clusters"))
ggsave(paste0(integrated_dir, "all_RNA_integrated_UMAP.pdf"), plot = u2, width = 16, height = 8)

u4 <- DimPlot(MultiOme_ATAC_integrated, reduction = "atac_integrated", group.by = c("orig.ident", "atac_integrated_clusters"))
ggsave(paste0(integrated_dir, "all_ATAC_integrated_UMAP.pdf"), plot = u4, width = 16, height = 8)

u5 <- DimPlot(MultiOme_RNA_integrated, reduction = "umap_wnn", group.by = c("orig.ident", "wsnn_clusters"))
ggsave(paste0(integrated_dir, project, "_combined_integrated_UMAP.pdf"), plot = u5, width = 16, height = 8)

ggarrange(u1, u3, ncol = 1, nrow = 2)
ggsave(paste0(integrated_dir, project, "_RNA_ATAC_unintegrated_UMAP.pdf"), width = 16, height = 16)

ggarrange(u2, u4, u5, ncol = 1, nrow = 3)
ggsave(paste0(integrated_dir, project, "_RNA_ATAC_WNN_integrated_UMAP.pdf"), width = 16, height = 24)

# rename combined object so that it makes more sense
MultiOme_RNA_ATAC_integrated <- MultiOme_RNA_integrated

#### Process MGATK files ####

# load in filtered mgatk rds files and create list of mgatk objects
MGATK_rds_list <- c()

for (sample in sample_list) {
  assign(paste0("MultiOme_MGATK_", sample), readRDS(paste0(rds_dir, sample, "_MGATK_filtered.rds")))
  MGATK_rds_list <- append(MGATK_rds_list, get(paste0("MultiOme_MGATK_", sample)))
}

# remove first object from list, as this is the object that the others will be merged to
y <- MGATK_rds_list[-1]

# merge into a single Seurat object 
MultiOme_MGATK_merged <- merge(MGATK_rds_list[[1]],
                               y = y,
                               add.cell.ids = sample_list,
                               project = paste0(project, "_MGATK"),
                               merge.data = TRUE)

# join layers
MultiOme_MGATK_merged[["mito"]] <- JoinLayers(MultiOme_MGATK_merged[["mito"]])

# import mgatk data into main object
MultiOme_RNA_ATAC_integrated[["mito"]] <- MultiOme_MGATK_merged[["mito"]]

# import MGATK metadata
MGATK_meta <- MultiOme_MGATK_merged@meta.data
MultiOme_RNA_ATAC_integrated <- AddMetaData(MultiOme_RNA_ATAC_integrated, metadata = MGATK_meta)

# rename object and remove old version
MultiOme_RNA_ATAC_MGATK_integrated <- MultiOme_RNA_ATAC_integrated
rm(MultiOme_RNA_ATAC_integrated)
rm(MultiOme_MGATK_merged)

# set default assay as mito
DefaultAssay(MultiOme_RNA_ATAC_MGATK_integrated) = "mito"

# create new metadata column (mtDNA_count) and plot this metric
MultiOme_RNA_ATAC_MGATK_integrated@meta.data$mtDNA_count = colSums(MultiOme_RNA_ATAC_MGATK_integrated, slot = "counts")

##### Calculate heteroplasmy from mgatk data #####

# assign depth cutoff
depth_cutoff <- 0

# define the heterplasmic sites
het_positions <- c("1781", "1866", "3009", "3823", "5019", "5024", "13614", "13715", "15200", "16232")

# define the reference base for each site (in same order)
ref_bases <- c("C", "A", "G", "T", "A", "C", "C", "C", "A", "A")

# define the mutant base for each site (in same order)
mut_bases <- c("T", "G", "T", "C", "G", "T", "T", "T", "G", "T")

# define the column numbers where the base counts will be found 
# order is: 1 = A-fwd, 2 = C-fwd, 3 = T-fwd, 4 = G-fwd, 5 = A-rev, 6 = C-rev, 7 = T-rev, 8 = G-rev
ref_bases_fwd <- c(2, 1, 4, 3, 1, 2, 2, 2, 1, 1)
ref_bases_rev <- c(6, 5, 8, 7, 5, 6, 6, 6, 5, 5)
mut_bases_fwd <- c(3, 4, 3, 2, 4, 3, 3, 3, 4, 3)
mut_bases_rev <- c(7, 8, 7, 6, 8, 7, 7, 7, 8, 7)

#generate dataframes containing base calls, depth and heteroplasmy for each mutation site and add heteroplasmy/depth to metadata
for(pos in het_positions) {
  index <- match(pos, het_positions)
  ref_base <- ref_bases[index]
  mut_base <- mut_bases[index]
  grep_pos <- paste0("-", pos, "-")
  print(grep_pos)

  pos_df <- as.data.frame(t(as.data.frame(MultiOme_RNA_ATAC_MGATK_integrated[["mito"]]$counts[grep(grep_pos, rownames(MultiOme_RNA_ATAC_MGATK_integrated[["mito"]]$counts)),])))
  pos_df$Depth <- rowSums(pos_df)
  pos_df$Allele_ref <- pos_df[,ref_bases_fwd[index]] + pos_df[,ref_bases_rev[index]]
  pos_df$Allele_mut <- pos_df[,mut_bases_fwd[index]] + pos_df[,mut_bases_rev[index]]

  pos_df$Heteroplasmy <- pos_df$Allele_mut / (pos_df$Allele_mut + pos_df$Allele_ref)
  pos_df$Heteroplasmy <- ifelse(pos_df$Depth >= depth_cutoff, pos_df$Heteroplasmy, NA)
  assign(paste0("MultiOme_RNA_ATAC_MGATK_integrated@meta.data$Heteroplasmy_", pos), pos_df$Heteroplasmy)
  assign(paste0("MultiOme_RNA_ATAC_MGATK_integrated@meta.data$Depth_", pos), pos_df$Depth)

  assign(paste0("mito_", pos), pos_df)
  write.csv(pos_df, paste0(mgatk_dir, project, "_mtDNA_heteroplasmy_depth_results_", pos, ".csv"))

}

# set default assay as mito
DefaultAssay(MultiOme_RNA_ATAC_MGATK_integrated) <- "mito"

# calculate heteroplasmy for all 5024 variants
mito_all_5024 <- data.frame(cell= rownames(mito_5024),
                            Depth = (mito_1781$Depth + mito_3009$Depth + mito_5024$Depth + mito_13614$Depth + mito_13715$Depth + mito_3823$Depth + mito_1866$Depth),
                            Allele_ref = (mito_1781$Allele_mut + mito_3009$Allele_mut + mito_5024$Allele_ref + mito_13614$Allele_mut + mito_13715$Allele_ref + mito_3823$Allele_mut + mito_1866$Allele_mut),
                            Allele_mut = (mito_1781$Allele_ref + mito_3009$Allele_ref + mito_5024$Allele_mut + mito_13614$Allele_ref + mito_13715$Allele_mut + mito_3823$Allele_ref + mito_1866$Allele_ref))
mito_all_5024$Heteroplasmy <- mito_all_5024$Allele_mut / (mito_all_5024$Allele_mut + mito_all_5024$Allele_ref)
mito_all_5024$Heteroplasmy <- ifelse(mito_all_5024$Depth >= 0, mito_all_5024$Heteroplasmy, NA)

write.csv(mito_all_5024, paste0(mgatk_dir, project, "_mtDNA_heteroplasmy_depth_results_all_5024.csv"))

# add heteroplasmy and depth for all 5024 variants to metadata
MultiOme_RNA_ATAC_MGATK_integrated@meta.data$Heteroplasmy_all_5024 <- mito_all_5024$Heteroplasmy
MultiOme_RNA_ATAC_MGATK_integrated@meta.data$Depth_all_5024 <- mito_all_5024$Depth

# save file as RDS
saveRDS(MultiOme_RNA_ATAC_MGATK_integrated, file = paste0(rds_dir, project, "_RNA_ATAC_MGATK_combined.rds"))

