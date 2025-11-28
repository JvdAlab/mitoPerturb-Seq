##### Guide detection 

# load packages
library(ComplexHeatmap)
library(cowplot)
library(dplyr)
library(fgsea)
library(ggalluvial)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(gprofiler2)
library(matrixStats)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(scales)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(Signac)
library(SingleCellExperiment)
library(SingleR)
library(viridis)

library(AnnotationDbi)
library(AnnotationHub)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(msigdbr)
library(org.Mm.eg.db)


#### Set Project Variables ####

# enter project name
Project <- "Sample1"

# set base directory
baseDir <- "/../../"
setwd(baseDir)


##### STEP ONE - Load Seurat object #####

# set data folder
rds_files <- "/../../rds_files/"

# load the Seurat object 
MultiOme_obj <- readRDS(paste0(rds_files,"Sample1_RNA_processed.rds"))
MultiOme_obj

# extract cell metadata 
cells_meta <- MultiOme_obj@meta.data
# create new column called cell_barcode (this is the rownames but without -1)
cells_meta$cell_barcode <- gsub( "-1", "", rownames(cells_meta))


##### STEP TWO - Guide detection #####

# set enrich data directory
enrich <- "/../../"

# enriched data 
Guides_sam_enrich <- read.csv(paste0(enrich,"Sample1_S1_L001_UMI_extract_NNNNNNNNNNNNNNNN_bwa_85bp_R2_aligned_guides_F4.sam"), sep="\t", header=FALSE)

head(Guides_sam_enrich, 5) 

# create new column called cell
# appended as last column
Guides_sam_enrich$cell <- gsub( ".*_", "", Guides_sam_enrich$V1)

# edit column names
colnames(Guides_sam_enrich)[1] <- "QNAME" # query template name
colnames(Guides_sam_enrich)[2] <- "flag" # bitwise flags
colnames(Guides_sam_enrich)[3] <- "guide" # reference sequence name, in this instance guide
colnames(Guides_sam_enrich)[4] <- "POS" # 1- based leftmost mapping position
colnames(Guides_sam_enrich)[5] <- "MAPQ" # mapping quality, 0-60 with 60 being high confidence
colnames(Guides_sam_enrich)[6] <- "CIGAR" # cigar string, use this to find end coordinate
colnames(Guides_sam_enrich)[7] <- "RNEXT" # ref name of the primary alignment of the next read 
colnames(Guides_sam_enrich)[8] <- "PNEXT" # position of the primary alignment of the next read
colnames(Guides_sam_enrich)[9] <- "TLEN" # observed template length, insert size if paired
colnames(Guides_sam_enrich)[10] <- "SEQ" # sequence, must equal sum of lengths of M/I/S/=/X in cigar

# note the number of non-assigned guides 
table(Guides_sam_enrich$guide) 

# order by MAPQ column
Guides_sam_enrich <- Guides_sam_enrich[order(Guides_sam_enrich$MAPQ, decreasing = TRUE),]

# retain only those with a MAPQ score greater than 30 
Guides_sam_enrich <- Guides_sam_enrich[Guides_sam_enrich$MAPQ >= 30,]

# remove those that have an NA
Guides_sam_enrich <- Guides_sam_enrich[!is.na(Guides_sam_enrich$guide),] 

# extract three columns of interest
Guides_df_enrich <- Guides_sam_enrich[,c("cell", "guide", "flag")]

# remove any guides that are 'empty' 
Guides_df_enrich <- Guides_df_enrich[Guides_df_enrich$guide != "",]

# add column called gene, this is the guide column minus the -qualifier, thereby collapsing guides
Guides_df_enrich$gene <-gsub( "-.*", "", Guides_df_enrich$guide)

# make guide matrix to extract top guide
Guides_mat_enrich <- as.data.frame(table(Guides_df_enrich[, c("cell", "guide")])) 

Guides_mat_enrich <- reshape2::dcast(Guides_mat_enrich, formula = cell ~ guide)
Guides_mat_enrich <- Guides_mat_enrich[,colnames(Guides_mat_enrich) != "mRFP"] 
rownames(Guides_mat_enrich) <- paste0(Guides_mat_enrich$cell, "-1")
Guides_mat_enrich <- Guides_mat_enrich[,colnames(Guides_mat_enrich) != "cell"]

# select only those cells in metadata
Guides_mat_filt_enrich <- Guides_mat_enrich[rownames(Guides_mat_enrich) %in% rownames(cells_meta),]
Guides_mat_filt_enrich[1:7,1:7]

Guides_mat_enrich2 <- Guides_mat_filt_enrich[rowSums(Guides_mat_filt_enrich) > 0,]

# max enrichment function
x <- c(NA, 0, 1, 0, 1, 2, 2,40,100)

which_max_ENRICHMENT <- function(x){
  # this function assigns a guide to each cell
  # input: x: a vector of occurrences of all guides in the given cell
  
  Max_indx <- which.max(x) # finds the position of the guide with highest frequency 
  Max_Val <- x[Max_indx] # value of guide with max frequency 
  
  SUM_X <- sum(x, na.rm = TRUE) # sum of all guide frequencies
  multi_guides <- Max_Val / SUM_X
  
  if(Max_Val < 2){ # alter this value to indicate minimum number of reads required to call specific guide
    return("unknown")
  } else if(multi_guides < 0.67){ # alter this value to indicate level of enrichment required to call specific guide
    return("multiple_guides")
  } else {
    return(names(Max_indx) ) # if you are confident enough, label as guide with max frequency
  }
}

# apply function - adds a column listing the results of the above function
Guides_mat_enrich2$top_guide <- apply(Guides_mat_enrich2, 1, which_max_ENRICHMENT)
Guides_mat_enrich2$cell <- rownames(Guides_mat_enrich2) # adds column with cell barcode
table(Guides_mat_enrich2$top_guide)


### multiple guides 

# creates a matrix containing cells that have multiple guides
guides_multiple_enrich <- Guides_mat_enrich2[Guides_mat_enrich2$top_guide == "multiple_guides",]

# set control_guides list and then remove those from the matrix (they will not be counted towards multiple guides)
control_guides <- c("Eomes-1", "Eomes-2", "Eomes-3", "Neurod1-2", "Neurod1-3", "Neurod1-4", "Nontargeting-1", "Nontargeting-2", "Nontargeting-3",
                    "Olig1-2", "Olig1-3", "Olig1-4", "Opn4-2", "Opn4-3", "Opn4-6", "Rgr-1", "Rgr-2", "Rgr-3", "Rrh-1", "Rrh-3", "Rrh-6")
guides_multiple_enrich <- guides_multiple_enrich[, !colnames(guides_multiple_enrich) %in% control_guides]

# re-run the enrichment function, excluding the final two columns  
guides_multiple_enrich$top_guide <- apply(guides_multiple_enrich[,-c(ncol(guides_multiple_enrich), ncol(guides_multiple_enrich)-1)], 1, which_max_ENRICHMENT)
# retain cell and guide column
guides_multiple_enrich <- guides_multiple_enrich[,c("cell", "top_guide")]
# remove those that are still assigned multiple guides
guides_multiple_enrich <- guides_multiple_enrich[guides_multiple_enrich$top_guide != "multiple_guides",]
# remove those that are unknown
guides_multiple_enrich <- guides_multiple_enrich[guides_multiple_enrich$top_guide != "unknown",]


# pull out top_guide
guide_anno_enrich <-data.frame(cell=rownames(Guides_mat_enrich2), top_guide = Guides_mat_enrich2$top_guide)
guide_anno_enrich <- na.omit(guide_anno_enrich)

# pull out top_guide from those that have multiple guides
guide_anno_enrich$top_guide_from_multiple <- guides_multiple_enrich[match(guide_anno_enrich$cell, guides_multiple_enrich$cell),]$top_guide
guide_anno_enrich$top_guide_final <- ifelse(!is.na(guide_anno_enrich$top_guide_from_multiple), guide_anno_enrich$top_guide_from_multiple, guide_anno_enrich$top_guide )

# add top_guide to metadata
MultiOme_obj@meta.data$top_guide_enrich <- guide_anno_enrich[match(rownames(MultiOme_obj@meta.data) , guide_anno_enrich$cell),]$top_guide
# add top_guide_final to metadata, this will take into account those with multiple guides
MultiOme_obj@meta.data$top_guide_final_enrich <- guide_anno_enrich[match(rownames(MultiOme_obj@meta.data) , guide_anno_enrich$cell),]$top_guide_final

# create 'gene' column - this is the guide column minus the -qualifier, thereby collapsing guides
MultiOme_obj@meta.data$top_guide_enrich_combined <- gsub( "-.*", "", MultiOme_obj@meta.data$top_guide_enrich)
MultiOme_obj@meta.data$top_guide_final_enrich_combined <- gsub( "-.*", "", MultiOme_obj@meta.data$top_guide_final_enrich)


##### STEP THREE - Save data  ##### 

# save file as RDS
saveRDS(MultiOme_obj, file="/../../Sample1_guide_detection/Sample1_GuideEnrich_Detection.rds")

