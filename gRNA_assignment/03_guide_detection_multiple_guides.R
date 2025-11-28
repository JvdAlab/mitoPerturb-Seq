##### Guide detection, multiple guides 

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
library(tidyverse)
library(viridis)

library(AnnotationDbi)
library(AnnotationHub)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(msigdbr)


#### Set Project Variables ####

# enter project name
project <- "project"

# results directory
baseDir <- "/../../"
setwd(baseDir)


##### STEP ONE - Load metadata #####

# metadata location
metadata_dir <- "/../../"

# metadata 
metadata <- read.csv(paste0(metadata_dir, "project_metadata.csv"))
head(metadata)

# add additional sample/cell id columns
metadata$sample_cell <- metadata$X
metadata <- metadata %>%
  separate(col = X,  
           into = c("sample", "cell_id"),
           sep = "_(?=[^_]+$)")  # separates at the last underscore
head(metadata)

# NB: this now contains a sample column which can be used for filtering

# create a subset of metadata for pilot1 and pilot2
metadata_sampleone <- metadata %>% dplyr::filter(sample == "Sample1")
rownames(metadata_sampleone) <- metadata_sampleone$cell_id

metadata_sampletwo <- metadata %>% dplyr::filter(sample == "Sample2")
rownames(metadata_sampletwo) <- metadata_sampletwo$cell_id


##### STEP TWO - Load and process Sample1 #####

# enrich data directory
enrichp1 <- "/../../"

# enriched data 
Guides_sam_enrichp1 <- read.csv(paste0(enrichp1,"Sample1_S1_L001_UMI_extract_NNNNNNNNNNNNNNNN_bwa_85bp_R2_aligned_guides_F4.sam"), sep="\t", header=FALSE)

head(Guides_sam_enrichp1, 5) 

# create new column called cell
# appended as last column
Guides_sam_enrichp1$cell <- gsub( ".*_", "", Guides_sam_enrichp1$V1)

# edit column names
colnames(Guides_sam_enrichp1)[1] <- "QNAME" # query template name
colnames(Guides_sam_enrichp1)[2] <- "flag" # bitwise flags
colnames(Guides_sam_enrichp1)[3] <- "guide" # reference sequence name, in this instance guide
colnames(Guides_sam_enrichp1)[4] <- "POS" # 1- based leftmost mapping position
colnames(Guides_sam_enrichp1)[5] <- "MAPQ" # mapping quality, 0-60 with 60 being high confidence
colnames(Guides_sam_enrichp1)[6] <- "CIGAR" # cigar string, use this to find end coordinate
colnames(Guides_sam_enrichp1)[7] <- "RNEXT" # ref name of the primary alignment of the next read 
colnames(Guides_sam_enrichp1)[8] <- "PNEXT" # position of the primary alignment of the next read
colnames(Guides_sam_enrichp1)[9] <- "TLEN" # observed template length, insert size if paired
colnames(Guides_sam_enrichp1)[10] <- "SEQ" # sequence, must equal sum of lengths of M/I/S/=/X in cigar

# note the number of non-assigned guides 
table(Guides_sam_enrichp1$guide) 

# order by MAPQ column
Guides_sam_enrichp1 <- Guides_sam_enrichp1[order(Guides_sam_enrichp1$MAPQ, decreasing = TRUE),]

# retain only those with a MAPQ score greater than 30 
Guides_sam_enrichp1 <- Guides_sam_enrichp1[Guides_sam_enrichp1$MAPQ >= 30,]

# remove those that have an NA
Guides_sam_enrichp1 <- Guides_sam_enrichp1[!is.na(Guides_sam_enrichp1$guide),] 

# extract three columns of interest
Guides_df_enrichp1 <- Guides_sam_enrichp1[,c("cell", "guide", "flag")]

# remove any guides that are 'empty' 
Guides_df_enrichp1 <- Guides_df_enrichp1[Guides_df_enrichp1$guide != "",]

# add column called gene, this is the guide column minus the -qualifier, thereby collapsing guides
Guides_df_enrichp1$gene <-gsub( "-.*", "", Guides_df_enrichp1$guide)

# make guide matrix to extract top guide
Guides_mat_enrichp1 <- as.data.frame(table(Guides_df_enrichp1[, c("cell", "guide")])) 

Guides_mat_enrichp1 <- reshape2::dcast(Guides_mat_enrichp1, formula = cell ~ guide)
Guides_mat_enrichp1 <- Guides_mat_enrichp1[,colnames(Guides_mat_enrichp1) != "mRFP"] 
rownames(Guides_mat_enrichp1) <- paste0(Guides_mat_enrichp1$cell, "-1")
Guides_mat_enrichp1 <- Guides_mat_enrichp1[,colnames(Guides_mat_enrichp1) != "cell"]

# select only those cells in metadata
Guides_mat_filt_enrichp1 <- Guides_mat_enrichp1[rownames(Guides_mat_enrichp1) %in% rownames(metadata_sampleone),]
Guides_mat_filt_enrichp1[1:7,1:7]


Guides_mat_enrich2p1 <- Guides_mat_filt_enrichp1[rowSums(Guides_mat_filt_enrichp1) > 0,]

# max enrichment function
x <- c(NA, 0, 1, 0, 1, 2, 2,40,100)
# this will list multiple guides and so all cells with a guide can removed
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
Guides_mat_enrich2p1$top_guide <- apply(Guides_mat_enrich2p1, 1, which_max_ENRICHMENT)
Guides_mat_enrich2p1$cell <- rownames(Guides_mat_enrich2p1) # adds column with cell barcode
table(Guides_mat_enrich2p1$top_guide)


### multiple guides 

# creates a matrix containing cells that have multiple guides
guides_multiple_enrichp1 <- Guides_mat_enrich2p1[Guides_mat_enrich2p1$top_guide == "multiple_guides",]

# set control_guides list and then remove those from the matrix (they will not be counted towards multiple guides)
control_guides <- c("Eomes-1", "Eomes-2", "Eomes-3", "Neurod1-2", "Neurod1-3", "Neurod1-4", "Nontargeting-1", "Nontargeting-2", "Nontargeting-3",
                    "Olig1-2", "Olig1-3", "Olig1-4", "Opn4-2", "Opn4-3", "Opn4-6", "Rgr-1", "Rgr-2", "Rgr-3", "Rrh-1", "Rrh-3", "Rrh-6")
guides_multiple_enrichp1 <- guides_multiple_enrichp1[, !colnames(guides_multiple_enrichp1) %in% control_guides]

# edited function so that it will return the names of all guides present, if multiple
# these will be a list guide,guide,guide,etc
which_max_ENRICHMENT <- function(x){
  Max_indx <- which.max(x)
  Max_Val <- x[Max_indx]
  
  SUM_X <- sum(x, na.rm = TRUE)
  multi_guides <- Max_Val / SUM_X
  
  if (Max_Val < 2) {
    return("unknown")
  } else if (multi_guides < 0.67) {
    # Return all guide names with non-zero values
    return(paste(names(x)[x > 0], collapse = ","))
  } else {
    return(names(x)[Max_indx])
  }
}

# re-run the enrichment function, excluding the final two columns  
# instead of multiple_guides, top_guide now contains which guides are observed
guides_multiple_enrichp1$top_guide <- apply(guides_multiple_enrichp1[,-c(ncol(guides_multiple_enrichp1), ncol(guides_multiple_enrichp1)-1)], 1, which_max_ENRICHMENT)
# remove those that are unknown
guides_multiple_enrichp1 <- guides_multiple_enrichp1[guides_multiple_enrichp1$top_guide != "unknown",]

# add dataset id so that the pilot 1 data can be merged with the pilot 2 data
guides_multiple_enrichp1 <- guides_multiple_enrichp1 %>%
  mutate(cell = paste0("sampleone_", cell))


##### STEP THREE - Load and process Sample2 #####

# enrich data directory
enrichp2 <- "/../../"

# enriched data 
Guides_sam_enrichp2 <- read.csv(paste0(enrichp2,"Sample2_S1_L001_UMI_extract_NNNNNNNNNNNNNNNN_bwa_85bp_R2_aligned_guides_F4.sam"), sep="\t", header=FALSE)

head(Guides_sam_enrichp2, 5) 

# create new column called cell
# appended as last column
Guides_sam_enrichp2$cell <- gsub( ".*_", "", Guides_sam_enrichp2$V1)

# edit column names
colnames(Guides_sam_enrichp2)[1] <- "QNAME" # query template name
colnames(Guides_sam_enrichp2)[2] <- "flag" # bitwise flags
colnames(Guides_sam_enrichp2)[3] <- "guide" # reference sequence name, in this instance guide
colnames(Guides_sam_enrichp2)[4] <- "POS" # 1- based leftmost mapping position
colnames(Guides_sam_enrichp2)[5] <- "MAPQ" # mapping quality, 0-60 with 60 being high confidence
colnames(Guides_sam_enrichp2)[6] <- "CIGAR" # cigar string, use this to find end coordinate
colnames(Guides_sam_enrichp2)[7] <- "RNEXT" # ref name of the primary alignment of the next read 
colnames(Guides_sam_enrichp2)[8] <- "PNEXT" # position of the primary alignment of the next read
colnames(Guides_sam_enrichp2)[9] <- "TLEN" # observed template length, insert size if paired
colnames(Guides_sam_enrichp2)[10] <- "SEQ" # sequence, must equal sum of lengths of M/I/S/=/X in cigar

# note the number of non-assigned guides 
table(Guides_sam_enrichp2$guide) 

# order by MAPQ column
Guides_sam_enrichp2 <- Guides_sam_enrichp2[order(Guides_sam_enrichp2$MAPQ, decreasing = TRUE),]

# retain only those with a MAPQ score greater than 30 
Guides_sam_enrichp2 <- Guides_sam_enrichp2[Guides_sam_enrichp2$MAPQ >= 30,]

# remove those that have an NA
Guides_sam_enrichp2 <- Guides_sam_enrichp2[!is.na(Guides_sam_enrichp2$guide),] 

# extract three columns of interest
Guides_df_enrichp2 <- Guides_sam_enrichp2[,c("cell", "guide", "flag")]

# remove any guides that are 'empty' 
Guides_df_enrichp2 <- Guides_df_enrichp2[Guides_df_enrichp2$guide != "",]

# add column called gene, this is the guide column minus the -qualifier, thereby collapsing guides
Guides_df_enrichp2$gene <-gsub( "-.*", "", Guides_df_enrichp2$guide)

# make guide matrix to extract top guide
Guides_mat_enrichp2 <- as.data.frame(table(Guides_df_enrichp2[, c("cell", "guide")])) 

Guides_mat_enrichp2 <- reshape2::dcast(Guides_mat_enrichp2, formula = cell ~ guide)
Guides_mat_enrichp2 <- Guides_mat_enrichp2[,colnames(Guides_mat_enrichp2) != "mRFP"] 
rownames(Guides_mat_enrichp2) <- paste0(Guides_mat_enrichp2$cell, "-1")
Guides_mat_enrichp2 <- Guides_mat_enrichp2[,colnames(Guides_mat_enrichp2) != "cell"]

# select only those cells in metadata
Guides_mat_filt_enrichp2 <- Guides_mat_enrichp2[rownames(Guides_mat_enrichp2) %in% rownames(metadata_sampletwo),]
Guides_mat_filt_enrichp2[1:7,1:7]

Guides_mat_enrich2p2 <- Guides_mat_filt_enrichp2[rowSums(Guides_mat_filt_enrichp2) > 0,]

# max enrichment function
x <- c(NA, 0, 1, 0, 1, 2, 2,40,100)
# this will list multiple guides and so all cells with a guide can removed
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
Guides_mat_enrich2p2$top_guide <- apply(Guides_mat_enrich2p2, 1, which_max_ENRICHMENT)
Guides_mat_enrich2p2$cell <- rownames(Guides_mat_enrich2p2) # adds column with cell barcode
table(Guides_mat_enrich2p2$top_guide)

### multiple guides 

# creates a matrix containing cells that have multiple guides
guides_multiple_enrichp2 <- Guides_mat_enrich2p2[Guides_mat_enrich2p2$top_guide == "multiple_guides",]

# set control_guides list and then remove those from the matrix (they will not be counted towards multiple guides)
control_guides <- c("Eomes-1", "Eomes-2", "Eomes-3", "Neurod1-2", "Neurod1-3", "Neurod1-4", "Nontargeting-1", "Nontargeting-2", "Nontargeting-3",
                    "Olig1-2", "Olig1-3", "Olig1-4", "Opn4-2", "Opn4-3", "Opn4-6", "Rgr-1", "Rgr-2", "Rgr-3", "Rrh-1", "Rrh-3", "Rrh-6")
guides_multiple_enrichp2 <- guides_multiple_enrichp2[, !colnames(guides_multiple_enrichp2) %in% control_guides]

# edited function so that it will return the names of all guides present, if multiple
# these will be a list guide,guide,guide,etc
which_max_ENRICHMENT <- function(x){
  Max_indx <- which.max(x)
  Max_Val <- x[Max_indx]
  
  SUM_X <- sum(x, na.rm = TRUE)
  multi_guides <- Max_Val / SUM_X
  
  if (Max_Val < 2) {
    return("unknown")
  } else if (multi_guides < 0.67) {
    # Return all guide names with non-zero values
    return(paste(names(x)[x > 0], collapse = ","))
  } else {
    return(names(x)[Max_indx])
  }
}

# re-run the enrichment function, excluding the final two columns  
# instead of multiple_guides, top_guide now contains which guides are observed
guides_multiple_enrichp2$top_guide <- apply(guides_multiple_enrichp2[,-c(ncol(guides_multiple_enrichp2), ncol(guides_multiple_enrichp2)-1)], 1, which_max_ENRICHMENT)
# remove those that are unknown
guides_multiple_enrichp2 <- guides_multiple_enrichp2[guides_multiple_enrichp2$top_guide != "unknown",]

# add dataset id so that the pilot 1 data can be merged with the pilot 2 data
guides_multiple_enrichp2 <- guides_multiple_enrichp2 %>%
  mutate(cell = paste0("sampletwo_", cell))


##### STEP FOUR - Join datasets together and filter #####

head(guides_multiple_enrichp1)
head(guides_multiple_enrichp2)

merged <- bind_rows(guides_multiple_enrichp1, guides_multiple_enrichp2)
head(merged)

###

# filter to keep only rows that contain multiple guides
# this removes all rows with only one guide
merged_multiple_only <- merged %>%
  dplyr::filter(str_detect(top_guide, ","))

# start from wide table of guide counts: merged_multiple_only
merged_multiple_countpercent_gene <- merged_multiple_only %>%
  # cell column needs to be first - this is the new cell column, includes sample ids
  relocate(cell) %>%
  
  # convert from wide (guide counts) to long format
  pivot_longer(
    cols = where(is.numeric),
    names_to = "guide",
    values_to = "count"
  ) %>%
  
  # collapse guides to gene names (remove "-guideNumber")
  mutate(gene = sub("-.*", "", guide)) %>%
  
  # summarise at the gene level: sum counts for all guides of the same gene
  group_by(cell, gene) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  
  # recalculate percentages from the gene-level counts
  group_by(cell) %>%
  mutate(percent = round(100 * count / sum(count), 1)) %>%
  arrange(cell, desc(percent), desc(count)) %>%
  ungroup() %>%
  
  # apply filtering: keep cells where the top gene % ≥ 20
  group_by(cell) %>%
  dplyr::filter(max(percent, na.rm = TRUE) >= 20) %>%
  ungroup() %>%
  
  # remove single-read genes (count = 1)
  dplyr::filter(count > 1) %>%
  
  # remove cells with only 1 gene after filtering
  group_by(cell) %>%
  dplyr::filter(n() > 1) %>%
  ungroup()

# rename Park2 = Prkn
merged_multiple_countpercent_gene <- merged_multiple_countpercent_gene %>%
  mutate(gene = ifelse(gene == "Park2", "Prkn", gene))

# remove guides that make up <10% of total reads per cell
merged_multiple_countpercent_gene <- merged_multiple_countpercent_gene %>% dplyr::filter(percent >= 10)


##### STEP FIVE - Grid of multiple guides, with counts #####

head(merged_multiple_countpercent_gene)

# create unique gene pairs per cell
gene_pairs <- merged_multiple_countpercent_gene %>%
  group_by(cell) %>%
  summarise(
    combos = list({
      genes <- sort(unique(gene))
      if (length(genes) > 1) {
        combn(genes, 2, simplify = TRUE) %>%      # always returns a matrix
          t() %>%                                 # each row = one pair
          as.data.frame(stringsAsFactors = FALSE) %>%
          setNames(c("gene1", "gene2"))
      } else {
        NULL
      }
    }),
    .groups = "drop"
  ) %>%
  dplyr::filter(!sapply(combos, is.null)) %>%   # keep only those with pairs
  unnest(combos) %>% 
  group_by(gene1, gene2) %>%
  summarise(n_cells = n(), .groups = "drop")

# Get all genes
all_genes <- sort(unique(c(gene_pairs$gene1, gene_pairs$gene2)))

# Create an empty matrix
gene_mat <- matrix(0, nrow = length(all_genes), ncol = length(all_genes),
                   dimnames = list(all_genes, all_genes))

# Fill in counts from gene_pairs
for(i in 1:nrow(gene_pairs)){
  g1 <- gene_pairs$gene1[i]
  g2 <- gene_pairs$gene2[i]
  n <- gene_pairs$n_cells[i]
  
  gene_mat[g1, g2] <- n
  gene_mat[g2, g1] <- n  # make it symmetric
}

# cluster to get row/col order
dist_mat <- dist(gene_mat)
clust <- hclust(dist_mat)
ord <- clust$order
mat_ord <- gene_mat[ord, ord, drop = FALSE]

# melt and keep only the upper triangle
df <- melt(mat_ord, varnames = c("row", "col"), value.name = "value")
df$keep <- upper.tri(mat_ord, diag = TRUE)[cbind(as.numeric(df$row), as.numeric(df$col))]
df <- df[df$keep, ]

ref_median <- 92  # median of gene numbers

# plot with ggplot (triangular heatmap)
ggplot(df, aes(x = row, y = col, fill = value)) +
  geom_tile(color = "grey80") +
  geom_text(aes(label = sprintf("%.0f", value)), size = 2) +
  scale_fill_gradientn(
    colors = c("white", "yellow", "red"),  # your 3-color scale
    limits = c(0, ref_median),             # fix scale 0 → 300
    oob = squish                           # values >300 shown as red
  ) +
  coord_equal() +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        panel.grid = element_blank()
  )

