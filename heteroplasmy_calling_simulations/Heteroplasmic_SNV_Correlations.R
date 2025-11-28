library(Seurat)
library(tidyverse)

##############################################################################################
# Correlation of heteroplasmy calls at m.5024C>T vs other SNV sites (Extended Data Figure 4H #
##############################################################################################

# Read in Seurat object containing experimental data with mgatk assay "mito"
experimental_data <- readRDS(".../path_to_seurat_object.rds")

DefaultAssay(experimental_data) <- "mito"

# extract 5024 data and calculate depth and alleles
mito_5024 <- as.data.frame(t(as.data.frame(experimental_data[["mito"]]$counts[grep( "-5024-", rownames(experimental_data[["mito"]]$counts)),])))
mito_5024$Depth5024 <- rowSums(mito_5024)
mito_5024$allele_C <- mito_5024[,2] + mito_5024[,6]
mito_5024$allele_T <- mito_5024[,3] + mito_5024[,7]
mito_5024$Heteroplasmy5024T <- NA

# calculate heteroplasmy - 5024
for(row in 1:nrow(mito_5024)){
  if(mito_5024[row, 9] > 0) {
    mito_5024[row, 12] = mito_5024[row, 11] / (mito_5024[row, 11] + mito_5024[row, 10])
  } 
}

# extract 13715 data and calculate depth and alleles
mito_13715 <- as.data.frame(t(as.data.frame(experimental_data[["mito"]]$counts[grep( "-13715-", rownames(experimental_data[["mito"]]$counts)),])))
mito_13715$Depth13715 <- rowSums(mito_13715)
mito_13715$allele_C <- mito_13715[,2] + mito_13715[,6]
mito_13715$allele_T <- mito_13715[,3] + mito_13715[,7]
mito_13715$Heteroplasmy13715T <- NA

# calculate heteroplasmy - 13715
for(row in 1:nrow(mito_13715)){
  if(mito_13715[row, 9] > 0) {
    mito_13715[row, 12] = mito_13715[row, 11] / (mito_13715[row, 11] + mito_13715[row, 10])
  } 
}

# extract 13614 data and calculate depth and alleles
mito_13614 <- as.data.frame(t(as.data.frame(experimental_data[["mito"]]$counts[grep( "-13614-", rownames(experimental_data[["mito"]]$counts)),])))
mito_13614$Depth13614 <- rowSums(mito_13614)
mito_13614$allele_C <- mito_13614[,2] + mito_13614[,6]
mito_13614$allele_T <- mito_13614[,3] + mito_13614[,7]
mito_13614$Heteroplasmy13614T <- NA

# calculate heteroplasmy - 13614
for(row in 1:nrow(mito_13614)){
  if(mito_13614[row, 9] > 0) {
    mito_13614[row, 12] = mito_13614[row, 11] / (mito_13614[row, 11] + mito_13614[row, 10])
  } 
}

# extract 1781 data and calculate depth and alleles
mito_1781 <- as.data.frame(t(as.data.frame(experimental_data[["mito"]]$counts[grep( "-1781-", rownames(experimental_data[["mito"]]$counts)),])))
mito_1781$Depth1781 <- rowSums(mito_1781)
mito_1781$allele_C <- mito_1781[,2] + mito_1781[,6]
mito_1781$allele_T <- mito_1781[,3] + mito_1781[,7]
mito_1781$Heteroplasmy1781T <- NA

# calculate heteroplasmy - 1781
for(row in 1:nrow(mito_1781)){
  if(mito_1781[row, 9] > 0) {
    mito_1781[row, 12] = mito_1781[row, 11] / (mito_1781[row, 11] + mito_1781[row, 10])
  } 
}

# extract 1866 data and calculate depth and alleles
mito_1866 <- as.data.frame(t(as.data.frame(experimental_data[["mito"]]$counts[grep( "-1866-", rownames(experimental_data[["mito"]]$counts)),])))
mito_1866$Depth1866 <- rowSums(mito_1866)
mito_1866$allele_A <- mito_1866[,1] + mito_1866[,5]
mito_1866$allele_G <- mito_1866[,4] + mito_1866[,8]
mito_1866$Heteroplasmy1866G <- NA

# calculate heteroplasmy - 1866
for(row in 1:nrow(mito_1866)){
  if(mito_1866[row, 9] > 0) {
    mito_1866[row, 12] = mito_1866[row, 11] / (mito_1866[row, 11] + mito_1866[row, 10])
  } 
}

# extract 3009 data and calculate depth and alleles
mito_3009 <- as.data.frame(t(as.data.frame(experimental_data[["mito"]]$counts[grep( "-3009-", rownames(experimental_data[["mito"]]$counts)),])))
mito_3009$Depth3009 <- rowSums(mito_3009)
mito_3009$allele_G <- mito_3009[,4] + mito_3009[,8]
mito_3009$allele_T <- mito_3009[,3] + mito_3009[,7]
mito_3009$Heteroplasmy3009T <- NA

# calculate heteroplasmy - 3009
for(row in 1:nrow(mito_3009)){
  if(mito_3009[row, 9] > 0) {
    mito_3009[row, 12] = mito_3009[row, 11] / (mito_3009[row, 11] + mito_3009[row, 10])
  } 
}

# extract 3823 data and calculate depth and alleles
mito_3823 <- as.data.frame(t(as.data.frame(experimental_data[["mito"]]$counts[grep( "-3823-", rownames(experimental_data[["mito"]]$counts)),])))
mito_3823$Depth3823 <- rowSums(mito_3823)
mito_3823$allele_C <- mito_3823[,2] + mito_3823[,6]
mito_3823$allele_T <- mito_3823[,3] + mito_3823[,7]
mito_3823$Heteroplasmy3823C <- NA

# calculate heteroplasmy - 3823
for(row in 1:nrow(mito_3823)){
  if(mito_3823[row, 9] > 0) {
    mito_3823[row, 12] = mito_3823[row, 10] / (mito_3823[row, 11] + mito_3823[row, 10])
  } 
}

# convert rownames to column for all preceding individual mito datasets
mito_5024_corr <- tibble::rownames_to_column(mito_5024, "cell")
mito_13715_corr <- tibble::rownames_to_column(mito_13715, "cell")
mito_13614_corr <- tibble::rownames_to_column(mito_13614, "cell")
mito_1781_corr <- tibble::rownames_to_column(mito_1781, "cell")
mito_1866_corr <- tibble::rownames_to_column(mito_1866, "cell")
mito_3009_corr <- tibble::rownames_to_column(mito_3009, "cell")
mito_3823_corr <- tibble::rownames_to_column(mito_3823, "cell")

# create a data frame that contains all heteroplasmy values calculated above for each allele
# keep only columns of interest
df_list <- list(mito_5024_corr, mito_13715_corr, mito_13614_corr, mito_1781_corr, mito_1866_corr, mito_3009_corr, mito_3823_corr)
het_table <- df_list %>% purrr::reduce(full_join, by="cell")
het_table <- het_table %>% dplyr::select("cell", "Heteroplasmy5024T", "Heteroplasmy13715T", "Heteroplasmy13614T",
                "Heteroplasmy1781T", "Heteroplasmy1866G", "Heteroplasmy3009T", 
                "Heteroplasmy3823C", "Depth5024")
het_table <- het_table %>% na.omit()

# Create empty list for correlation results
het_results <- list()

# Calculate correlations for each SNV compared to m.5024C>T with depth cut off between 0-150 at increments of 5
for(depth in 1:31){
  vec <- numeric(8)
  min_coverage <- (depth -1) * 5
  vec[[1]] <- min_coverage
  het_table_filtered <- het_table %>% dplyr::filter(Depth5024 >= min_coverage)
  vec[[2]] <- (nrow(het_table_filtered) / nrow(het_table)) * 100
  vec[[3]] <- cor(het_table_filtered$Heteroplasmy5024T, het_table_filtered$Heteroplasmy13715T)
  vec[[4]] <- cor(het_table_filtered$Heteroplasmy5024T, het_table_filtered$Heteroplasmy13614T)
  vec[[5]] <- cor(het_table_filtered$Heteroplasmy5024T, het_table_filtered$Heteroplasmy1781T)
  vec[[6]] <- cor(het_table_filtered$Heteroplasmy5024T, het_table_filtered$Heteroplasmy1866G)
  vec[[7]] <- cor(het_table_filtered$Heteroplasmy5024T, het_table_filtered$Heteroplasmy3009T)
  vec[[8]] <- cor(het_table_filtered$Heteroplasmy5024T, het_table_filtered$Heteroplasmy3823C)
  het_results[[depth]] <- vec
}

het_results_df <- as.data.frame(do.call(rbind, het_results))
