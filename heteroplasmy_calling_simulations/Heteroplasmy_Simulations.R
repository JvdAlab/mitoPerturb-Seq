library(tidyverse)
library(Seurat)

# Based on experimental data - copy number accounts for the fact that there are 7 SNVs on each mtDNA genome, so the theoretical
# total number of bases at which heteroplasmy can be measured is the copy number multiplied by 7
copy_number <- 1750 * 7
het_mean <- 0.58
het_sd <- 0.11
depth_cutoff <- 20

####################################################################################
# Simulation of sampling heteroplasmy at different depths (Extended Data Figure 2B #
####################################################################################

# Set sampling depth
depth = 100

# Generate normally distributed set of heteroplasmy values for 2000 cells
simulation_df <- as.data.frame(rnorm(2000, het_mean, het_sd))

# Add sampling depth
simulation_df$X <- depth
colnames(simulation_df) <- c("het", "depth")

# Add empty column for heteroplasmy calls
simulation_df$call <- NA

for(cell in 1:nrow(simulation_df)) {
  mut <- abs(round(simulation_df[cell, 1] * copy_number)) # Calculates number of mutant alleles
  copies <- vector("list", copy_number) # Creates empty list with a position for each mtDNA copy

  # These for loops populate the "copies" vector with the correct number of wt and mut copies to sample from
  for(pos in 1:mut){
    copies[[pos]] = "mut"
  }
  for(pos in (mut1+1):copy_number) {
    copies[[pos]] = "wt"
  }
  
  sample_size <- simulation_df[cell, 2] # Sample size (set by "depth" variable)
  sample <- sample(copies, sample_size) # Randomly sample from the total mtDNA population based on the sample size
  het <- sum(sample == "mut") / sample_size # Calculate the heteroplasmy of the sampled molecules
  simulation_df[cell, 3] <- het # Add result to dataframe
}

########################################################################################
# Simulation of heteroplasmy calls in MitoPerturb-Seq dataset (Extended Data Figure 2C #
########################################################################################

# Read in Seurat object containing experimental data
experimental_data <- readRDS(".../path_to_seurat_object.rds")
meta_data <- experimental_data@meta.data

# Extract coverage depth data at heteroplasmic site(s) from relevant column in meta data and filter out any cells with 0 coverage
meta_depths <- as.data.frame(meta_data$mtDNA_depth)
colnames(meta_depths) <- "depth"
meta_depths <- meta_depths %>% dplyr::filter(depth > 0)

# Generate normally distributed set of heteroplasmy values for the correct number of cells
simulation_df <- as.data.frame(rnorm(nrow(meta_depths), het_mean, het_sd))

# Combine depth and heteroplasmy values and add an empty column to which heteroplasmy calls will be added
simulation_df <- merge(simulation_df, meta_depths, by = "row.names")
colnames(simulation_df) <- c("row_name", "het", "depth")
simulation_df <- simulation_df %>% dplyr::select(het, depth)
simulation_df$call <- NA

# Generate heteroplasmy calls based on sampling from the total mtDNA population at the corresponding mtDNA coverage depth
for(cell in 1:nrow(simulation_df)) {
  mut <- abs(round(simulation_df[cell, 1] * copy_number)) # Calculates number of mutant alleles
  copies <- vector("list", copy_number) # Creates empty list with a position for each mtDNA copy
  
  # These for loops populate the "copies" vector with the correct number of wt and mut copies to sample from
  for(pos in 1:mut){
    copies[[pos]] = "mut"
  }
  for(pos in (mut1+1):copy_number) {
    copies[[pos]] = "wt"
  }
  
  sample_size <- simulation_df[cell, 2] # Sample size is the coverage depth for the corresponding cell
  sample <- sample(copies, sample_size) # Randomly sample from the total mtDNA population based on the sample size
  het <- sum(sample == "mut") / sample_size # Calculate the heteroplasmy of the sampled molecules
  simulation_df[cell, 3] <- het # Add result to dataframe
}

# Colour based on depth cutoff and split into separate dataframes
simulation_df <- simulation_df %>% arrange(depth)
simulation_df$colour <- NA
simulation_df <- simulation_df %>% mutate(colour = case_when(depth < depth_cutoff ~ "red", depth >= depth_cutoff ~ "black"))

simulation_depth_over_cutoff <- simulation_df %>% dplyr::filter(depth >= depth_cutoff)
simulation_depth_under_cutoff <- simulation_df %>% dplyr::filter(depth < depth_cutoff)


