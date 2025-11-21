##### tricycle analysis

# load packages
library(cowplot)
library(data.table)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(RColorBrewer)
library(scCustomize)
library(Seurat)
library(SeuratWrappers)
library(Signac)
library(tidyverse)
library(tricycle)

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
MultiOme_obj_int <- readRDS(paste0(integrate,"project_RNA_ATAC_MGATK_GUIDES_combined.rds"))
MultiOme_obj_int

# create guide column
MultiOme_obj_int@meta.data$guide <- MultiOme_obj_int@meta.data$top_guide_final_enrich
MultiOme_obj_int@meta.data$guide <- ifelse(MultiOme_obj_int@meta.data$top_guide_final_enrich %in% c("Eomes-1", "Eomes-2", "Eomes-3", "Neurod1-2", "Neurod1-3", "Neurod1-4", 
                                                                                                    "Nontargeting-1", "Nontargeting-2", "Nontargeting-3", "Olig1-2", "Olig1-3", 
                                                                                                    "Olig1-4", "Opn4-2", "Opn4-3", "Opn4-6", "Rgr-1", "Rgr-2", "Rgr-3", "Rrh-1", 
                                                                                                    "Rrh-3", "Rrh-6"), "NT", MultiOme_obj_int@meta.data$top_guide_final_enrich)

# create gene column
MultiOme_obj_int@meta.data$gene <- gsub( "-.*", "", MultiOme_obj_int@meta.data$guide)

# create chosen_control column grouping NT/controls/10 genes
MultiOme_obj_int@meta.data$chosen_controls <- ifelse(MultiOme_obj_int@meta.data$gene %in% c("NT", "Akap1", "Atg5", "Dnm1l", "Mtfp1", "Nnt", "Oma1", "Park2", "Pink1", "Ppargc1a", "Snx9"), "nt_plus_ten_gene", 
                                                     ifelse(MultiOme_obj_int@meta.data$gene %in% c("multiple_guides", "unknown"), "unkmul", "opa1_polg_tfam")) 


# set results folder 
tricycle_folder <- "/../../tricycle/"


##### STEP TWO - Infer cell cycle position #####

# set default assay
DefaultAssay(MultiOme_obj_int) <- "RNA"

# need to join layers in order to use RNA assay 
MultiOme_obj_int[["RNA"]] <- JoinLayers(MultiOme_obj_int[["RNA"]])

# run tricycle - this function runs 'estimate_cycle_position' 
MultiOme_obj_int <- Runtricycle(object = MultiOme_obj_int, 
                                slot = "data", 
                                reduction.name = "tricycleEmbedding", 
                                reduction.key = "tricycleEmbedding_", 
                                gname = NULL, 
                                gname.type = "SYMBOL", # make sure this matches current gene names 
                                species = "mouse")

# create a meta data column for 0-2 pi
MultiOme_obj_int@meta.data$pi <- MultiOme_obj_int@meta.data$tricyclePosition/3.14
 

##### STEP THREE - Assess performance #####

# example one gene - Top2a
fit_top2a <- fit_periodic_loess(MultiOme_obj_int$tricyclePosition, 
                                plot.df$Top2a, 
                                plot = TRUE, 
                                x_lab = "Cell cycle position θ", 
                                y_lab = "log2(Top2a)", 
                                fig.title = paste0("Expression of Top2a along θ (n=", ncol(MultiOme_obj_int),")"))

fit_top2a$fig + theme_bw(base_size = 14)


# example multiple genes
# extract the expression level of chosen genes, tricyclePosition and pi
genes_to_plot <- FetchData(object = MultiOme_obj_int, vars = c("gene", "chosen_controls", "tricyclePosition", "pi", 
                                                               "Anln", "Aurka", "Bub1b", "Ccna2", "Ccnb1", "Ccnb2", "Ccnd1", "Ccne1", "Ccne2", "Ccnf", "Cdc6", "Cdc20", 
                                                               "Cdk1", "Cdk6", "Cdkn1a", "Cdkn3", "Cdt1", "Cenpe", "Dhfr", "E2f1", "Gins2", "Kpna2", "Mcm6", "Melk", 
                                                               "Nusap1", "Orc1", "Pbk", "Pcna", "Plk1", "Rfc4", "Rpa1", "Rrm2", "Smc2", "Tk1", "Top2a", "Tpx2"))

# save data 
write.csv(genes_to_plot, paste(tricycle_folder, project, "_tricycle_genes_to_plot_dataframe.csv", sep = ""), row.names = TRUE)

# define a list of columns to pivot
columns_to_pivot <- c("Anln", "Aurka", "Bub1b", "Ccna2", "Ccnb1", "Ccnb2", "Ccnd1", "Ccne1", "Ccne2", "Ccnf", "Cdc6", "Cdc20", 
                      "Cdk1", "Cdk6", "Cdkn1a", "Cdkn3", "Cdt1", "Cenpe", "Dhfr", "E2f1", "Gins2", "Kpna2", "Mcm6", "Melk", 
                      "Nusap1", "Orc1", "Pbk", "Pcna", "Plk1", "Rfc4", "Rpa1", "Rrm2", "Smc2", "Tk1", "Top2a", "Tpx2")

# pivot the data so that all genes can be plotted on the same figure
genes_to_plot_long <- genes_to_plot %>% pivot_longer(cols = all_of(columns_to_pivot), 
                                                     names_to = "Genes", 
                                                     values_to = "Expression")

# add additional column 
genes_to_plot_long$Phase <- ifelse(genes_to_plot_long$Genes %in% c("Ccnd1", "Cdk6", "Cdkn1a"), "G1",
                                   ifelse(genes_to_plot_long$Genes %in% c("Cdc6", "Cdt1", "E2f1", "Gins2", "Mcm6", "Orc1"), "G1/S", 
                                          ifelse(genes_to_plot_long$Genes %in% c("Ccne1", "Ccne2", "Dhfr", "Pcna", "Rfc4", "Rpa1", "Rrm2", "Tk1"), "S",
                                                 ifelse(genes_to_plot_long$Genes %in% c("Anln", "Ccnf", "Cdk1", "Kpna2", "Melk", "Pbk", "Smc2", "Top2a"), "G2", 
                                                        ifelse(genes_to_plot_long$Genes %in% c("Aurka", "Bub1b", "Ccna2", "Ccnb1", "Cenpe", "Nusap1", "Plk1", "Tpx2"), "G2/M", "M"))))) # M = Ccnb2, Cdc20, Cdkn3

phase_cols_indiv_genes <- c("Ccnd1" = "#009E73", "Cdk6" = "#009E73", "Cdkn1a" = "#009E73", 
                            "Cdc6" = "#56B4E9", "Cdt1" = "#56B4E9", "E2f1" = "#56B4E9", "Gins2" = "#56B4E9", "Mcm6" = "#56B4E9", "Orc1" = "#56B4E9",
                            "Ccne1" = "#0072B2", "Ccne2"= "#0072B2", "Dhfr"= "#0072B2", "Pcna"= "#0072B2", "Rfc4"= "#0072B2", "Rpa1"= "#0072B2", "Rrm2"= "#0072B2", "Tk1"= "#0072B2",
                            "Anln" = "#E69F00", "Ccnf" = "#E69F00", "Cdk1" = "#E69F00", "Kpna2" = "#E69F00", "Melk" = "#E69F00", "Pbk" = "#E69F00", "Smc2" = "#E69F00", "Top2a" = "#E69F00",
                            "Aurka" = "#D55E00", "Bub1b" = "#D55E00", "Ccna2" = "#D55E00", "Ccnb1" = "#D55E00", "Cenpe" = "#D55E00", "Nusap1" = "#D55E00", "Plk1" = "#D55E00", "Tpx2" = "#D55E00",
                            "Ccnb2" = "#F0E442", "Cdc20" = "#F0E442", "Cdkn3" = "#F0E442")

# set order of genes
genes_to_plot_long <- genes_to_plot_long %>% 
  mutate(Genes=factor(Genes)) %>%
  mutate(Genes=fct_relevel(Genes,c("Ccnd1", "Cdk6", "Cdkn1a", 
                                   "Cdc6", "Cdt1", "E2f1", "Gins2", "Mcm6", "Orc1",
                                   "Ccne1", "Ccne2", "Dhfr", "Pcna", "Rfc4", "Rpa1", "Rrm2", "Tk1",
                                   "Anln", "Ccnf", "Cdk1", "Kpna2", "Melk", "Pbk", "Smc2", "Top2a",
                                   "Aurka", "Bub1b", "Ccna2", "Ccnb1", "Cenpe", "Nusap1", "Plk1", "Tpx2",
                                   "Ccnb2", "Cdc20", "Cdkn3"))) %>%
  arrange(Genes)

# set order of phases
genes_to_plot_long <- genes_to_plot_long %>% 
  mutate(Phase=factor(Phase)) %>%
  mutate(Phase=fct_relevel(Phase,c("G1", "G1/S", "S", "G2", "G2/M", "M"))) %>%
  arrange(Phase)

dplyr::filter(genes_to_plot_long, chosen_controls == "nt_plus_ten_gene") %>% # extract category of interest
  ggplot(aes(x = pi, y = Expression, colour = Genes)) + 
  geom_smooth(method = "loess", se = FALSE) +
  scale_x_continuous(limits = c(0,2),
                     breaks = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2),
                     labels = c("0", "0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75", "2")) +
  facet_wrap(~ Genes, scales = "free_y") + # plot by gene
  scale_colour_manual(values = phase_cols_indiv_genes) +
  geom_vline(xintercept = 0.5, linetype = "dashed") + # amended start of S
  geom_vline(xintercept = 1, linetype = "dashed") + # pi start of G2M
  geom_vline(xintercept = 1.5, linetype = "dashed") + # 1.5 pi middle of M
  geom_vline(xintercept = 1.75, linetype = "dashed") + # start of G1
  labs(title = "", 
       x = "Cell Cycle Position (radians)", 
       y = "Gene Expression") +
  theme_minimal() + NoLegend()

