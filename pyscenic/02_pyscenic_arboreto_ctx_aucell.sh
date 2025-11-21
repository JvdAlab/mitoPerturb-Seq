#!/bin/sh
#SBATCH -p cluster7
#SBATCH --time=96:00:00
#SBATCH -w sledge14
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -J pyscenic
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL,BEGIN,END

#SBATCH -o /../../slurm-%j.out
#SBATCH -e /../../slurm-%j.err

module load miniconda/2024-02-20
conda deactivate
conda deactivate
conda activate /../../.conda/envs/pyscenic

cd /../../

/../../.conda/envs/pyscenic/bin/arboreto_with_multiprocessing.py \
project_rna_assay_filtered.loom \
allTFs_hg38.txt \
--method grnboost2 \
--output project_rna_assay_filtered_adjacencies.csv \
--seed 777


db_names=("pyscenic_databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" "pyscenic_databases/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather" "pyscenic_databases/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" "pyscenic_databases/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather")

echo ${db_names[@]}

pyscenic ctx \
project_rna_assay_filtered_adjacencies.csv \
${db_names[@]} \
--annotations_fname pyscenic_databases/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname project_rna_assay_filtered.loom \
--output project_rna_assay_filtered_ctxreg.csv \
--mask_dropouts \
--num_workers 32 


pyscenic aucell \
project_rna_assay_filtered.loom \
project_rna_assay_filtered_ctxreg.csv \
--output project_rna_assay_pyscenic_output_filtered.loom \
--num_workers 32 
