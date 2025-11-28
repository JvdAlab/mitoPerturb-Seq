#!/bin/bash

# Generate Combined Heatmap for Multiple DEG Sets

source "../config/config.sh"

# Script-specific parameters
TSS_UPSTREAM=5000
TSS_DOWNSTREAM=5000

# Define input and output
BW_FILE="${WORK_DIR}/ATF4_rpm-norm_avg.bw"
DEG_DIR="${BASE_DIR}/deg_genelists"
HEATMAP_DIR="${DEG_DIR}/deeptools_heatmaps"
mkdir -p "${HEATMAP_DIR}"

# Define DEG gene sets
declare -a DEG_NAMES=("POLG" "OPA1" "TFAM")

# Build file paths
combined_bed_files=()
combined_labels=()

for i in "${!DEG_NAMES[@]}"; do
    name="${DEG_NAMES[$i]}"
    overlap_file="${DEG_DIR}/${name}_ATF4_overlap.bed"
    deg_file="${DEG_DIR}/${name}_genes.bed"

    if [[ -f "${overlap_file}" ]]; then
        combined_bed_files+=("${overlap_file}")
        combined_labels+=("${name}_overlap")
    fi

    if [[ -f "${deg_file}" ]]; then
        combined_bed_files+=("${deg_file}")
        combined_labels+=("${name}_DEG")
    fi
done

sample_name=$(basename "${BW_FILE}" | cut -d '-' -f1)

sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=combined_heatmap
#SBATCH --output=${HEATMAP_DIR}/combined_heatmap_%j.out
#SBATCH --error=${HEATMAP_DIR}/combined_heatmap_%j.err
#SBATCH --partition=<PARTITION>
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Compute matrix
computeMatrix reference-point \
    --scoreFileName "${BW_FILE}" \
    --regionsFileName ${combined_bed_files[@]} \
    --beforeRegionStartLength ${TSS_UPSTREAM} \
    --afterRegionStartLength ${TSS_DOWNSTREAM} \
    --referencePoint TSS \
    --numberOfProcessors 4 \
    --skipZeros \
    --outFileName "${HEATMAP_DIR}/combined_genes_matrix.gz"

# Plot combined heatmap (SVG)
plotHeatmap \
    --matrixFile "${HEATMAP_DIR}/combined_genes_matrix.gz" \
    --outFileName "${HEATMAP_DIR}/combined_genes_heatmap.svg" \
    --outFileNameMatrix "${HEATMAP_DIR}/combined_genes_matrix.tab" \
    --regionsLabel ${combined_labels[@]} \
    --samplesLabel "${sample_name}" \
    --colorMap viridis \
    --legendLocation upper-right \
    --heatmapWidth 10 \
    --heatmapHeight 20 \
    --xAxisLabel "Distance from TSS (bp)" \
    --refPointLabel "TSS" \
    --perGroup

# Plot combined heatmap (PNG)
plotHeatmap \
    --matrixFile "${HEATMAP_DIR}/combined_genes_matrix.gz" \
    --outFileName "${HEATMAP_DIR}/combined_genes_heatmap.png" \
    --regionsLabel ${combined_labels[@]} \
    --samplesLabel "${sample_name}" \
    --colorMap viridis \
    --legendLocation upper-right \
    --heatmapWidth 10 \
    --heatmapHeight 20 \
    --xAxisLabel "Distance from TSS (bp)" \
    --refPointLabel "TSS" \
    --perGroup

# Profile plot
plotProfile \
    --matrixFile "${HEATMAP_DIR}/combined_genes_matrix.gz" \
    --outFileName "${HEATMAP_DIR}/combined_genes_profile.svg" \
    --regionsLabel ${combined_labels[@]} \
    --samplesLabel "${sample_name}" \
    --xAxisLabel "Distance from TSS (bp)" \
    --refPointLabel "TSS" \
    --plotHeight 7 \
    --plotWidth 12 \
    --perGroup

EOT
