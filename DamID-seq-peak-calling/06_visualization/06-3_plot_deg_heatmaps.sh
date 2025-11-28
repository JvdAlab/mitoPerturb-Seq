#!/bin/bash

# Generate deepTools Heatmaps at DEG Regions

source "../config/config.sh"

# Script-specific parameters
TSS_UPSTREAM=5000
TSS_DOWNSTREAM=5000

# Define input files
BW_FILE="${WORK_DIR}/ATF4_rpm-norm_avg.bw"
DEG_DIR="${BASE_DIR}/deg_genelists"
HEATMAP_DIR="${DEG_DIR}/deeptools_heatmaps"
mkdir -p "${HEATMAP_DIR}"

# Find DEG BED files
DEG_bed_files=($(find "${DEG_DIR}" -type f -name "*_genes.bed"))
OVERLAP_bed_files=($(find "${DEG_DIR}" -type f -name "*_overlap.bed"))

sample_name=$(basename "${BW_FILE}" | cut -d '-' -f1)

# Submit jobs for each DEG set
for bed_file in "${DEG_bed_files[@]}"; do
    bed_name=$(basename "${bed_file}" | cut -d'_' -f1)

    # Find corresponding overlap file
    overlap_bed_file=""
    for overlap_file in "${OVERLAP_bed_files[@]}"; do
        if [[ "${overlap_file}" == *"${bed_name}"* ]]; then
            overlap_bed_file="${overlap_file}"
            break
        fi
    done

    sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=heatmap_${bed_name}
#SBATCH --output=${HEATMAP_DIR}/${bed_name}_heatmap_%j.out
#SBATCH --error=${HEATMAP_DIR}/${bed_name}_heatmap_%j.err
#SBATCH --partition=<PARTITION>
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Prepare regions list
if [[ -n "${overlap_bed_file}" && -f "${overlap_bed_file}" ]]; then
    REGIONS_FILES=("${overlap_bed_file}" "${bed_file}")
    REGION_LABELS=("${bed_name}_ATF4_overlap" "${bed_name}_all_DEGs")
else
    REGIONS_FILES=("${bed_file}")
    REGION_LABELS=("${bed_name}_DEGs")
fi

# Compute signal matrix
computeMatrix reference-point \
    --scoreFileName "${BW_FILE}" \
    --regionsFileName "\${REGIONS_FILES[@]}" \
    --beforeRegionStartLength ${TSS_UPSTREAM} \
    --afterRegionStartLength ${TSS_DOWNSTREAM} \
    --referencePoint TSS \
    --numberOfProcessors 4 \
    --skipZeros \
    --outFileName "${HEATMAP_DIR}/${bed_name}_matrix.gz"

# Plot heatmap (SVG)
plotHeatmap \
    --matrixFile "${HEATMAP_DIR}/${bed_name}_matrix.gz" \
    --outFileName "${HEATMAP_DIR}/${bed_name}_heatmap.svg" \
    --outFileNameMatrix "${HEATMAP_DIR}/${bed_name}_matrix.tab" \
    --colorMap viridis \
    --samplesLabel "${sample_name}" \
    --regionsLabel "\${REGION_LABELS[@]}" \
    --xAxisLabel "Distance from TSS (bp)" \
    --refPointLabel "TSS" \
    --heatmapHeight 15 \
    --heatmapWidth 3

# Plot heatmap (PNG)
plotHeatmap \
    --matrixFile "${HEATMAP_DIR}/${bed_name}_matrix.gz" \
    --outFileName "${HEATMAP_DIR}/${bed_name}_heatmap.png" \
    --colorMap viridis \
    --samplesLabel "${sample_name}" \
    --regionsLabel "\${REGION_LABELS[@]}" \
    --xAxisLabel "Distance from TSS (bp)" \
    --refPointLabel "TSS" \
    --heatmapHeight 15 \
    --heatmapWidth 3

# Profile plot
plotProfile \
    --matrixFile "${HEATMAP_DIR}/${bed_name}_matrix.gz" \
    --outFileName "${HEATMAP_DIR}/${bed_name}_profile.svg" \
    --samplesLabel "${sample_name}" \
    --regionsLabel "\${REGION_LABELS[@]}" \
    --xAxisLabel "Distance from TSS (bp)" \
    --refPointLabel "TSS" \
    --plotHeight 7 \
    --plotWidth 10

EOT

done
