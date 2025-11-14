#!/bin/bash

# IDR Analysis - Round 1: Technical Replicates

# Source config file for environment variables
source "../config/config.sh"

# Script-specific parameters
IDR_THRESHOLD_ROUND1=0.01
IDR_RANK="signal.value"

mkdir -p "${IDR_DIR}"

# Get all broadPeak files
peak_files=$(find "${WORK_DIR}" -name "*.broadPeak" ! -name "*_sorted.broadPeak" ! -name "*_idr_*")

# Sort all peak files by coordinates
for file in ${peak_files[@]}; do
    sorted_file="${file%.broadPeak}_sorted.broadPeak"
    awk '{$1=$1}1' OFS='\t' "$file" | sort -t$'\t' -k1,1 -k2,2n > "$sorted_file"
done

sorted_peak_files=$(find "${WORK_DIR}" -name "*_sorted.broadPeak")

# Group sorted files by Dam-fusion replicate prefix
declare -A sorted_peak_groups

for file in ${sorted_peak_files[@]}; do
    prefix=$(basename "$file" | cut -d'_' -f1)
    sorted_peak_groups[$prefix]+="$file "
done

# Run IDR on each group of sorted files
for prefix in "${!sorted_peak_groups[@]}"; do
    files=(${sorted_peak_groups[$prefix]})

    if [[ ${#files[@]} -eq 2 ]]; then
        echo "Running IDR on files $(basename "${files[0]}") and $(basename "${files[1]}")"
        sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=idr_${prefix}
#SBATCH --output=${IDR_DIR}/${prefix}_idr1_${IDR_THRESHOLD_ROUND1}_%j.out
#SBATCH --error=${IDR_DIR}/${prefix}_idr1_${IDR_THRESHOLD_ROUND1}_%j.err
#SBATCH --partition=<PARTITION>
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

idr --input-file-type broadPeak \
    --rank ${IDR_RANK} \
    --samples "${files[0]}" "${files[1]}" \
    --output-file "${IDR_DIR}/${prefix}_idr1_${IDR_THRESHOLD_ROUND1}.broadPeak" \
    --output-file-type broadPeak \
    --plot \
    --log-output-file "${IDR_DIR}/${prefix}_idr1_${IDR_THRESHOLD_ROUND1}.log"
EOT
    else
        echo "Error: Expected 2 files for prefix ${prefix}, but found ${#files[@]}"
    fi
done
