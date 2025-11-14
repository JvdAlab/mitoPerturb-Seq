#!/bin/bash

# IDR Analysis - Round 2: Biological Replicates (All Pairwise Combinations)

# Source config file for environment variables
source "../config/config.sh"

# Script-specific parameters
IDR_THRESHOLD_ROUND1=0.01
IDR_THRESHOLD_ROUND2=0.01
IDR_RANK="signal.value"

# Get all Round 1 IDR peak files
sample_peak_files=$(find "${IDR_DIR}" -name "*_idr1_${IDR_THRESHOLD_ROUND1}.broadPeak" ! -name "*_sorted.broadPeak")

# Sort Round 1 IDR files by coordinates
for file in ${sample_peak_files[@]}; do
    sorted_file="${file%.broadPeak}_sorted.broadPeak"
    awk '{$1=$1}1' OFS='\t' "$file" | sort -t$'\t' -k1,1 -k2,2n > "$sorted_file"
done

# Get sorted files
sorted_peak_files=$(find "${IDR_DIR}" -name "*_idr1_${IDR_THRESHOLD_ROUND1}_sorted.broadPeak")

# Convert to array
sample_peak_files_array=($sorted_peak_files)

# Function to extract the sample identifier from filename
extract_sample_id() {
    basename "$1" | cut -d'_' -f1
}

# Generate all pairwise combinations for Round 2 IDR
for ((i=0; i<${#sample_peak_files_array[@]}; i++)); do
    for ((j=i+1; j<${#sample_peak_files_array[@]}; j++)); do
        file1=${sample_peak_files_array[$i]}
        file2=${sample_peak_files_array[$j]}

        sample1=$(extract_sample_id "$file1")
        sample2=$(extract_sample_id "$file2")

        # Create consistent output naming
        if [[ "$sample1" < "$sample2" ]]; then
            output_prefix="${sample1}_${sample2}"
        else
            output_prefix="${sample2}_${sample1}"
        fi

        echo "Submitting Round 2 IDR for ${sample1} and ${sample2}"

        sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=idr2_${sample1}_${sample2}
#SBATCH --output=${IDR_DIR}/${output_prefix}_idr2_${IDR_THRESHOLD_ROUND2}_%j.out
#SBATCH --error=${IDR_DIR}/${output_prefix}_idr2_${IDR_THRESHOLD_ROUND2}_%j.err
#SBATCH --partition=<PARTITION>
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

idr --input-file-type broadPeak \
    --rank ${IDR_RANK} \
    --samples "${file1}" "${file2}" \
    --output-file "${IDR_DIR}/${output_prefix}_idr2_${IDR_THRESHOLD_ROUND2}.broadPeak" \
    --output-file-type broadPeak \
    --plot \
    --log-output-file "${IDR_DIR}/${output_prefix}_idr2_${IDR_THRESHOLD_ROUND2}.log"
EOT

    done
done
