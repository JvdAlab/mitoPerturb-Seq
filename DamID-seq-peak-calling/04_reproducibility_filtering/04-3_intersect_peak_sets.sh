#!/bin/bash

# Multi-Intersection Analysis and Filtering of IDR Peaks

# Source config file for environment variables
source "../config/config.sh"

# Script-specific parameters
IDR_THRESHOLD_ROUND2=0.01
MIN_PEAK_SUPPORT=2

# Find all Round 2 IDR peak files
idr_outfiles=$(find "${IDR_DIR}" -type f -name "*_idr2_${IDR_THRESHOLD_ROUND2}.broadPeak")

# Sort the IDR output files
sorted_idr_outfiles=()

while IFS= read -r file; do
    sorted_file="${file%.broadPeak}_sorted.broadPeak"
    # Ensure tab-delimited format and sort by the first two columns
    awk '{$1=$1}1' OFS='\t' "$file" | sort -t$'\t' -k1,1 -k2,2n > "$sorted_file"
    sorted_idr_outfiles+=("$sorted_file")
done <<< "$idr_outfiles"

# Extract prefixes for multi-intersection
prefixes=($(for file in ${sorted_idr_outfiles[@]}; do basename "$file" | cut -d'_' -f1; done))

# Perform multi-intersection
output_multiinter="${IDR_DIR}/ATF4_idr2_${IDR_THRESHOLD_ROUND2}_peaks_multiinter.bed"

bedtools multiinter \
    -header \
    -names "${prefixes[@]}" \
    -i "${sorted_idr_outfiles[@]}" > "${output_multiinter}"

# Filter multi-intersection results
OUTPUT_DIR="${IDR_DIR}/filtered_peaks"
mkdir -p "${OUTPUT_DIR}"

OVERLAP_ALL="${OUTPUT_DIR}/ATF4_peaks_all_replicates.bed"
OVERLAP_MIN="${OUTPUT_DIR}/ATF4_peaks_min${MIN_PEAK_SUPPORT}_reps_overlap.bed"
OVERLAP_MIN_MERGED="${OUTPUT_DIR}/ATF4_peaks_min${MIN_PEAK_SUPPORT}_reps_overlap_merged.bed"

# Filter for peaks in all comparisons
num_replicates=$(head -n 1 "${output_multiinter}" | awk '{print NF - 4}')
awk -v n="${num_replicates}" 'BEGIN {OFS="\t"} NR > 1 && $4 == n {print $1, $2, $3}' \
    "${output_multiinter}" > "${OVERLAP_ALL}"

# Filter for peaks in at least MIN_PEAK_SUPPORT comparisons
awk -v min="${MIN_PEAK_SUPPORT}" 'BEGIN {OFS="\t"} NR > 1 && $4 >= min {print $1, $2, $3}' \
    "${output_multiinter}" > "${OVERLAP_MIN}"

# Merge overlapping peaks
bedtools merge -i "${OVERLAP_MIN}" > "${OVERLAP_MIN_MERGED}"
