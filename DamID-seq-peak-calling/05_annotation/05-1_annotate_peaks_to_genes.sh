#!/bin/bash

# Annotate Peaks to Nearest Genes using HOMER

source "../config/config.sh"

# Script-specific parameters
MIN_PEAK_SUPPORT=2

# Define input and output
PEAK_FILE="${IDR_DIR}/filtered_peaks/ATF4_peaks_min${MIN_PEAK_SUPPORT}_reps_overlap_merged.bed"
OUTPUT_DIR="${IDR_DIR}/annotation"
mkdir -p "${OUTPUT_DIR}"
ANNOTATION_FILE="${OUTPUT_DIR}/ATF4_peaks_annotated.txt"

# Run HOMER annotatePeaks.pl
annotatePeaks.pl \
    "${PEAK_FILE}" \
    ${REF_GENOME_NAME} \
    > "${ANNOTATION_FILE}"
