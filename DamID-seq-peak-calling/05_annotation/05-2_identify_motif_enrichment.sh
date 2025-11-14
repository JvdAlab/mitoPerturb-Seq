#!/bin/bash

# Motif Enrichment Analysis using HOMER

source "../config/config.sh"

# Script-specific parameters
MIN_PEAK_SUPPORT=2
MOTIF_SIZE=200

# Define input and output
PEAK_FILE="${IDR_DIR}/filtered_peaks/ATF4_peaks_min${MIN_PEAK_SUPPORT}_reps_overlap_merged.bed"
OUTPUT_DIR="${WORK_DIR}/motif_analysis"
mkdir -p "${OUTPUT_DIR}"

# Run HOMER findMotifsGenome.pl
findMotifsGenome.pl \
    "${PEAK_FILE}" \
    ${REF_GENOME_NAME} \
    "${OUTPUT_DIR}" \
    -size ${MOTIF_SIZE} \
    -mask
