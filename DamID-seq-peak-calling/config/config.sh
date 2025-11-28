#!/bin/bash

################################################################################
# DamID-seq Analysis Pipeline Configuration
#
# This file defines shared paths, infrastructure settings, and environment variables.
#
# AUTOMATION:
# When sourced, this script will automatically:
# 1. Define all necessary environment variables.
# 2. Create the directory structure (WORK_DIR, REF_GENOME_DIR, IDR_DIR) if missing.
#
# INSTRUCTIONS:
# 1. Update variables in the "USER CONFIGURATION" section.
# 2. Source this file in every script: source ../config/config.sh
################################################################################

# ==============================================================================
# 1. USER CONFIGURATION (EDIT THESE)
# ==============================================================================

# Base directory for the entire analysis
export BASE_DIR="/path/to/your/analysis_folder"

# Filename of your reference genome (Must be placed inside REF_GENOME_DIR)
export REF_GENOME="mm10.fa.gz"

# Path to the DamID-seq pipeline executable
# Ensure this points to the damidseq_pipeline.pl script
export DAMIDSEQ_PIPELINE="/home/<USERNAME>/bin/damidseq_pipeline-1.5.3/damidseq_pipeline.pl"

# Effective Genome size for MACS3 (default is mm10: 2.72e9)
# hs: 2.7e9, mm: 1.87e9, ce: 9e7, dm: 1.2e8
export GENOME_SIZE=2725537669

# ==============================================================================
# 2. DERIVED PATHS (DO NOT EDIT USUALLY)
# ==============================================================================

# Working directory for DamID-seq intermediate files
export WORK_DIR="${BASE_DIR}/damidseq_wd"

# Directory where reference genome and indices are stored
export REF_GENOME_DIR="${BASE_DIR}/ref_genomes"

# Directory for IDR analysis and final results
export IDR_DIR="${WORK_DIR}/IDR_analysis"

# Extract the Genome Name (removes .fa, .fasta, .gz extensions)
# Used to name indexes and fragments automatically
# e.g., "mm10.fa.gz" -> "mm10"
REF_GENOME_NAME=$(basename "${REF_GENOME}")
REF_GENOME_NAME=${REF_GENOME_NAME%.gz}
REF_GENOME_NAME=${REF_GENOME_NAME%.fa}
REF_GENOME_NAME=${REF_GENOME_NAME%.fasta}
export REF_GENOME_NAME

# Full path to GATC fragment file
export GATC_FRAG_FILE="${REF_GENOME_DIR}/${REF_GENOME_NAME}.GATC.gff"

# Full path to Bowtie2 Index base
export BOWTIE2_INDEX_DIR="${REF_GENOME_DIR}/${REF_GENOME_NAME}"

# ==============================================================================
# 3. DIRECTORY INITIALIZATION
# ==============================================================================

# Automatically create directories if they don't exist
# This prevents "directory not found" errors in downstream scripts
if [ ! -d "${WORK_DIR}" ]; then
    echo "[Config] Creating working directory: ${WORK_DIR}"
    mkdir -p "${WORK_DIR}"
fi

if [ ! -d "${REF_GENOME_DIR}" ]; then
    echo "[Config] Creating reference directory: ${REF_GENOME_DIR}"
    mkdir -p "${REF_GENOME_DIR}"
fi

if [ ! -d "${IDR_DIR}" ]; then
    echo "[Config] Creating IDR analysis directory: ${IDR_DIR}"
    mkdir -p "${IDR_DIR}"
fi

# ==============================================================================
# 4. CONFIGURATION SUMMARY
# ==============================================================================
# Optional: Uncomment for debugging
# echo "Configuration loaded. Root: $BASE_DIR | Genome: $REF_GENOME_NAME"