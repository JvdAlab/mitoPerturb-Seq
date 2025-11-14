#!/bin/bash

################################################################################
# DamID-seq Analysis Pipeline Configuration
#
# This file contains only shared paths and infrastructure settings.
# Script-specific parameters are defined within each individual script.
#
# INSTRUCTIONS:
# 1. Update all paths marked with <PLACEHOLDER>
# 2. Source this file in your scripts: source ../config/config.sh
################################################################################

# ==============================================================================
# DIRECTORY PATHS
# ==============================================================================

# Base directory for analysis
export BASE_DIR="path-to-analysis-directory"

# Working directory containing DamID-seq data subdirectories
export WORK_DIR="${BASE_DIR}/damidseq_wd"

# Reference genome directory
export REF_GENOME_DIR="${BASE_DIR}/ref_genomes"

# Directory for IDR analysis
export IDR_DIR="${WORK_DIR}/IDR_analysis"

# ==============================================================================
# REFERENCE GENOME SETTINGS
# ==============================================================================

# Reference genome
export REF_GENOME="mm10.fa.gz"

# GATC fragment annotation file
export GATC_FRAG_FILE="${REF_GENOME_DIR}/${REF_GENOME_NAME}.GATC.gff"

# Bowtie2 index directory
export BOWTIE2_INDEX_DIR="${REF_GENOME_DIR}/${REF_GENOME_NAME}"

# Genome size (for MACS3)
export GENOME_SIZE=2725537669

# ==============================================================================
# SOFTWARE PATHS
# ==============================================================================

# DamID-seq pipeline path
export DAMIDSEQ_PIPELINE="/home/<USERNAME>/bin/damidseq_pipeline-1.5.3"
