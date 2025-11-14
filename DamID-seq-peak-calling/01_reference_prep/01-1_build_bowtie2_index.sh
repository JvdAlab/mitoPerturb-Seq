#!/bin/bash

# Build Bowtie2 index for reference genome

# Source config file for environment variables
source "../config/config.sh"

# Uncompress genome if needed
if [[ $REF_GENOME == *.gz ]]; then
    echo "Decompressing the reference genome file..."
    gunzip -c "${REF_GENOME_DIR}/${REF_GENOME}" > "${REF_GENOME_DIR}/${REF_GENOME%.gz}"

    # Update variables
    REF_GENOME=${REF_GENOME%.gz}
    REF_GENOME_INDEX=${REF_GENOME%.fa}

# Build Bowtie2 index
echo "Creating the Bowtie2 index..."
bowtie2-build \
    "${REF_GENOME_DIR}/${REF_GENOME}" \
    "${REF_GENOME_DIR}/${REF_GENOME_INDEX}"
