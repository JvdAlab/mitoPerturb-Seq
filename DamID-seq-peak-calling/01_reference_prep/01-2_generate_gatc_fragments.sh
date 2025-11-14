#!/bin/bash

# Generate GATC fragment annotation

# Source config file for environment variables
source "../config/config.sh"

perl gatc.track.maker.pl \
    --name="${REF_GENOME_DIR}/${REF_GENOME%.fa.gz}" \
    "${REF_GENOME_DIR}/${REF_GENOME}"
