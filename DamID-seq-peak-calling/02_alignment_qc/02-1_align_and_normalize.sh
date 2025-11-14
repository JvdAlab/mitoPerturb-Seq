#!/bin/bash

# Run DamID-seq Pipeline for Alignment and Normalization

# Source config files for environment variables
source "../config/config.sh"

# Parameters for damidseq_pipeline
GATC_BIN_SIZE=300
NORM_METHOD="rpm"

# Script specific variables
GATC_FRAG_FILE="${REF_GENOME_DIR}/${REF_GENOME_NAME}.GATC.gff"
BOWTIE2_INDEX_DIR="${REF_GENOME_DIR}/${REF_GENOME}"

# Get list of sample directories
directories=()
for dir in "${WORK_DIR}"/*/ ; do
    if [ -d "$dir" ]; then
        directories+=("$dir")
    fi
done

# Submit SLURM job for each directory
for dir_index in "${!directories[@]}"; do
    selected_dir=${directories[$dir_index]}

    echo "Submitting job for data directory: $selected_dir"
    sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=damidseq_${dir_index}
#SBATCH --output=${selected_dir}/damidseq_%j.out
#SBATCH --error=${selected_dir}/damidseq_%j.err
#SBATCH --partition=<PARTITION>
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

cd "${selected_dir}"

perl ${DAMIDSEQ_PIPELINE} \
    --paired \
    --gatc_frag_file="${GATC_FRAG_FILE}" \
    --bowtie2_genome_dir="${BOWTIE2_INDEX_DIR}" \
    --bins=${GATC_BIN_SIZE} \
    --norm_method=${NORM_METHOD}

EOT

done
