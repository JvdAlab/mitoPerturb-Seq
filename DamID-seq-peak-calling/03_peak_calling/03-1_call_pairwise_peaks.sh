#!/bin/bash

# Call Peaks with MACS3 (Pairwise Comparisons)

# Source config file for environment variables
source "../config/config.sh"

# MACS3 parameters
MACS3_QVALUE=0.05
MACS3_EXTSIZE=300
MACS3_SHIFT=0
MACS3_BW=300
MACS3_MFOLD_LOW=5
MACS3_MFOLD_HIGH=50

# List the directories in the work directory
directories=()
for dir in "${WORK_DIR}"/*; do
    if [ -d "${dir}" ]; then
        directories+=("${dir}")
    fi
done

# Iterate over each directory
for dir_index in "${!directories[@]}"; do
    dir="${directories[$dir_index]}"
    echo "Processing directory: ${dir}"

    dam_only_file=""
    dam_fusion_file=""

    # Identify Dam-only and Dam-fusion BAM files
    for bam_file in "${dir}"/*.bam; do
        if [[ $(basename "${bam_file}") == *"DAM"* ]]; then
            dam_only_file="${bam_file}"
        else
            dam_fusion_file="${bam_file}"
        fi
    done

    if [ -z "${dam_only_file}" ] || [ -z "${dam_fusion_file}" ]; then
        echo "Error: Could not find DAM-only or DAM-fusion bam file in directory: ${dir}"
        continue
    else
    # Get the prefix for the peak file
        f_peak_prefix="$(basename "${dam_fusion_file}" .bam)_vs_$(basename "${dam_only_file}" .bam)"
    fi

    echo "Submitting job for directory: ${dir}"

    # Submit SLURM job for MACS3
    sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=macs3_pairwise_${dir_index}
#SBATCH --output=${dir}/macs3_pairwise_%j.out
#SBATCH --error=${dir}/macs3_pairwise_%j.err
#SBATCH --partition=<PARTITION>
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

macs3 callpeak \
    --broad \
    --format BAMPE \
    --treatment "${dam_fusion_file}" \
    --control "${dam_only_file}" \
    --name "${f_peak_prefix}" \
    --outdir "${dir}" \
    --gsize ${GENOME_SIZE} \
    --keep-dup all \
    --nomodel \
    --extsize ${MACS3_EXTSIZE} \
    --shift ${MACS3_SHIFT} \
    --bw ${MACS3_BW} \
    --qvalue ${MACS3_QVALUE} \
    --mfold ${MACS3_MFOLD_LOW} ${MACS3_MFOLD_HIGH} &>> ${dir}/${f_peak_prefix}.macs3.log

EOT

done
