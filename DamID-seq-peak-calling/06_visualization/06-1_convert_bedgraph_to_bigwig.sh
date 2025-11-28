#!/bin/bash

# Process bedGraph Files for Visualization

source "../config/config.sh"

sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=bedgraph_to_bw
#SBATCH --output=${WORK_DIR}/bedgraph_to_bw_%j.out
#SBATCH --error=${WORK_DIR}/bedgraph_to_bw_%j.err
#SBATCH --partition=<PARTITION>
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Find directories
directories=()
for dir in "${WORK_DIR}"/*/ ; do
    if [[ -d "\$dir" && "\$dir" != *"DAM"* ]]; then
        directories+=("\${dir%/}")
    fi
done

# Find bedGraph files
bedgraph_files=()
for dir in "\${directories[@]}"; do
    for file in "\${dir}"/*.gatc.bedgraph; do
        if [ -f "\$file" ]; then
            bedgraph_files+=("\$file")
        fi
    done
done

# Un-log bedGraph files
for file in "\${bedgraph_files[@]}"; do
    awk -v OFS='\t' '{print \$1, \$2, \$3, 2^\$4}' "\$file" > "\${file%.bedgraph}_unlog.bedgraph"
done

# Generate chromosome sizes file
if [[ ! -f "\${REF_GENOME_DIR}/\${REF_GENOME_FASTA}.chrom.sizes" ]]; then
    samtools faidx "\${REF_GENOME_DIR}/\${REF_GENOME_FASTA}"
    cut -f1,2 "\${REF_GENOME_DIR}/\${REF_GENOME_FASTA}.fai" > "\${REF_GENOME_DIR}/\${REF_GENOME_FASTA}.chrom.sizes"
fi

CHROM_SIZES="\${REF_GENOME_DIR}/\${REF_GENOME_FASTA}.chrom.sizes"

# Convert bedGraph to bigWig
unlog_bw_files=()
for dir in "\${directories[@]}"; do
    for bedgraph_file in "\$dir"/*_unlog.bedgraph; do
        if [ -f "\$bedgraph_file" ]; then
            bw_file="\${bedgraph_file%.bedgraph}.bw"
            sort -k1,1 -k2,2n "\$bedgraph_file" > "\${bedgraph_file}.sorted"
            bedGraphToBigWig "\${bedgraph_file}.sorted" "\${CHROM_SIZES}" "\$bw_file"
            unlog_bw_files+=("\$bw_file")
            rm "\${bedgraph_file}.sorted"
        fi
    done
done

# Average bigWig files across replicates
OUTPUT_BW="\${WORK_DIR}/ATF4_rpm-norm_avg_unlog.bw"

bigwigAverage \
    --bigwigs \${unlog_bw_files[@]} \
    --outFileName "\${OUTPUT_BW}" \
    --outFileFormat bigwig \
    --numberOfProcessors 4

EOT
