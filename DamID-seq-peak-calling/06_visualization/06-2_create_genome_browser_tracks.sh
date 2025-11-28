#!/bin/bash

# Generate Genome Browser Tracks with pyGenomeTracks

source "../config/config.sh"

# Script-specific parameters
MIN_PEAK_SUPPORT=2

# Define input and output
BW_FILE="${WORK_DIR}/ATF4_rpm-norm_avg_unlog.bw"
OUTPUT_DIR="${WORK_DIR}/genome_tracks"
mkdir -p "${OUTPUT_DIR}"
TRACKS_INI="${OUTPUT_DIR}/ATF4_tracks.ini"

# Create template tracks.ini if doesn't exist
if [[ ! -f "${TRACKS_INI}" ]]; then
    cat > "${TRACKS_INI}" <<'INIFILE'
[x-axis]
where = top

[ATF4 Signal]
file = BIGWIG_FILE_PLACEHOLDER
height = 4
title = ATF4 DamID Signal
color = darkblue
min_value = 0
max_value = auto
file_type = bigwig

[ATF4 Peaks]
file = PEAKS_FILE_PLACEHOLDER
height = 1
title = ATF4 Peaks
color = red
file_type = bed
display = collapsed

[genes]
file = mm10
height = 3
title = Genes
file_type = gtf
gene_rows = 10
INIFILE

    sed -i "s|BIGWIG_FILE_PLACEHOLDER|${BW_FILE}|g" "${TRACKS_INI}"
    sed -i "s|PEAKS_FILE_PLACEHOLDER|${IDR_DIR}/filtered_peaks/ATF4_peaks_min${MIN_PEAK_SUPPORT}_support_merged.bed|g" "${TRACKS_INI}"
fi

# Define regions of interest
declare -a REGIONS=(
    "chr8:27150000-27320000:Eif4ebp1"
    "chr18:42202150-42299457:Lars"
    "chr8:111006992-111076403:Aars"
    "chr4:3790187-3928053:Rps20"
    "chr7:44792519-44836396:Atf5"
    "chr6:7668013-7713186:Asns"
)

sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=pygenometracks
#SBATCH --output=${OUTPUT_DIR}/pygenometracks_%j.out
#SBATCH --error=${OUTPUT_DIR}/pygenometracks_%j.err
#SBATCH --partition=<PARTITION>
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Process each region
for region_info in "${REGIONS[@]}"; do
    IFS=':' read -ra PARTS <<< "\$region_info"
    chr="\${PARTS[0]}"
    coords="\${PARTS[1]}"
    gene="\${PARTS[2]}"
    region="\${chr}:\${coords}"

    pyGenomeTracks \
        --tracks "${TRACKS_INI}" \
        --region "\${region}" \
        --dpi 320 \
        --out "${OUTPUT_DIR}/ATF4_tracks_\${gene}.svg"

    pyGenomeTracks \
        --tracks "${TRACKS_INI}" \
        --region "\${region}" \
        --dpi 320 \
        --out "${OUTPUT_DIR}/ATF4_tracks_\${gene}.png"
done

EOT
