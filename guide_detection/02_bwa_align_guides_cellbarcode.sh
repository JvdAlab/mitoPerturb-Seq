#!/bin/sh
#SBATCH -p cluster7
#SBATCH --ntasks=12
#SBATCH --mem=64gb
#SBATCH --time=05:00:00
#SBATCH -J bwa_align
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL,BEGIN,END

#SBATCH -o /../../slurm-%j.out
#SBATCH -e /../../slurm-%j.err

export OMP_NUM_THREADS=8

# load the necessary modules
module load anaconda/3
module load samtools

# choose appropriate reference file
GENOME_FASTA="/../../CROPseq_85bp_reference.fa"

# set home/data directories
HOME_DIR="/../../"
DATA_DIR="/../../"


: '
STEP 1 - ALIGNMENT
'

for SAMPLE_SRR in ${DATA_DIR}/Sample1_R2_001_UMI_extract_NNNNNNNNNNNNNNNN.fastq.gz ; do    

    echo "STEP 1 - commencing"
    echo $SAMPLE_SRR

    bwa mem $GENOME_FASTA ${SAMPLE_SRR} | samtools sort -o ${SAMPLE_SRR/R2_001_UMI_extract_NNNNNNNNNNNNNNNN.fastq.gz/UMI_extract_NNNNNNNNNNNNNNNN_bwa_85bp_R2_aligned_guides.bam}

    samtools view -F 4  ${SAMPLE_SRR/R2_001_UMI_extract_NNNNNNNNNNNNNNNN.fastq.gz/UMI_extract_NNNNNNNNNNNNNNNN_bwa_85bp_R2_aligned_guides.bam} >  ${SAMPLE_SRR/R2_001_UMI_extract_NNNNNNNNNNNNNNNN.fastq.gz/UMI_extract_NNNNNNNNNNNNNNNN_bwa_85bp_R2_aligned_guides_F4.sam}

    echo "Done."
    echo "STEP 1 - complete"

: '
done
'


done

