#!/bin/sh
#SBATCH -p cluster7
#SBATCH --ntasks=12
#SBATCH --mem=64gb
#SBATCH --time=05:00:00
#SBATCH -J umi_tools
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL,BEGIN,END

#SBATCH -o /../../slurm-%j.out
#SBATCH -e /../../slurm-%j.err

export OMP_NUM_THREADS=8 

# load conda environment
module load anaconda/3

source activate /../../micromamba/envs/umi-tools

# set correct path to UMI tools
umi_tools extract -I /../../Sample1_R1_001.fastq.gz \
				--bc-pattern=NNNNNNNNNNNNNNNN \
				--log=/../../UMI_extract_cellbarcode.log \
				-S /../../Sample1_R1_001_UMI_extract_NNNNNNNNNNNNNNNN.fastq.gz \
				--read2-in=/../../Sample1_R2_001.fastq.gz \
				--read2-out=/../../Sample1_R2_001_UMI_extract_NNNNNNNNNNNNNNNN.fastq.gz
