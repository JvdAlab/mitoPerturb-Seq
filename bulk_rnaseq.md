Run Trimmomatic (v0.39) in Paired End mode, with Read 1 and Read 2 fastqs as inputs, outputting trimmed paired and unpaired fastqs for each read. \
ILLUMINACLIP trimming step performed using the TruSeq3-PE-2.fa available here https://github.com/usadellab/Trimmomatic/tree/main/adapters

```
trimmomatic PE \
-threads 24 \
Sample_Read_1.fastq \
Sample_Read_2.fastq \
Sample_Read_1_trimmed_paired.fastq \
Sample_Read_1_trimmed_unpaired.fastq \
Sample_Read_2_trimmed_paired.fastq \
Sample_Read_2_trimmed_unpaired.fastq \
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True \
LEADING:3 \
TRAILING:3 \
MINLEN:36
```

Run RUM (v2.0.4) in Align mode with trimmed Read 1 and Read 2 fastqs as inputs

```
rum_runner align \
--index-dir ../rum_indexes/mm10 \
--output Sample_trimmed_aligned \
--name Sample_trimmed \
--chunks 32 \
--variable-length-reads \
Sample_Read_1_trimmed_paired.fastq Sample_Read_2_trimmed_paired.fastq
```

Run Htseq-count in Union mode with aligned RUM.bam as imput. Output .sam file of alignments is optional

```
htseq-count \
-f bam \
-r name \
-m union \
--additional-attr gene_name \
-o Sample_counts.sam \
-s no \
Sample_trimmed_aligned/RUM.bam \
Reference_genes.gtf > Sample_counts.txt
```
