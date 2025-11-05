```
trimmomatic PE \
-threads 24 \
Read_1.fastq \
Read_2.fastq \
Read_1_trimmed_paired.fastq \
Read_1_trimmed_unpaired.fastq \
Read_2_trimmed_paired.fastq \
Read_2_trimmed_unpaired.fastq \
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True \
LEADING:3 \
TRAILING:3 \
MINLEN:36
```
