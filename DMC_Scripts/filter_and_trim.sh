#3_filtering_and_trimming.sh
#!/bin/bash

############################################
#     Bioinformatic pipeline               #
#     ITS amplicon sequences               #
#               STATS                      #
############################################
module load vsearch
# STEP 3 - FILTERING AND TRIMMING
vsearch -fastq_filter trimmed_combined.fastq -fastq_maxee 1.0 -fastq_trunclen 263 -fastq_maxns 0 -fastaout filtered_R1.fasta -fastqout filtered_R1.fastq
vsearch -fastq_filter filtered_R1.fastq -fastq_stripleft 44 -fastaout trimmed_R1.fasta -fastqout trimmed_R1.fastq
module load fastqc
fastqc trimmed_R1.fastq
